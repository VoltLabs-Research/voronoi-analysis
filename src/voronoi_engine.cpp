#include <volt/voronoi_engine.h>
#include <volt/voronoi_cell.h>
#include <spdlog/spdlog.h>
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/blocked_range.h>
#include <algorithm>
#include <cmath>
#include <vector>

namespace Volt{

using Volt::Particles::ParticleProperty;

VoronoiAnalysisEngine::VoronoiAnalysisEngine(
    ParticleProperty* positions,
    const SimulationCell& simCell,
    ParticleProperty* radii,
    ParticleProperty* selection,
    const Params& params
)
    : _positions(positions),
      _simCell(simCell),
      _radii(radii),
      _selection(selection),
      _params(params),
      _effectiveCutoff(params.cutoff)
{
    const size_t n = _positions ? _positions->size() : 0;

    // Use the 5-arg DataType-explicit constructor for UserProperty outputs
    // because the 4-arg constructor pins UserProperty to DataType::Void.
    _atomicVolumes = std::make_shared<ParticleProperty>(
        n, Volt::Particles::DataType::Double, 1, sizeof(double), true);
    _coordNumbers = std::make_shared<ParticleProperty>(
        n, ParticleProperty::CoordinationProperty, 0, true);
    _maxFaceOrder = std::make_shared<ParticleProperty>(
        n, Volt::Particles::DataType::Int, 1, sizeof(int), true);

    const int indexComponents = std::max(1, _params.edgeCount - 2);
    _voronoiIndex = std::make_shared<ParticleProperty>(
        n, Volt::Particles::DataType::Int,
        static_cast<size_t>(indexComponents),
        static_cast<size_t>(indexComponents) * sizeof(int), true);
}

double VoronoiAnalysisEngine::estimateAutoCutoff(){
    // Estimate a cutoff from the mean volume per atom (sphere-equivalent
    // radius) and scale up generously so the bootstrap cube safely contains
    // every Voronoi vertex. For a perfect crystal this is ~2-3x the nearest
    // neighbor distance, which is the typical recommendation for voro++.
    const size_t n = _positions ? _positions->size() : 0;
    if(n < 2) return 3.0;

    const double cellVol = _simCell.volume3D();
    if(!(cellVol > 0.0) || !std::isfinite(cellVol)){
        return 3.0;
    }
    const double volPerAtom = cellVol / static_cast<double>(n);
    // Radius of sphere with volume volPerAtom: r = (3V/4pi)^(1/3)
    const double rSphere = std::cbrt(3.0 * volPerAtom / (4.0 * M_PI));
    // Nearest-neighbor distance in a dense packing is ~2*rSphere; add a 2x
    // safety margin so the bootstrap cube is comfortably larger than the
    // true Voronoi cell.
    const double cutoff = 4.0 * rSphere;
    return std::max(cutoff, 2.5);
}

void VoronoiAnalysisEngine::perform(){
    if(!_positions || _positions->size() == 0) return;

    if(!(_effectiveCutoff > 0.0)){
        _effectiveCutoff = estimateAutoCutoff();
        spdlog::info("Voronoi: auto-selected cutoff = {:.4f}", _effectiveCutoff);
    }

    CutoffNeighborFinder finder;
    if(!finder.prepare(_effectiveCutoff, _positions, _simCell)){
        spdlog::error("Voronoi: failed to build neighbor list");
        return;
    }

    const size_t n = _positions->size();

    oneapi::tbb::parallel_for(
        oneapi::tbb::blocked_range<size_t>(0, n),
        [this, &finder](const oneapi::tbb::blocked_range<size_t>& range){
            // The finder is thread-safe for concurrent Query construction
            // because Query stores all iteration state.
            CutoffNeighborFinder& sharedFinder = const_cast<CutoffNeighborFinder&>(finder);
            for(size_t i = range.begin(); i < range.end(); ++i){
                processAtom(i, sharedFinder);
            }
        }
    );
}

void VoronoiAnalysisEngine::processAtom(size_t atomIndex, CutoffNeighborFinder& finder){
    // Respect onlySelected: atoms outside the selection get zeroed output.
    if(_selection){
        if(_selection->getInt(atomIndex) == 0){
            _atomicVolumes->setDouble(atomIndex, 0.0);
            _coordNumbers->setInt(atomIndex, 0);
            _maxFaceOrder->setInt(atomIndex, 0);
            const size_t comps = _voronoiIndex->componentCount();
            for(size_t c = 0; c < comps; ++c){
                _voronoiIndex->setIntComponent(atomIndex, c, 0);
            }
            return;
        }
    }

    // Collect neighbors as (distance, index, delta) triples. We sort by
    // distance because closer planes tend to carve away more of the cell
    // first, which reduces the per-clip work.
    struct Neigh {
        double distSq;
        Vector3 delta;
        int neighIndex;
    };
    std::vector<Neigh> neighbors;
    neighbors.reserve(64);

    for(CutoffNeighborFinder::Query q(finder, atomIndex); !q.atEnd(); q.next()){
        Neigh nb;
        nb.distSq = q.distanceSquared();
        nb.delta = q.delta();
        nb.neighIndex = static_cast<int>(q.current());
        neighbors.push_back(nb);
    }

    std::sort(neighbors.begin(), neighbors.end(),
        [](const Neigh& a, const Neigh& b){ return a.distSq < b.distSq; });

    // Seed the bootstrap cube. Half-side = effectiveCutoff guarantees the
    // initial cube encloses every neighbor-bisector plane we are about to
    // apply.
    VoronoiCell cell;
    cell.initCube(_effectiveCutoff);

    const double rI = (_params.useRadii && _radii) ? _radii->getDouble(atomIndex) : 0.0;
    const double rISq = rI * rI;

    for(const Neigh& nb : neighbors){
        if(cell.isEmpty()) break;

        const Vector3& d = nb.delta;
        const double dSq = nb.distSq;
        if(dSq <= 0.0) continue;

        // Plane normal in the atom-local frame is the neighbor direction.
        // Placement of the plane along that direction:
        //   - plain Voronoi: midpoint => plane at d/2 along |d|
        //   - radical (power) Voronoi: offset by the radii difference.
        double t = 0.5;
        if(_params.useRadii && _radii){
            const double rJ = _radii->getDouble(static_cast<size_t>(nb.neighIndex));
            t = 0.5 + (rISq - rJ * rJ) / (2.0 * dSq);
        }

        // Plane equation in atom-local coords: n.x = d (where n = unit(d))
        // Using non-normalized d simplifies to d.x = t * dSq.
        const double planeD = t * dSq;
        cell.clip(d, planeD, nb.neighIndex);
    }

    // Evaluate cell properties.
    const double volume = cell.volume();
    _atomicVolumes->setDouble(atomIndex, volume);

    // Count neighbor-generated faces that survived the thresholds.
    int coordNum = 0;
    int maxOrder = 0;
    const int indexComponents = static_cast<int>(_voronoiIndex->componentCount());
    std::vector<int> indexCounts(static_cast<size_t>(indexComponents), 0);

    for(const VoronoiCell::Face& face : cell.faces()){
        if(face.neighborIndex < 0) continue; // bootstrap face remnant
        const int order = VoronoiCell::faceOrder(face);
        if(order < 3) continue;

        const double area = cell.faceArea(face);
        if(area < _params.faceThreshold) continue;

        const double minEdge = cell.faceMinEdgeLength(face);
        if(minEdge < _params.edgeThreshold) continue;

        ++coordNum;
        if(order > maxOrder) maxOrder = order;

        if(_params.computeIndices){
            // voronoi_index[k] counts faces of order (k+3), up to edgeCount.
            const int bin = order - 3;
            if(bin >= 0 && bin < indexComponents){
                indexCounts[static_cast<size_t>(bin)] += 1;
            }
        }
    }

    _coordNumbers->setInt(atomIndex, coordNum);
    _maxFaceOrder->setInt(atomIndex, maxOrder);

    for(int c = 0; c < indexComponents; ++c){
        _voronoiIndex->setIntComponent(atomIndex,
            static_cast<size_t>(c),
            _params.computeIndices ? indexCounts[static_cast<size_t>(c)] : 0);
    }
}

}
