#include <volt/voronoi_engine.h>
#include <volt/voronoi_cell.h>
#include <spdlog/spdlog.h>
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/blocked_range.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <atomic>

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

    _atomicVolumes = std::make_shared<ParticleProperty>(
        n, Volt::Particles::DataType::Double, 1, sizeof(double), true);
    _coordNumbers = std::make_shared<ParticleProperty>(
        n, ParticleProperty::CoordinationProperty, 0, true);
    _cavityRadii = std::make_shared<ParticleProperty>(
        n, Volt::Particles::DataType::Double, 1, sizeof(double), true);
    _faceIndices.resize(n);
}

double VoronoiAnalysisEngine::estimateAutoCutoff(){
    const size_t n = _positions ? _positions->size() : 0;
    if(n < 2) return 3.0;

    const double cellVol = _simCell.volume3D();
    if(!(cellVol > 0.0) || !std::isfinite(cellVol)){
        return 3.0;
    }
    const double volPerAtom = cellVol / static_cast<double>(n);
    const double rSphere = std::cbrt(3.0 * volPerAtom / (4.0 * PI));
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
            CutoffNeighborFinder& sharedFinder = const_cast<CutoffNeighborFinder&>(finder);
            for(size_t i = range.begin(); i < range.end(); ++i){
                processAtom(i, sharedFinder);
            }
        }
    );

    computeAggregates();
}

void VoronoiAnalysisEngine::processAtom(size_t atomIndex, CutoffNeighborFinder& finder){
    if(_selection){
        if(_selection->getInt(atomIndex) == 0){
            _atomicVolumes->setDouble(atomIndex, 0.0);
            _coordNumbers->setInt(atomIndex, 0);
            _cavityRadii->setDouble(atomIndex, 0.0);
            return;
        }
    }

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

    VoronoiCell cell;
    cell.initCube(_effectiveCutoff);

    const double rI = (_params.useRadii && _radii) ? _radii->getDouble(atomIndex) : 0.0;
    const double rISq = rI * rI;

    for(const Neigh& nb : neighbors){
        if(cell.isEmpty()) break;

        const Vector3& d = nb.delta;
        const double dSq = nb.distSq;
        if(dSq <= 0.0) continue;

        double t = 0.5;
        if(_params.useRadii && _radii){
            const double rJ = _radii->getDouble(static_cast<size_t>(nb.neighIndex));
            t = 0.5 + (rISq - rJ * rJ) / (2.0 * dSq);
        }

        const double planeD = t * dSq;
        cell.clip(d, planeD, nb.neighIndex);
    }

    const double volume = cell.volume();
    _atomicVolumes->setDouble(atomIndex, volume);
    _cavityRadii->setDouble(atomIndex, cell.cavityRadius());

    int coordNum = 0;
    std::vector<int> localFaceOrders;

    for(const VoronoiCell::Face& face : cell.faces()){
        if(face.neighborIndex < 0) continue;
        const int order = VoronoiCell::faceOrder(face);
        if(order < 3) continue;

        const double area = cell.faceArea(face);
        if(area < _params.faceThreshold) continue;

        const double minEdge = cell.faceMinEdgeLength(face);
        if(minEdge < _params.edgeThreshold) continue;

        ++coordNum;
        localFaceOrders.push_back(order);
    }

    std::sort(localFaceOrders.begin(), localFaceOrders.end());
    _faceIndices[atomIndex] = std::move(localFaceOrders);
    _coordNumbers->setInt(atomIndex, coordNum);
}

void VoronoiAnalysisEngine::computeAggregates(){
    const size_t n = _faceIndices.size();
    long long totalFaceCount = 0;
    long long totalFaceOrderSum = 0;
    int maxOrder = 0;
    int validCount = 0;

    for(size_t i = 0; i < n; ++i){
        const auto& fi = _faceIndices[i];
        if(fi.empty()) continue;
        ++validCount;
        for(int o : fi){
            totalFaceOrderSum += o;
            ++totalFaceCount;
            if(o > maxOrder) maxOrder = o;
        }
    }

    _meanFaceOrder = totalFaceCount > 0
        ? static_cast<double>(totalFaceOrderSum) / static_cast<double>(totalFaceCount)
        : 0.0;
    _maxFaceOrderObserved = maxOrder;
    _polyhedraMeshCount = validCount;
}

}
