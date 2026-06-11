#include <volt/voronoi_service.h>
#include <volt/voronoi_engine.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/analysis_result.h>
#include <volt/utilities/parquet_atom_writer.h>
#include <volt/utilities/json_utils.h>
#include <spdlog/spdlog.h>
#include <algorithm>
#include <string>
#include <vector>

namespace Volt{

using namespace Volt::Particles;

VoronoiService::VoronoiService()
    : _cutoff(0.0),
      _edgeThreshold(0.0),
      _faceThreshold(0.0),
      _useRadii(false),
      _onlySelected(false) {}

void VoronoiService::setCutoff(double cutoff){ _cutoff = cutoff; }
void VoronoiService::setEdgeThreshold(double v){ _edgeThreshold = v; }
void VoronoiService::setFaceThreshold(double v){ _faceThreshold = v; }
void VoronoiService::setUseRadii(bool v){ _useRadii = v; }
void VoronoiService::setOnlySelected(bool v){ _onlySelected = v; }

namespace {

// Resolve an optional scalar double column from the LAMMPS frame.
std::shared_ptr<ParticleProperty> tryBindDoubleColumn(
    const LammpsParser::Frame& frame,
    std::initializer_list<const char*> names
){
    for(const char* name : names){
        const auto* col = frame.findAtomProperty(name);
        if(col && col->dataType == DataType::Double && !col->doubles.empty()){
            auto prop = std::make_shared<ParticleProperty>();
            prop->bindExternalData(
                const_cast<void*>(col->constData()),
                col->doubles.size(),
                DataType::Double,
                1,
                sizeof(double),
                nullptr
            );
            return prop;
        }
    }
    return nullptr;
}

std::shared_ptr<ParticleProperty> tryBindIntSelection(const LammpsParser::Frame& frame){
    // Accept either "Selection" or "selection" columns, int-typed.
    for(const char* name : {"Selection", "selection"}){
        const auto* col = frame.findAtomProperty(name);
        if(col && col->dataType == DataType::Int && !col->ints.empty()){
            auto prop = std::make_shared<ParticleProperty>();
            prop->bindExternalData(
                const_cast<void*>(col->constData()),
                col->ints.size(),
                DataType::Int,
                1,
                sizeof(int),
                nullptr
            );
            return prop;
        }
    }
    return nullptr;
}

}

json VoronoiService::compute(const LammpsParser::Frame& frame, const std::string& outputBase){
    auto startTime = std::chrono::high_resolution_clock::now();

    if(frame.natoms <= 0)
        return AnalysisResult::failure("Invalid number of atoms");

    if(!FrameAdapter::validateSimulationCell(frame.simulationCell))
        return AnalysisResult::failure("Invalid simulation cell");

    auto positions = FrameAdapter::createPositionPropertyShared(frame);
    if(!positions)
        return AnalysisResult::failure("Failed to create position property");

    std::shared_ptr<ParticleProperty> radii;
    if(_useRadii){
        radii = tryBindDoubleColumn(frame, {"Radius", "radius", "c_radius"});
        if(!radii){
            spdlog::warn("Voronoi: --useRadii requested but no Radius column found; falling back to plain Voronoi");
        }
    }

    std::shared_ptr<ParticleProperty> selection;
    if(_onlySelected){
        selection = tryBindIntSelection(frame);
        if(!selection){
            spdlog::warn("Voronoi: --onlySelected requested but no Selection column found; treating all atoms as selected");
        }
    }

    VoronoiAnalysisEngine::Params params;
    params.cutoff = _cutoff;
    params.edgeThreshold = _edgeThreshold;
    params.faceThreshold = _faceThreshold;
    params.useRadii = _useRadii && radii != nullptr;

    spdlog::info(
        "Starting Voronoi tessellation (cutoff={}, edgeThr={}, faceThr={}, radii={})",
        params.cutoff, params.edgeThreshold, params.faceThreshold,
        params.useRadii ? "true" : "false"
    );

    VoronoiAnalysisEngine engine(
        positions.get(),
        frame.simulationCell,
        radii.get(),
        selection.get(),
        params
    );
    engine.perform();

    auto volumes = engine.atomicVolumes();
    auto coords = engine.coordNumbers();
    auto cavityRadii = engine.cavityRadii();

    // Aggregate stats.
    double volumeSum = 0.0;
    long long coordSum = 0;
    const int atomCount = frame.natoms;
    for(int i = 0; i < atomCount; ++i){
        volumeSum += volumes->getDouble(i);
        coordSum += coords->getInt(i);
    }
    const double meanVolume = atomCount > 0 ? volumeSum / atomCount : 0.0;
    const double meanCoord = atomCount > 0 ? static_cast<double>(coordSum) / atomCount : 0.0;

    // Summary table (written to <base>_voronoi.parquet) -----------------------
    // Scalar aggregates only; per-atom data now lives in <base>_atoms.parquet.
    json result;
    result["main_listing"] = {
        {"total_atoms", atomCount},
        {"mean_atomic_volume", meanVolume},
        {"mean_coordination", meanCoord},
        {"effective_cutoff", engine.effectiveCutoff()},
        {"edge_threshold", params.edgeThreshold},
        {"face_threshold", params.faceThreshold},
        {"use_radii", params.useRadii},
    };
    AnalysisResult::addTiming(result, startTime);
    result["is_failed"] = false;

    if(!outputBase.empty()){
        const std::string voronoiPath = outputBase + "_voronoi.parquet";
        if(JsonUtils::writeJsonToParquet(result, voronoiPath)){
            spdlog::info("Voronoi summary parquet written to {}", voronoiPath);
        }else{
            spdlog::warn("Could not write Voronoi summary parquet: {}", voronoiPath);
        }

        // Per-atom table (<base>_atoms.parquet) via the canonical streaming writer.
        // streamAtomsToParquet is used directly (instead of serializePluginOutput)
        // so a StructureIdResolver can pin structure_id = coordination, matching the
        // legacy export and the "Coordination_<k>" bucket grouping.
        const std::string atomsPath = outputBase + "_atoms.parquet";
        streamAtomsToParquet(
            atomsPath,
            frame,
            [&coords](std::size_t i){
                return "Coordination_" + std::to_string(coords->getInt(i));
            },
            [&](ColumnarAtomWriter& w, std::size_t i){
                w.field("coordination",  coords->getInt(i));
                w.field("atomic_volume", volumes->getDouble(i));
                w.field("cavity_radius", cavityRadii->getDouble(i));
            },
            [&coords](std::size_t i){
                return coords->getInt(i);
            }
        );
        spdlog::info("Voronoi atoms parquet written to {}", atomsPath);
    }

    return result;
}

}
