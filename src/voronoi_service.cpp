#include <volt/voronoi_service.h>
#include <volt/voronoi_analysis_log.h>
#include <volt/voronoi_engine.h>
#include <volt/core/frame_adapter.h>
#include <volt/core/analysis_result.h>
#include <volt/core/progress.h>
#include <volt/plugin/output_files.h>
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
        radii = FrameAdapter::findColumnShared(frame, {"Radius", "radius", "c_radius"}, DataType::Double);
        if(!radii){
            spdlog::warn("Voronoi: --useRadii requested but no Radius column found; falling back to plain Voronoi");
        }
    }

    std::shared_ptr<ParticleProperty> selection;
    if(_onlySelected){
        selection = FrameAdapter::findColumnShared(frame, {"Selection", "selection"}, DataType::Int);
        if(!selection){
            spdlog::warn("Voronoi: --onlySelected requested but no Selection column found; treating all atoms as selected");
        }
    }

    VoronoiAnalysisEngine::Params params;
    params.cutoff = _cutoff;
    params.edgeThreshold = _edgeThreshold;
    params.faceThreshold = _faceThreshold;
    params.useRadii = _useRadii && radii != nullptr;

    VoronoiAnalysisEngine engine(
        positions.get(),
        frame.simulationCell,
        radii.get(),
        selection.get(),
        params
    );
    VoronoiAnalysisLog log;
    ProgressSteps steps({
        {"Computing the Voronoi tessellation", 3.0},
        {"Writing outputs", 2.0},
    });
    engine.perform(log);
    log.dump(engine);
    steps.next();

    auto volumes = engine.atomicVolumes();
    auto coords = engine.coordNumbers();
    auto cavityRadii = engine.cavityRadii();

    double volumeSum = 0.0;
    long long coordSum = 0;
    const int atomCount = frame.natoms;
    for(int i = 0; i < atomCount; ++i){
        volumeSum += volumes->getDouble(i);
        coordSum += coords->getInt(i);
    }
    const double meanVolume = atomCount > 0 ? volumeSum / atomCount : 0.0;
    const double meanCoord = atomCount > 0 ? static_cast<double>(coordSum) / atomCount : 0.0;

    json result;
    result["main_listing"] = {
        {"total_atoms", atomCount},
        {"mean_atomic_volume", meanVolume},
        {"mean_coordination", meanCoord},
        {"effective_cutoff", engine.effectiveCutoff()},
        {"edge_threshold", params.edgeThreshold},
        {"face_threshold", params.faceThreshold},
        {"use_radii", params.useRadii},
        {"mean_face_order", engine.meanFaceOrder()},
        {"max_face_order_observed", engine.maxFaceOrderObserved()},
        {"polyhedra_mesh_count", engine.polyhedraMeshCount()},
    };
    AnalysisResult::addTiming(result, startTime);
    result["is_failed"] = false;

    if(!outputBase.empty()){
        log.dump(JsonUtils::writeSummaryTable(result, outputBase, "voronoi", "Voronoi summary"));

        const WrittenTable atoms{"Voronoi atoms", outputPath(outputBase, "atoms")};
        const auto& faceIndices = engine.faceIndices();
        streamAtomsToParquet(
            atoms.path,
            frame,
            [&coords](std::size_t i){
                return "Coordination_" + std::to_string(coords->getInt(i));
            },
            [&](ColumnarAtomWriter& w, std::size_t i){
                w.field("coordination",  coords->getInt(i));
                w.field("atomic_volume", volumes->getDouble(i));
                w.field("cavity_radius", cavityRadii->getDouble(i));
                if(i < faceIndices.size()){
                    const auto& fi = faceIndices[i];
                    std::vector<double> fiDouble(fi.begin(), fi.end());
                    w.field("face_indices", fiDouble);
                    w.field("max_face_order", static_cast<std::int64_t>(fi.empty() ? 0 : fi.back()));
                }
            },
            /*resolveStructureId=*/{},
            /*includeStructureColumns=*/false
        );
        log.dump(atoms);
    }

    return result;
}

}
