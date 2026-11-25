#include <volt/voronoi_analysis_log.h>
#include <spdlog/spdlog.h>

namespace Volt{

void VoronoiAnalysisLog::dump(const VoronoiAnalysisEngine::AutoCutoff& cutoff){
    spdlog::info("Voronoi: auto-selected cutoff = {:.4f}", cutoff.radius);
}

void VoronoiAnalysisLog::dump(const VoronoiAnalysisEngine& engine){
    const std::string elapsed = phaseElapsed();
    const size_t n = engine.positions()->size();
    const ParticleProperty* selection = engine.selection();
    size_t cellCount = n;
    if(selection){
        cellCount = 0;
        for(size_t i = 0; i < n; ++i){
            if(selection->getInt(i) != 0) ++cellCount;
        }
    }
    spdlog::info("Computed {} Voronoi cells in {}", cellCount, elapsed);

    const size_t validCells = static_cast<size_t>(engine.polyhedraMeshCount());
    if(validCells < cellCount){
        spdlog::warn("{} of {} Voronoi cells have no neighbor faces left (coordination 0)", cellCount - validCells, cellCount);
    }

    const ParticleProperty& volumes = *engine.atomicVolumes();
    const ParticleProperty& coords = *engine.coordNumbers();
    double volumeSum = 0.0;
    long long coordSum = 0;
    for(size_t i = 0; i < n; ++i){
        volumeSum += volumes.getDouble(i);
        coordSum += coords.getInt(i);
    }
    const double meanVolume = n > 0 ? volumeSum / n : 0.0;
    const double meanCoord = n > 0 ? static_cast<double>(coordSum) / n : 0.0;

    const double cellVolume = engine.cell().volume3D();
    spdlog::info(
        "Voronoi statistics: mean atomic volume {:.4g} (cells fill {:.2f}% of the simulation cell), "
        "mean coordination {:.2f}, mean face order {:.2f}, max face order {}",
        meanVolume, cellVolume > 0.0 ? 100.0 * volumeSum / cellVolume : 0.0, meanCoord,
        engine.meanFaceOrder(), engine.maxFaceOrderObserved()
    );
}

}
