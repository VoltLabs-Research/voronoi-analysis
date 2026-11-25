#pragma once

#include <volt/core/run_log.h>
#include <volt/voronoi_engine.h>

namespace Volt{

// The execution log of one tessellation: each phase hands it what it built.
class VoronoiAnalysisLog : public RunLog{
public:
    using RunLog::dump;

    void dump(const VoronoiAnalysisEngine::AutoCutoff& cutoff);
    // The cells the engine computed, the ones left without faces, and their statistics.
    void dump(const VoronoiAnalysisEngine& engine);
};

}
