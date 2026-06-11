#pragma once

#include <volt/core/volt.h>
#include <volt/core/lammps_parser.h>
#include <volt/core/particle_property.h>
#include <volt/voronoi_engine.h>
#include <nlohmann/json.hpp>
#include <string>

namespace Volt{

using json = nlohmann::json;

/**
 * High-level Voronoi analysis driver.
 *
 * Wraps VoronoiAnalysisEngine with LAMMPS-frame IO and per-atom Parquet export
 * (a summary table plus the AtomisticExporter grouping-by-coordination atoms
 * table written via streamAtomsToParquet).
 */
class VoronoiService{
public:
    VoronoiService();

    void setCutoff(double cutoff);
    void setEdgeThreshold(double edgeThreshold);
    void setFaceThreshold(double faceThreshold);
    void setComputeIndices(bool computeIndices);
    void setEdgeCount(int edgeCount);
    void setUseRadii(bool useRadii);
    void setOnlySelected(bool onlySelected);

    json compute(
        const LammpsParser::Frame& frame,
        const std::string& outputBase = ""
    );

private:
    double _cutoff;
    double _edgeThreshold;
    double _faceThreshold;
    bool _computeIndices;
    int _edgeCount;
    bool _useRadii;
    bool _onlySelected;
};

}
