#include <volt/plugin/plugin_entry.h>
#include <volt/voronoi_service.h>

using namespace Volt;
using namespace Volt::Plugin;
using S = VoronoiService;

static const std::vector<OptionBinding<S>> bindings = {
    opt("--cutoff",         "Neighbor search cutoff (0 = auto-estimate)",   0.0,   &S::setCutoff),
    opt("--edgeThreshold",  "Minimum face edge length to count a neighbor", 0.0,   &S::setEdgeThreshold),
    opt("--faceThreshold",  "Minimum face area to count a neighbor",        0.0,   &S::setFaceThreshold),
    opt("--computeIndices", "Emit per-atom Voronoi index vector",           false, &S::setComputeIndices),
    opt("--edgeCount",      "Max polygon order tracked in the index",       6,     &S::setEdgeCount),
    opt("--useRadii",       "Radical/power Voronoi from Radius column",      false, &S::setUseRadii),
    opt("--onlySelected",   "Only analyze atoms with non-zero Selection",   false, &S::setOnlySelected),
};

// NOTE: --threads is intentionally NOT a binding; pluginMain owns TBB global_control.
VOLT_SERVICE_PLUGIN("volt-voronoi", "Voronoi Analysis", S, bindings)
