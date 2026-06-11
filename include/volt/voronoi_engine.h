#pragma once

#include <volt/core/particle_property.h>
#include <volt/core/simulation_cell.h>
#include <volt/analysis/cutoff_neighbor_finder.h>
#include <memory>

namespace Volt{

using Volt::Particles::ParticleProperty;
using Volt::Particles::SimulationCell;

/**
 * Per-atom Voronoi tessellation engine.
 *
 * For each central atom i we:
 *   1. Build a neighbor list within `cutoff` (auto-estimated from the
 *      nearest-neighbor distance if the user-provided cutoff is <= 0).
 *   2. Seed a cubic cell of half-side `cutoff` around atom i.
 *   3. Clip the cell by the perpendicular bisector (or radical plane when
 *      useRadii is true) of each neighbor vector.
 *   4. Count neighbor faces, optionally filtered by face area
 *      (faceThreshold) and minimum edge length (edgeThreshold), producing
 *      a topological coordination number.
 *   5. Compute the cell volume (atomic volume) and the cavity radius (the
 *      largest empty sphere centred on the atom that fits inside the cell).
 */
class VoronoiAnalysisEngine{
public:
    struct Params{
        double cutoff = 0.0;
        double edgeThreshold = 0.0;
        double faceThreshold = 0.0;
        bool useRadii = false;
    };

    VoronoiAnalysisEngine(
        ParticleProperty* positions,
        const SimulationCell& simCell,
        ParticleProperty* radii,
        ParticleProperty* selection,
        const Params& params
    );

    // Runs the tessellation. After this returns all output properties are
    // fully populated.
    void perform();

    // The cutoff that was actually used (either user-provided or
    // auto-estimated from the mean nearest-neighbor distance).
    double effectiveCutoff() const noexcept { return _effectiveCutoff; }

    // Outputs (allocated by the engine, always non-null after construction).
    std::shared_ptr<ParticleProperty> atomicVolumes() const { return _atomicVolumes; }
    std::shared_ptr<ParticleProperty> coordNumbers() const { return _coordNumbers; }
    std::shared_ptr<ParticleProperty> cavityRadii() const { return _cavityRadii; }

    ParticleProperty* positions() const { return _positions; }
    const SimulationCell& cell() const { return _simCell; }
    const Params& params() const { return _params; }

private:
    // Estimates a reasonable cutoff from the first-shell distance when the
    // user left `cutoff` at 0.
    double estimateAutoCutoff();

    // Processes a single central atom, writing results into the output
    // properties.
    void processAtom(size_t atomIndex, CutoffNeighborFinder& finder);

    ParticleProperty* _positions;
    SimulationCell _simCell;
    ParticleProperty* _radii;
    ParticleProperty* _selection;
    Params _params;
    double _effectiveCutoff;

    std::shared_ptr<ParticleProperty> _atomicVolumes;
    std::shared_ptr<ParticleProperty> _coordNumbers;
    std::shared_ptr<ParticleProperty> _cavityRadii;
};

}
