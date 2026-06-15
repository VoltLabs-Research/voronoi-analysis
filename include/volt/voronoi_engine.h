#pragma once

#include <volt/core/particle_property.h>
#include <volt/core/simulation_cell.h>
#include <volt/analysis/cutoff_neighbor_finder.h>
#include <memory>
#include <vector>

namespace Volt{

using Volt::Particles::ParticleProperty;
using Volt::Particles::SimulationCell;

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

    void perform();

    double effectiveCutoff() const noexcept { return _effectiveCutoff; }

    std::shared_ptr<ParticleProperty> atomicVolumes() const { return _atomicVolumes; }
    std::shared_ptr<ParticleProperty> coordNumbers() const { return _coordNumbers; }
    std::shared_ptr<ParticleProperty> cavityRadii() const { return _cavityRadii; }

    // Per-atom face-order lists (index vector). Each element is a sorted list of
    // face orders (3,4,5,...) from the atom's valid Voronoi faces.
    const std::vector<std::vector<int>>& faceIndices() const { return _faceIndices; }

    // Aggregate face-order statistics (populated after perform()).
    double meanFaceOrder() const { return _meanFaceOrder; }
    int maxFaceOrderObserved() const { return _maxFaceOrderObserved; }
    int polyhedraMeshCount() const { return _polyhedraMeshCount; }

    ParticleProperty* positions() const { return _positions; }
    const SimulationCell& cell() const { return _simCell; }
    const Params& params() const { return _params; }

private:
    double estimateAutoCutoff();
    void processAtom(size_t atomIndex, CutoffNeighborFinder& finder);
    void computeAggregates();

    ParticleProperty* _positions;
    SimulationCell _simCell;
    ParticleProperty* _radii;
    ParticleProperty* _selection;
    Params _params;
    double _effectiveCutoff;

    std::shared_ptr<ParticleProperty> _atomicVolumes;
    std::shared_ptr<ParticleProperty> _coordNumbers;
    std::shared_ptr<ParticleProperty> _cavityRadii;

    // Per-atom face-order index vectors (one entry per atom, size = n)
    std::vector<std::vector<int>> _faceIndices;

    double _meanFaceOrder = 0.0;
    int _maxFaceOrderObserved = 0;
    int _polyhedraMeshCount = 0;
};

}
