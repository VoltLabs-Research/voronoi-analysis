#pragma once

#include <volt/core/volt.h>
#include <volt/math/vector3.h>
#include <vector>
#include <cstddef>

namespace Volt{

/**
 * A convex polytope represented as a vertex/face (indexed) mesh, suitable for
 * iterative half-space clipping. The cell is initialized as a cube and then
 * clipped against the perpendicular bisector of each central-neighbor vector
 * to produce the central atom's Voronoi region.
 *
 * This is a minimal, self-contained Voronoi implementation (no external
 * voro++ dependency). It is accurate enough for typical MD environments but
 * not intended for publication-grade topology statistics in degenerate
 * geometries - use OVITO's voro++-backed modifier for those cases.
 */
class VoronoiCell{
public:
    // A single polygonal face of the cell.
    struct Face{
        // Ordered vertex indices (CCW when viewed from outside the cell).
        std::vector<int> vertices;
        // Outward-pointing normal of the plane that introduced this face.
        // (For the six bootstrap cube faces this is the axis-aligned normal.)
        Vector3 normal;
        // Index of the neighbor atom that produced this face, or -1 for
        // bootstrap/bounding-box faces that never got replaced by a neighbor.
        int neighborIndex;
    };

    // Constructs an empty cell. Call initCube() before using.
    VoronoiCell() = default;

    // Initializes this cell as an axis-aligned cube of half-side `halfSide`
    // centered on the origin (the central atom is always treated as the
    // origin in clip() coordinates).
    void initCube(double halfSide);

    // Clips the current cell by the half-space
    //     n . x <= d
    // (i.e. keeps points with n.x <= d; discards anything on the other side).
    //
    // `neighborIndex` is stored on the newly-introduced face so we can later
    // map faces back to neighbor atoms. If the entire cell lies inside the
    // half-space this is a no-op; if the entire cell lies outside the cell
    // becomes empty.
    //
    // Returns true if the cell still has volume after the clip.
    bool clip(const Vector3& n, double d, int neighborIndex);

    // Returns the total volume of the cell by fan-triangulating each face
    // from its first vertex and summing signed tetrahedra from the origin.
    double volume() const;

    // Returns true if the cell currently has zero vertices or faces.
    bool isEmpty() const noexcept {
        return _vertices.empty() || _faces.empty();
    }

    const std::vector<Vector3>& vertices() const noexcept { return _vertices; }
    const std::vector<Face>& faces() const noexcept { return _faces; }

    // Returns the area of a face. Fan-triangulated from the first vertex.
    double faceArea(const Face& face) const;

    // Returns the minimum edge length of a face (used for edgeThreshold).
    double faceMinEdgeLength(const Face& face) const;

    // Returns the number of edges (== number of vertices) of a face.
    static int faceOrder(const Face& face) noexcept {
        return static_cast<int>(face.vertices.size());
    }

private:
    std::vector<Vector3> _vertices;
    std::vector<Face> _faces;
};

}
