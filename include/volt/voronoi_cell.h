#pragma once

#include <volt/core/volt.h>
#include <volt/math/vector3.h>
#include <vector>
#include <cstddef>

namespace Volt{

class VoronoiCell{
public:
    struct Face{
        std::vector<int> vertices;
        Vector3 normal;
        int neighborIndex;
    };

    VoronoiCell() = default;

    void initCube(double halfSide);

    bool clip(const Vector3& n, double d, int neighborIndex);

    double volume() const;

    double cavityRadius() const;

    bool isEmpty() const noexcept {
        return _vertices.empty() || _faces.empty();
    }

    const std::vector<Vector3>& vertices() const noexcept { return _vertices; }
    const std::vector<Face>& faces() const noexcept { return _faces; }

    double faceArea(const Face& face) const;

    double faceMinEdgeLength(const Face& face) const;

    static int faceOrder(const Face& face) noexcept {
        return static_cast<int>(face.vertices.size());
    }

private:
    std::vector<Vector3> _vertices;
    std::vector<Face> _faces;
};

}
