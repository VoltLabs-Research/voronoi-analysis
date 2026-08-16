#include <volt/voronoi_cell.h>
#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <unordered_map>
#include <utility>

namespace Volt{

constexpr double kEps = 1e-10;

void VoronoiCell::initCube(double halfSide){
    _vertices.clear();
    _faces.clear();
    _vertices.reserve(8);
    _faces.reserve(6);

    _vertices.emplace_back(-halfSide, -halfSide, -halfSide);
    _vertices.emplace_back( halfSide, -halfSide, -halfSide);
    _vertices.emplace_back( halfSide,  halfSide, -halfSide);
    _vertices.emplace_back(-halfSide,  halfSide, -halfSide);
    _vertices.emplace_back(-halfSide, -halfSide,  halfSide);
    _vertices.emplace_back( halfSide, -halfSide,  halfSide);
    _vertices.emplace_back( halfSide,  halfSide,  halfSide);
    _vertices.emplace_back(-halfSide,  halfSide,  halfSide);

    auto addFace = [this](std::initializer_list<int> verts,
                          const Vector3& n){
        Face f;
        f.vertices.assign(verts.begin(), verts.end());
        f.normal = n;
        f.neighborIndex = -1;
        _faces.push_back(std::move(f));
    };

    addFace({0, 3, 2, 1}, Vector3(0.0, 0.0, -1.0));
    addFace({4, 5, 6, 7}, Vector3(0.0, 0.0, 1.0));
    addFace({0, 1, 5, 4}, Vector3(0.0, -1.0, 0.0));
    addFace({2, 3, 7, 6}, Vector3(0.0, 1.0, 0.0));
    addFace({0, 4, 7, 3}, Vector3(-1.0, 0.0, 0.0));
    addFace({1, 2, 6, 5}, Vector3(1.0, 0.0, 0.0));
}

bool VoronoiCell::clip(const Vector3& n, double d, int neighborIndex){
    if(_vertices.empty() || _faces.empty()){
        return false;
    }

    const size_t nVerts = _vertices.size();

    enum Side : std::uint8_t { INSIDE = 0, ON = 1, OUTSIDE = 2 };
    std::vector<double> distances(nVerts);
    std::vector<Side> sides(nVerts);
    int insideCount = 0;
    int outsideCount = 0;
    int onCount = 0;
    for(size_t i = 0; i < nVerts; ++i){
        const double s = n.dot(_vertices[i]) - d;
        distances[i] = s;
        if(s > kEps){ sides[i] = OUTSIDE; ++outsideCount; }
        else if(s < -kEps){ sides[i] = INSIDE; ++insideCount; }
        else{ sides[i] = ON; ++onCount; }
    }

    if(outsideCount == 0){
        return true;
    }

    if(insideCount == 0){
        _vertices.clear();
        _faces.clear();
        return false;
    }

    std::vector<Vector3> newVertices;
    newVertices.reserve(_vertices.size() + 8);

    std::vector<int> vertexRemap(nVerts, -1);
    for(size_t i = 0; i < nVerts; ++i){
        if(sides[i] != OUTSIDE){
            vertexRemap[i] = static_cast<int>(newVertices.size());
            newVertices.push_back(_vertices[i]);
        }
    }

    std::unordered_map<std::uint64_t, int> edgeIntersection;
    auto edgeKey = [](int a, int b) -> std::uint64_t {
        if(a > b) std::swap(a, b);
        return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(a)) << 32)
             | static_cast<std::uint64_t>(static_cast<std::uint32_t>(b));
    };

    auto getOrInsertIntersection = [&](int a, int b) -> int {
        const auto key = edgeKey(a, b);
        auto it = edgeIntersection.find(key);
        if(it != edgeIntersection.end()) return it->second;

        const double da = distances[a];
        const double db = distances[b];
        const double t = da / (da - db);
        const Vector3& pa = _vertices[a];
        const Vector3& pb = _vertices[b];
        Vector3 p(pa.x() + t * (pb.x() - pa.x()),
                  pa.y() + t * (pb.y() - pa.y()),
                  pa.z() + t * (pb.z() - pa.z()));
        const int newIdx = static_cast<int>(newVertices.size());
        newVertices.push_back(p);
        edgeIntersection.emplace(key, newIdx);
        return newIdx;
    };

    std::vector<Face> newFaces;
    newFaces.reserve(_faces.size() + 1);

    struct CapEdge { int from; int to; };
    std::vector<CapEdge> capEdges;
    capEdges.reserve(_faces.size());

    for(const Face& face : _faces){
        const int m = static_cast<int>(face.vertices.size());
        std::vector<int> resultLoop;
        resultLoop.reserve(m + 2);

        int entryVertex = -1;
        int exitVertex = -1;

        for(int i = 0; i < m; ++i){
            const int cur = face.vertices[i];
            const int nxt = face.vertices[(i + 1) % m];
            const Side sCur = sides[cur];
            const Side sNxt = sides[nxt];

            if(sCur != OUTSIDE){
                resultLoop.push_back(vertexRemap[cur]);
            }

            if(sCur == INSIDE && sNxt == OUTSIDE){
                const int cut = getOrInsertIntersection(cur, nxt);
                resultLoop.push_back(cut);
                exitVertex = cut;
            }else if(sCur == OUTSIDE && sNxt == INSIDE){
                const int cut = getOrInsertIntersection(cur, nxt);
                resultLoop.push_back(cut);
                entryVertex = cut;
            }else if(sCur == INSIDE && sNxt == ON){
                exitVertex = vertexRemap[nxt];
            }else if(sCur == ON && sNxt == OUTSIDE){
                exitVertex = vertexRemap[cur];
            }else if(sCur == OUTSIDE && sNxt == ON){
                entryVertex = vertexRemap[nxt];
            }else if(sCur == ON && sNxt == INSIDE){
                entryVertex = vertexRemap[cur];
            }else if(sCur == ON && sNxt == ON){
                exitVertex = vertexRemap[cur];
                entryVertex = vertexRemap[nxt];
            }
        }

        std::vector<int> cleaned;
        cleaned.reserve(resultLoop.size());
        for(int v : resultLoop){
            if(cleaned.empty() || cleaned.back() != v){
                cleaned.push_back(v);
            }
        }
        if(cleaned.size() >= 2 && cleaned.front() == cleaned.back()){
            cleaned.pop_back();
        }
        if(cleaned.size() < 3){
            continue;
        }

        Face out;
        out.vertices = std::move(cleaned);
        out.normal = face.normal;
        out.neighborIndex = face.neighborIndex;
        newFaces.push_back(std::move(out));

        if(entryVertex >= 0 && exitVertex >= 0 && entryVertex != exitVertex){
            capEdges.push_back({entryVertex, exitVertex});
        }
    }

    _vertices = std::move(newVertices);
    _faces = std::move(newFaces);

    if(capEdges.size() >= 3){
        std::unordered_map<int, int> fromMap;
        fromMap.reserve(capEdges.size() * 2);
        for(const CapEdge& e : capEdges){
            fromMap.emplace(e.from, e.to);
        }

        std::vector<int> capLoop;
        capLoop.reserve(capEdges.size());
        int start = capEdges.front().from;
        int current = start;
        for(size_t step = 0; step < capEdges.size(); ++step){
            capLoop.push_back(current);
            auto it = fromMap.find(current);
            if(it == fromMap.end()) break;
            current = it->second;
            if(current == start) break;
        }

        std::vector<int> cleaned;
        cleaned.reserve(capLoop.size());
        for(int v : capLoop){
            if(cleaned.empty() || cleaned.back() != v) cleaned.push_back(v);
        }
        if(cleaned.size() >= 2 && cleaned.front() == cleaned.back()){
            cleaned.pop_back();
        }

        if(cleaned.size() >= 3){
            Face cap;
            cap.vertices = std::move(cleaned);
            cap.normal = n;
            cap.neighborIndex = neighborIndex;
            _faces.push_back(std::move(cap));
        }
    }

    return !_faces.empty() && !_vertices.empty();
}

double VoronoiCell::volume() const {
    if(_vertices.empty() || _faces.empty()) return 0.0;

    double v = 0.0;
    for(const Face& face : _faces){
        const int m = static_cast<int>(face.vertices.size());
        if(m < 3) continue;
        const Vector3& v0 = _vertices[face.vertices[0]];
        for(int i = 1; i + 1 < m; ++i){
            const Vector3& v1 = _vertices[face.vertices[i]];
            const Vector3& v2 = _vertices[face.vertices[i + 1]];
            v += v0.dot(v1.cross(v2));
        }
    }
    return std::abs(v) / 6.0;
}

double VoronoiCell::cavityRadius() const {
    if(_vertices.empty() || _faces.empty()) return 0.0;

    double minDist = std::numeric_limits<double>::infinity();
    for(const Face& face : _faces){
        const int m = static_cast<int>(face.vertices.size());
        if(m < 3) continue;
        const Vector3& v0 = _vertices[face.vertices[0]];

        Vector3 normal(0.0, 0.0, 0.0);
        for(int i = 1; i + 1 < m; ++i){
            const Vector3& v1 = _vertices[face.vertices[i]];
            const Vector3& v2 = _vertices[face.vertices[i + 1]];
            normal += (v1 - v0).cross(v2 - v0);
        }
        const double nLen = normal.length();
        if(nLen <= kEps) continue;

        const double dist = std::abs(v0.dot(normal) / nLen);
        if(dist < minDist) minDist = dist;
    }

    return std::isfinite(minDist) ? minDist : 0.0;
}

double VoronoiCell::faceArea(const Face& face) const {
    const int m = static_cast<int>(face.vertices.size());
    if(m < 3) return 0.0;
    Vector3 acc(0.0, 0.0, 0.0);
    const Vector3& v0 = _vertices[face.vertices[0]];
    for(int i = 1; i + 1 < m; ++i){
        const Vector3& v1 = _vertices[face.vertices[i]];
        const Vector3& v2 = _vertices[face.vertices[i + 1]];
        acc += (v1 - v0).cross(v2 - v0);
    }
    return 0.5 * acc.length();
}

double VoronoiCell::faceMinEdgeLength(const Face& face) const {
    const int m = static_cast<int>(face.vertices.size());
    if(m < 2) return 0.0;
    double minLen = std::numeric_limits<double>::infinity();
    for(int i = 0; i < m; ++i){
        const Vector3& a = _vertices[face.vertices[i]];
        const Vector3& b = _vertices[face.vertices[(i + 1) % m]];
        const double len = (b - a).length();
        if(len < minLen) minLen = len;
    }
    return std::isfinite(minLen) ? minLen : 0.0;
}

}
