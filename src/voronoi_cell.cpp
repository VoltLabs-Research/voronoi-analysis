#include <volt/voronoi_cell.h>
#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <unordered_map>
#include <utility>

namespace Volt{

namespace {
// Tolerance for coplanarity / edge-crossing decisions. All clipping is done
// in the central-atom-local frame where typical distances are O(few Å).
constexpr double kEps = 1e-10;
}

void VoronoiCell::initCube(double halfSide){
    _vertices.clear();
    _faces.clear();
    _vertices.reserve(8);
    _faces.reserve(6);

    // Eight cube corners, labeled so that the face lists below are CCW when
    // viewed from outside the cube.
    _vertices.emplace_back(-halfSide, -halfSide, -halfSide); // 0
    _vertices.emplace_back( halfSide, -halfSide, -halfSide); // 1
    _vertices.emplace_back( halfSide,  halfSide, -halfSide); // 2
    _vertices.emplace_back(-halfSide,  halfSide, -halfSide); // 3
    _vertices.emplace_back(-halfSide, -halfSide,  halfSide); // 4
    _vertices.emplace_back( halfSide, -halfSide,  halfSide); // 5
    _vertices.emplace_back( halfSide,  halfSide,  halfSide); // 6
    _vertices.emplace_back(-halfSide,  halfSide,  halfSide); // 7

    auto addFace = [this](std::initializer_list<int> verts,
                          const Vector3& n){
        Face f;
        f.vertices.assign(verts.begin(), verts.end());
        f.normal = n;
        f.neighborIndex = -1;
        _faces.push_back(std::move(f));
    };

    // -Z face (viewed from outside, i.e. from -Z):
    addFace({0, 3, 2, 1}, Vector3(0.0, 0.0, -1.0));
    // +Z face
    addFace({4, 5, 6, 7}, Vector3(0.0, 0.0, 1.0));
    // -Y face
    addFace({0, 1, 5, 4}, Vector3(0.0, -1.0, 0.0));
    // +Y face
    addFace({2, 3, 7, 6}, Vector3(0.0, 1.0, 0.0));
    // -X face
    addFace({0, 4, 7, 3}, Vector3(-1.0, 0.0, 0.0));
    // +X face
    addFace({1, 2, 6, 5}, Vector3(1.0, 0.0, 0.0));
}

bool VoronoiCell::clip(const Vector3& n, double d, int neighborIndex){
    if(_vertices.empty() || _faces.empty()){
        return false;
    }

    const size_t nVerts = _vertices.size();

    // Classify each vertex by signed distance to the plane n.x = d.
    // Side: INSIDE (s < -eps, strictly kept), OUTSIDE (s > eps, dropped),
    // ON (|s| <= eps, lies on plane, kept and also part of cap boundary).
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

    // No vertex strictly outside: plane does not cut the cell.
    if(outsideCount == 0){
        return true;
    }

    // No vertex strictly inside (all outside or on plane): cell becomes empty.
    if(insideCount == 0){
        _vertices.clear();
        _faces.clear();
        return false;
    }

    // We rebuild the vertex list and face list. Vertices classified as
    // INSIDE or ON are kept (ON verts lie on the new cap plane and become
    // shared between old faces and the new cap face). Each edge that has
    // one INSIDE endpoint and one OUTSIDE endpoint produces exactly one new
    // intersection vertex; cached by undirected edge key so adjacent faces
    // share the same vertex index.
    std::vector<Vector3> newVertices;
    newVertices.reserve(_vertices.size() + 8);

    // old-vertex-index -> new-vertex-index (for vertices we keep as-is).
    std::vector<int> vertexRemap(nVerts, -1);
    for(size_t i = 0; i < nVerts; ++i){
        if(sides[i] != OUTSIDE){
            vertexRemap[i] = static_cast<int>(newVertices.size());
            newVertices.push_back(_vertices[i]);
        }
    }

    // Edge -> intersection-vertex cache. Key is packed (min,max) of the two
    // old vertex indices. Only used for IN<->OUT crossings.
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
        // Only called for a strict IN<->OUT crossing so (db - da) has the
        // same sign as (db - da) and is >= 2*kEps in magnitude.
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

    // Rebuild each face. A face that lies entirely on the outside disappears.
    std::vector<Face> newFaces;
    newFaces.reserve(_faces.size() + 1);

    // For the cap polygon we collect one directed segment per clipped face.
    // On the retained polygon we walk CCW from the face's outward normal,
    // and the cap-plane segment runs exit -> entry. The cap polygon (when
    // walked CCW from its own outward normal +n) traverses the shared edge
    // in the opposite direction, so we store it as {entry, exit}. Chained
    // together these yield the cap polygon in CCW order.
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

            // Emit `cur` into the retained loop if it is kept.
            if(sCur != OUTSIDE){
                resultLoop.push_back(vertexRemap[cur]);
            }

            // Classify the transition from cur to nxt and emit any new
            // intersection vertex plus record entry/exit.
            if(sCur == INSIDE && sNxt == OUTSIDE){
                // IN -> OUT: true crossing. Intersection = exit.
                const int cut = getOrInsertIntersection(cur, nxt);
                resultLoop.push_back(cut);
                exitVertex = cut;
            }else if(sCur == OUTSIDE && sNxt == INSIDE){
                // OUT -> IN: true crossing. Intersection = entry.
                const int cut = getOrInsertIntersection(cur, nxt);
                resultLoop.push_back(cut);
                entryVertex = cut;
            }else if(sCur == INSIDE && sNxt == ON){
                // IN -> ON: no new intersection (nxt is on plane and will
                // be pushed in the next iteration). nxt is the exit.
                exitVertex = vertexRemap[nxt];
            }else if(sCur == ON && sNxt == OUTSIDE){
                // ON -> OUT: cur is already in resultLoop. cur is the exit.
                exitVertex = vertexRemap[cur];
            }else if(sCur == OUTSIDE && sNxt == ON){
                // OUT -> ON: nxt will be pushed next. nxt is the entry.
                entryVertex = vertexRemap[nxt];
            }else if(sCur == ON && sNxt == INSIDE){
                // ON -> IN: cur is already in resultLoop. cur is the entry.
                entryVertex = vertexRemap[cur];
            }else if(sCur == ON && sNxt == ON){
                // Edge lies entirely on the clip plane. Use it directly as
                // a cap edge for this face (going cur -> nxt along the plane
                // is the cap boundary segment). Record exit=cur, entry=nxt;
                // the cap traversal direction will follow from stitching.
                exitVertex = vertexRemap[cur];
                entryVertex = vertexRemap[nxt];
            }
            // IN->IN, OUT->OUT, ON->same are handled implicitly (no action).
        }

        // Collapse consecutive duplicates and close-duplicate (wraparound).
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
            // Degenerate face (tangential touch or sliver). Drop it.
            continue;
        }

        Face out;
        out.vertices = std::move(cleaned);
        out.normal = face.normal;
        out.neighborIndex = face.neighborIndex;
        newFaces.push_back(std::move(out));

        if(entryVertex >= 0 && exitVertex >= 0 && entryVertex != exitVertex){
            // Within a retained face, the cap-plane segment is traversed
            // exit -> entry (as we walk the retained polygon CCW). For the
            // cap polygon to be CCW from its own outward normal (+n), the
            // shared edge must be traversed in the OPPOSITE direction on the
            // cap: entry -> exit. So push {entry, exit} for stitching.
            capEdges.push_back({entryVertex, exitVertex});
        }
    }

    _vertices = std::move(newVertices);
    _faces = std::move(newFaces);

    // Stitch cap edges into a single closed polygon. Each intersection /
    // on-plane vertex that participates in the cap appears in exactly one
    // capEdge as "from" and exactly one as "to", so the chain is well-defined.
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
            // The stitched loop follows `{entry -> exit}` chaining, which
            // matches the cap polygon walked CCW as seen from +n (the cap's
            // outward normal). This keeps every face in the cell oriented
            // CCW from its stored outward normal, which is the invariant
            // relied upon by subsequent clips' entry/exit bookkeeping.
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

    // Sum signed tetrahedra from the origin (central atom) to each
    // fan-triangulated face triangle. The sign cancels correctly for a
    // closed convex polytope regardless of face orientation.
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
