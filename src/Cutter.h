
#pragma once

#include "geometrycentral/surface/manifold_surface_mesh.h"

#include <array>
#include <map>
#include <queue>
#include <set>
#include <vector>

#include "Utils.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace CEPS {
namespace ImplementationDetails {

// ==========================================================================
//                         Compute Mesh Cuts
// ==========================================================================

std::set<Edge> goodCuts(ManifoldSurfaceMesh& mesh,
                        const std::vector<Vertex>& cones,
                        const EdgeData<double>& edgeCost);

std::set<Edge> goodCuts(ManifoldSurfaceMesh& mesh,
                        const std::vector<Vertex>& cones,
                        const EdgeData<double>& edgeCost,
                        std::set<Edge> startingCut);

std::set<Edge> bestCuts(ManifoldSurfaceMesh& mesh,
                        const std::vector<Vertex>& cones,
                        const EdgeData<double>& edgeCost);

std::set<Edge> bestCuts(ManifoldSurfaceMesh& mesh,
                        const std::vector<Vertex>& cones,
                        const EdgeData<double>& edgeCost,
                        std::set<Edge> startingCut);

std::set<Edge> cutOutFront(ManifoldSurfaceMesh& mesh,
                           const std::set<Vertex>& frontVertices,
                           const std::vector<Vertex>& boundaryVertices);

std::set<Edge> shortestPath(ManifoldSurfaceMesh& mesh,
                            const EdgeData<double>& weight, Vertex src,
                            Vertex dst);

VertexData<double> distances(ManifoldSurfaceMesh& mesh,
                             const EdgeData<double>& weight, Vertex src);

// Compute vertex indices of corners in cut mesh
// Useful for building Laplacian on cut mesh
// Returns the new indices, and the number of distinct indices (i.e. number of
// vertices in the cut mesh)
std::tuple<CornerData<size_t>, size_t> indexCorners(ManifoldSurfaceMesh& mesh,
                                                    const std::set<Edge>& cut);


// ==========================================================================
//                 Homology Generators via Tree-Cotree
// ==========================================================================

// Encodes the primal spanning tree via a VertexData<Halfedge> primalParent
// where primalParent[v].vertex() is v's parent, and primalParent[v].edge() is
// the edge from the parent to v

// Similarly, encodes the dual spanning tree via a Faedge> dualParent
// where dualParent[f].face() is f's parent and dualParent[f].edge() is the
// (primal) edge between f and its parent

VertexData<Halfedge> buildPrimalSpanningTree(ManifoldSurfaceMesh& mesh);
FaceData<Halfedge>
buildDualSpanningTree(ManifoldSurfaceMesh& mesh,
                      const VertexData<Halfedge>& primalParent);

VertexData<Halfedge> buildPrimalSpanningTree(ManifoldSurfaceMesh& mesh,
                                             const EdgeData<double>& edgeCost);
FaceData<Halfedge>
buildDualSpanningTree(ManifoldSurfaceMesh& mesh,
                      const VertexData<Halfedge>& primalParent,
                      const EdgeData<double>& edgeCost);

std::vector<std::vector<Edge>>
computeHomologyGenerators(ManifoldSurfaceMesh& mesh);

// Compute greedy homology generators according to
// http://ddg.dreamhosters.com/CS1762014/wp-content/uploads/2014/02/Paper-Presentation.pdf
// Slide 26 [Erickson?]

std::vector<std::vector<Edge>>
computeHomologyGenerators(ManifoldSurfaceMesh& mesh,
                          const EdgeData<double>& edgeCost);

bool inPrimalSpanningTree(Edge e, const VertexData<Halfedge>& primalParent);
bool inDualSpanningTree(Edge e, const FaceData<Halfedge>& dualParent);

std::vector<Edge> primalTreeLoop(Edge e,
                                 const VertexData<Halfedge>& primalParent);


// ==========================================================================
//                       Data Structure for Cuts
// ==========================================================================
// Stolen from boundary-first-flattening
typedef std::tuple<double, Vertex, Halfedge> VertexEntry;
class DisjointSets {
  public:
    // constructor
    DisjointSets(size_t n_) : parent(n_ + 1), rank(n_ + 1), marked(n_ + 1) {
        // initialize all elements to be in different sets and to have rank 0
        for (size_t i = 0; i <= n_; i++) {
            rank[i]   = 0;
            parent[i] = i;
            marked[i] = false;
        }
    }

    // find parent of element x
    int find(size_t x) {
        if (x != parent[x]) parent[x] = find(parent[x]);
        return parent[x];
    }

    // union by rank
    // if either set is marked, then the result is marked
    void merge(size_t x, size_t y) {
        x = find(x);
        y = find(y);

        // smaller tree becomes a subtree of the larger tree
        if (rank[x] > rank[y])
            parent[y] = x;
        else
            parent[x] = y;

        if (rank[x] == rank[y]) rank[y]++;

        // if either set was marked, both are marked
        if (marked[x] || marked[y]) {
            marked[x] = true;
            marked[y] = true;
        }
    }

    // mark set
    void mark(size_t x) { marked[find(x)] = true; }

    // unmark set
    void unmark(size_t x) { marked[find(x)] = false; }

    // check if set is marked
    bool isMarked(size_t x) { return marked[find(x)]; }

  private:
    // members
    std::vector<size_t> parent;
    std::vector<size_t> rank;
    std::vector<bool> marked;
};

} // namespace ImplementationDetails
} // namespace CEPS
