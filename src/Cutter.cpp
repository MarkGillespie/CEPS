#include "Cutter.h"

namespace CEPS {
namespace ImplementationDetails {
// ==========================================================================
//                         Compute Mesh Cuts
// ==========================================================================

std::set<Edge> goodCuts(ManifoldSurfaceMesh& mesh,
                        const std::vector<Vertex>& cones,
                        const EdgeData<double>& edgeCost) {
    return goodCuts(mesh, cones, edgeCost, {});
}


std::set<Edge> goodCuts(ManifoldSurfaceMesh& mesh,
                        const std::vector<Vertex>& cones,
                        const EdgeData<double>& edgeCost,
                        std::set<Edge> startingCut) {

    // Add in homology generators to cut
    std::vector<std::vector<Edge>> homologyGenerators =
        computeHomologyGenerators(mesh, edgeCost);

    for (const std::vector<Edge>& generator : homologyGenerators) {
        for (Edge e : generator) startingCut.insert(e);
    }

    std::set<Vertex> unusedCones(std::begin(cones), std::end(cones));
    for (Edge e : startingCut) {
        unusedCones.erase(e.halfedge().vertex());
        unusedCones.erase(e.halfedge().twin().vertex());
    }

    std::set<Edge> cutEdges = startingCut;


    // initialize data structures
    size_t nV = (size_t)mesh.nVertices();
    DisjointSets ds(nV);
    Halfedge nullHe;
    VertexData<Halfedge> parent(mesh, nullHe);
    std::priority_queue<VertexEntry, std::vector<VertexEntry>,
                        std::greater<VertexEntry>>
        pq;

    VertexData<size_t> vIdx = mesh.getVertexIndices();

    // add starting cut to pq
    for (Edge e : startingCut) {
        ds.merge(vIdx[src(e)], vIdx[dst(e)]);
        std::array<Vertex, 2> neigh{src(e), dst(e)};
        for (Vertex v : neigh) {
            for (Halfedge vHe : v.outgoingHalfedges()) {
                pq.push(VertexEntry(edgeCost[vHe.edge()], vHe.twin().vertex(),
                                    vHe));
            }
            ds.mark(vIdx[v]);
        }
    }

    // add cone neighbors to pq
    for (size_t i = 0; i < (size_t)cones.size(); i++) {
        for (Halfedge he : cones[i].outgoingHalfedges()) {
            pq.push(VertexEntry(edgeCost[he.edge()], he.twin().vertex(), he));
        }

        // mark set
        ds.mark(vIdx[cones[i]]);
    }

    // construct approximate steiner tree
    while (!pq.empty()) {
        VertexEntry entry = pq.top();
        pq.pop();

        double weight = std::get<0>(entry);
        Halfedge he   = std::get<2>(entry);
        Vertex v1     = std::get<1>(entry);
        Vertex v2     = he.vertex();

        if (ds.find(vIdx[v1]) != ds.find(vIdx[v2])) {
            // if merging two marked sets, then mark edges that connect
            // these two regions as cut edges
            if (ds.isMarked(vIdx[v1]) && ds.isMarked(vIdx[v2])) {
                // one side
                Halfedge currHe = parent[v1];
                while (currHe != nullHe) {
                    cutEdges.insert(currHe.edge());
                    currHe = parent[currHe.vertex()];
                }

                // bridge
                cutEdges.insert(he.edge());

                // other side
                currHe = parent[v2];
                while (currHe != nullHe) {
                    cutEdges.insert(currHe.edge());
                    currHe = parent[currHe.vertex()];
                }
            }

            // record potential parent and merge sets
            parent[v1] = he;
            ds.merge(vIdx[v1], vIdx[v2]);

            for (Halfedge vHe : v1.outgoingHalfedges()) {
                pq.push(VertexEntry(edgeCost[vHe.edge()] + weight,
                                    vHe.twin().vertex(), vHe));
            }
        }
    }

    return cutEdges;
}

std::set<Edge> bestCuts(ManifoldSurfaceMesh& mesh,
                        const std::vector<Vertex>& cones,
                        const EdgeData<double>& edgeCost) {
    return bestCuts(mesh, cones, edgeCost, {});
}

std::set<Edge> bestCuts(ManifoldSurfaceMesh& mesh,
                        const std::vector<Vertex>& cones,
                        const EdgeData<double>& edgeCost,
                        std::set<Edge> startingCut) {

    std::set<Vertex> unusedCones(std::begin(cones), std::end(cones));
    for (Edge e : startingCut) {
        unusedCones.erase(e.halfedge().vertex());
        unusedCones.erase(e.halfedge().twin().vertex());
    }

    std::set<Edge> cutEdges = startingCut;
    std::set<Vertex> cutVertices;
    if (cutEdges.empty()) {
        cutVertices.insert(cones[0]);
        unusedCones.erase(cones[0]);
    } else {
        for (Edge e : cutEdges) {
            cutVertices.insert(e.halfedge().vertex());
            cutVertices.insert(e.halfedge().twin().vertex());
        }
    }

    std::map<Vertex, VertexData<double>> pathDistances;
    pathDistances[cones[0]] = distances(mesh, edgeCost, cones[0]);

    while (!unusedCones.empty()) {
        for (Vertex v : cutVertices) {
            if (pathDistances.find(v) == pathDistances.end()) {
                pathDistances[v] = distances(mesh, edgeCost, v);
            }
        }

        Vertex cheapestDst;
        Vertex cheapestDstSrc;
        double cheapestDstCost = std::numeric_limits<double>::infinity();
        for (Vertex dst : unusedCones) {
            Vertex bestSrc;
            double bestCost = std::numeric_limits<double>::infinity();
            for (Vertex src : cutVertices) {
                if (pathDistances[src][dst] < bestCost) {
                    bestCost = pathDistances[src][dst];
                    bestSrc  = src;
                }
            }

            if (bestCost < cheapestDstCost) {
                cheapestDstCost = bestCost;
                cheapestDst     = dst;
                cheapestDstSrc  = bestSrc;
            }
        }

        std::set<Edge> pathEdges =
            shortestPath(mesh, edgeCost, cheapestDst, cheapestDstSrc);

        for (Edge e : pathEdges) {
            cutEdges.insert(e);
            cutVertices.insert(e.halfedge().vertex());
            cutVertices.insert(e.halfedge().twin().vertex());
        }
        unusedCones.erase(cheapestDst);
    }

    return cutEdges;
}

std::set<Edge> cutOutFront(ManifoldSurfaceMesh& mesh,
                           const std::set<Vertex>& frontVertices,
                           const std::vector<Vertex>& boundaryVertices) {
    EdgeData<double> edgeCost(mesh, 0.01);
    for (Vertex v : frontVertices) {
        for (Edge e : v.adjacentEdges()) {
            edgeCost[e] *= 100;
        }
    }

    return goodCuts(mesh, boundaryVertices, edgeCost);
}

std::set<Edge> shortestPath(ManifoldSurfaceMesh& mesh,
                            const EdgeData<double>& weight, Vertex src,
                            Vertex dst) {
    VertexData<double> dist(mesh, std::numeric_limits<double>::infinity());
    VertexData<Halfedge> parent(mesh);
    VertexData<char> visited(mesh, false);
    dist[src] = 0;

    std::priority_queue<std::pair<double, Vertex>> visitQueue;
    for (Vertex v : mesh.vertices()) {
        visitQueue.push(std::make_pair(-dist[v], v));
    }

    while (!visitQueue.empty()) {
        Vertex curr = visitQueue.top().second;
        visitQueue.pop();
        if (visited[curr] || std::isinf(dist[curr]) || curr == dst) continue;

        for (Halfedge he : curr.incomingHalfedges()) {
            Vertex opp = he.vertex();
            if (!visited[opp]) {
                double cost = dist[curr] + weight[he.edge()];
                if (cost < dist[opp]) {
                    dist[opp]   = cost;
                    parent[opp] = he.twin();
                    visitQueue.push(std::make_pair(-cost, opp));
                }
            }
        }

        visited[curr] = true;
    }

    std::set<Edge> path;
    Vertex curr = dst;
    while (curr != src) {
        path.insert(parent[curr].edge());
        curr = parent[curr].vertex();
    }

    return path;
}

VertexData<double> distances(ManifoldSurfaceMesh& mesh,
                             const EdgeData<double>& weight, Vertex src) {
    VertexData<double> dist(mesh, std::numeric_limits<double>::infinity());
    VertexData<char> visited(mesh, false);
    dist[src] = 0;

    std::priority_queue<std::pair<double, Vertex>> visitQueue;
    for (Vertex v : mesh.vertices()) {
        visitQueue.push(std::make_pair(-dist[v], v));
    }

    while (!visitQueue.empty()) {
        Vertex curr = visitQueue.top().second;
        visitQueue.pop();
        if (visited[curr] || std::isinf(dist[curr])) continue;

        for (Halfedge he : curr.incomingHalfedges()) {
            Vertex opp = he.vertex();
            if (!visited[opp]) {
                double cost = dist[curr] + weight[he.edge()];
                if (cost < dist[opp]) {
                    dist[opp] = cost;
                    visitQueue.push(std::make_pair(-cost, opp));
                }
            }
        }

        visited[curr] = true;
    }

    return dist;
}

// Stolen from BFF
// Bff.cpp > BFFData::index Wedges();
std::tuple<CornerData<size_t>, size_t> indexCorners(ManifoldSurfaceMesh& mesh,
                                                    const std::set<Edge>& cut) {
    VertexData<bool> interior(mesh, true);
    EdgeData<bool> onCut(mesh, false);
    for (Edge e : cut) {
        interior[src(e)] = false;
        interior[dst(e)] = false;
        onCut[e]         = true;
    }

    size_t iV = 0;
    CornerData<size_t> cornerIndices(mesh);

    // Index interior vertices
    for (Vertex v : mesh.vertices()) {
        if (interior[v]) {
            for (Corner c : v.adjacentCorners()) cornerIndices[c] = iV;
            iV++;
        }
    }

    // Index boundary (cut) vertices
    for (Edge e : cut) {
        std::array<Halfedge, 2> hedges{e.halfedge(), e.halfedge().twin()};
        for (Halfedge he : hedges) {
            Halfedge curr = he;
            do {
                curr                         = curr.twin().next();
                cornerIndices[curr.corner()] = iV;
            } while (!onCut[curr.edge()]);
            iV++;
        }
    }

    return std::make_tuple(cornerIndices, iV);
}

// ==========================================================================
//                 Homology Generators via Tree-Coterie
// ==========================================================================
VertexData<Halfedge> buildPrimalSpanningTree(ManifoldSurfaceMesh& mesh) {

    // Initialize each vertex as its own parent
    VertexData<Halfedge> primalParent(mesh);
    for (Vertex v : mesh.vertices()) primalParent[v] = v.halfedge();

    VertexData<bool> inTree(mesh, false);

    for (Vertex v : mesh.vertices()) {
        if (!inTree[v]) {
            inTree[v] = true;

            Vertex root = v;
            std::queue<Vertex> toVisit;
            toVisit.push(root);

            while (!toVisit.empty()) {
                Vertex v = toVisit.front();
                toVisit.pop();

                for (Halfedge he : v.outgoingHalfedges()) {
                    Vertex w = he.twin().vertex();
                    if (primalParent[w].vertex() == w && w != root) {
                        primalParent[w] = he;
                        toVisit.push(w);
                        inTree[w] = true;
                    }
                }
            }
        }
    }

    return primalParent;
}

FaceData<Halfedge>
buildDualSpanningTree(ManifoldSurfaceMesh& mesh,
                      const VertexData<Halfedge>& primalParent) {
    // Initialize each face as its own parent
    FaceData<Halfedge> dualParent(mesh);
    for (Face f : mesh.faces()) dualParent[f] = f.halfedge();

    FaceData<bool> inTree(mesh, false);

    for (Face f : mesh.faces()) {
        if (!inTree[f]) {
            inTree[f] = true;

            Face root = f;
            std::queue<Face> toVisit;
            toVisit.push(root);

            while (!toVisit.empty()) {
                Face f = toVisit.front();
                toVisit.pop();

                for (Halfedge he : f.adjacentHalfedges()) {
                    Face g = he.twin().face();
                    if (!inPrimalSpanningTree(he.edge(), primalParent) &&
                        dualParent[g].face() == g && g != root) {
                        dualParent[g] = he;
                        toVisit.push(g);
                        inTree[g] = true;
                    }
                }
            }
        }
    }

    return dualParent;
}

VertexData<Halfedge> buildPrimalSpanningTree(ManifoldSurfaceMesh& mesh,
                                             const EdgeData<double>& edgeCost) {
    // Initialize each vertex as its own parent
    VertexData<Halfedge> primalParent(mesh);
    for (Vertex v : mesh.vertices()) primalParent[v] = v.halfedge();

    // Distance estimates to build a shortest-paths tree
    VertexData<double> dist(mesh, std::numeric_limits<double>::infinity());


    VertexData<bool> inTree(mesh, false);

    for (Vertex v : mesh.vertices()) {
        if (!inTree[v]) {
            inTree[v] = true;

            Vertex root = v;
            std::priority_queue<std::pair<double, Vertex>> toVisit;
            toVisit.push({0, root});
            dist[root] = 0;

            while (!toVisit.empty()) {
                Vertex v;
                std::tie(std::ignore, v) = toVisit.top();
                toVisit.pop();

                for (Halfedge he : v.outgoingHalfedges()) {
                    Vertex w = he.twin().vertex();
                    dist[w]  = fmin(dist[w], dist[v] + edgeCost[he.edge()]);
                    if (primalParent[w].vertex() == w && w != root) {
                        primalParent[w] = he;

                        // std::priority_queue is a max-heap, so we sort by
                        // negative cost to pop the closest vertex next
                        toVisit.push({-dist[w], w});
                        inTree[w] = true;
                    }
                }
            }
        }
    }

    return primalParent;
}

FaceData<Halfedge>
buildDualSpanningTree(ManifoldSurfaceMesh& mesh,
                      const VertexData<Halfedge>& primalParent,
                      const EdgeData<double>& edgeCost) {
    // Initialize each face as its own parent
    FaceData<Halfedge> dualParent(mesh);
    for (Face f : mesh.faces()) dualParent[f] = f.halfedge();

    // The cost of an edge not in the primal spanning tree is the length
    // of the shortest loop it makes with the primal spanning
    EdgeData<double> edgeLoopCost(mesh, 0);
    for (Edge e : mesh.edges()) {
        if (!inPrimalSpanningTree(e, primalParent)) {
            for (Edge e : primalTreeLoop(e, primalParent)) {
                edgeLoopCost[e] += edgeCost[e];
            }
        }
    }

    FaceData<bool> inTree(mesh, false);

    for (Face f : mesh.faces()) {
        if (!inTree[f]) {
            inTree[f] = true;

            Face root = f;
            std::priority_queue<std::pair<double, Face>> toVisit;
            toVisit.push({0, root});

            while (!toVisit.empty()) {
                Face f;
                std::tie(std::ignore, f) = toVisit.top();
                toVisit.pop();

                for (Halfedge he : f.adjacentHalfedges()) {
                    Face g = he.twin().face();
                    if (!inPrimalSpanningTree(he.edge(), primalParent) &&
                        dualParent[g].face() == g && g != root) {
                        dualParent[g] = he;
                        toVisit.push({edgeLoopCost[he.edge()], g});
                        inTree[g] = true;
                    }
                }
            }
        }
    }

    return dualParent;
}

std::vector<std::vector<Edge>>
computeHomologyGenerators(ManifoldSurfaceMesh& mesh) {

    VertexData<Halfedge> primalParent = buildPrimalSpanningTree(mesh);
    FaceData<Halfedge> dualParent = buildDualSpanningTree(mesh, primalParent);

    std::vector<std::vector<Edge>> generators;
    for (Edge e : mesh.edges()) {
        if (!inPrimalSpanningTree(e, primalParent) &&
            !inDualSpanningTree(e, dualParent)) {

            generators.push_back(primalTreeLoop(e, primalParent));
        }
    }

    return generators;
}

std::vector<std::vector<Edge>>
computeHomologyGenerators(ManifoldSurfaceMesh& mesh,
                          const EdgeData<double>& edgeCost) {

    VertexData<Halfedge> primalParent = buildPrimalSpanningTree(mesh, edgeCost);
    FaceData<Halfedge> dualParent =
        buildDualSpanningTree(mesh, primalParent, edgeCost);

    std::vector<std::vector<Edge>> generators;
    for (Edge e : mesh.edges()) {
        if (!inPrimalSpanningTree(e, primalParent) &&
            !inDualSpanningTree(e, dualParent)) {

            generators.push_back(primalTreeLoop(e, primalParent));
        }
    }

    return generators;
}

bool inPrimalSpanningTree(Edge e, const VertexData<Halfedge>& primalParent) {
    Vertex v = e.halfedge().vertex();
    Vertex w = e.halfedge().twin().vertex();
    return primalParent[v].edge() == e || primalParent[w].edge() == e;
}

bool inDualSpanningTree(Edge e, const FaceData<Halfedge>& dualParent) {
    Face f = e.halfedge().face();
    Face g = e.halfedge().twin().face();
    return dualParent[f].edge() == e || dualParent[g].edge() == e;
}

std::vector<Edge> primalTreeLoop(Edge e,
                                 const VertexData<Halfedge>& primalParent) {
    Halfedge he = e.halfedge();

    // Walk up along primal tree to extract generator
    Halfedge curr = primalParent[he.twin().vertex()];
    std::vector<Edge> forwardPath{he.edge(), curr.edge()};
    do {
        curr = primalParent[curr.vertex()];
        forwardPath.push_back(curr.edge());
    } while (primalParent[curr.vertex()].vertex() != curr.vertex());

    curr = primalParent[he.vertex()];
    std::vector<Edge> backwardPath{curr.edge()};
    do {
        curr = primalParent[curr.vertex()];
        backwardPath.push_back(curr.edge());
    } while (primalParent[curr.vertex()].vertex() != curr.vertex());

    // Identify and remove common edges at the ends of both paths
    size_t nShared = 0;
    size_t nF      = forwardPath.size() - 1;
    size_t nB      = backwardPath.size() - 1;
    while (nShared < nF && nShared < nB &&
           forwardPath[nF - nShared] == backwardPath[nB - nShared]) {
        nShared++;
    }

    // Remove last nShared elements from paths
    // https://stackoverflow.com/questions/34452139/how-to-remove-several-elements-from-the-end-of-stdvector
    forwardPath.resize(forwardPath.size() - nShared);
    backwardPath.resize(backwardPath.size() - nShared);


    std::vector<Edge> loop;
    // First go backwards along backwardPath
    loop.insert(std::end(loop), std::rbegin(backwardPath),
                std::rend(backwardPath));
    // Then go forwards along forwardPath
    loop.insert(std::end(loop), std::begin(forwardPath), std::end(forwardPath));
    return loop;
}
} // namespace ImplementationDetails
} // namespace CEPS
