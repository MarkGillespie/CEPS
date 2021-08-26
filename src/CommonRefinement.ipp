namespace CEPS {

namespace ImplementationDetails {

using std::cerr;
using std::cout;
using std::endl;


template <typename T1, typename T2>
std::tuple<SimplePolygonMesh, std::vector<T1>, std::vector<T2>, VertexData<int>>
computeCommonRefinement(Triangulation& Ta, Triangulation& Tb, Triangulation& Tc,
                        const VertexData<T1>& initialData,
                        const CornerData<T2>& finalData,
                        const std::vector<char>& frontFace, bool verbose) {

    auto isFrontFace = [&](size_t supportingFace) {
        return frontFace.empty() || frontFace[supportingFace];
    };

    std::vector<T1> vertices;
    std::vector<std::vector<size_t>> polygons;
    std::vector<T2> interpolatedData;
    VertexData<int> refinedVertices(*Ta.mesh);

    if (verbose) std::cout << "\t tracing Tc over Tb" << std::endl;
    EdgeData<MeshPath> edgesOnFinalTriangulation =
        traceTransposedGeodesicTriangulation(Tb, Tc, GeometryType::HYPERBOLIC,
                                             verbose);

    if (verbose) std::cout << "\t tracing Ta over Tb" << std::endl;
    EdgeData<MeshPath> edgesOnOriginalTriangulation =
        ImplementationDetails::transpose(
            traceGeodesicTriangulation(Ta, Tb, GeometryType::EUCLIDEAN,
                                       verbose),
            Ta, Tb);

    // Make CornerData on original mesh for interpolation
    CornerData<DirectSum<T1, double>> posAndParentFace(*Ta.mesh);

    FaceData<size_t> fIdx = Ta.mesh->getFaceIndices();
    for (Corner c : Ta.mesh->corners()) {
        posAndParentFace[c] = DirectSum<T1, double>(initialData[c.vertex()],
                                                    (double)fIdx[c.face()]);
    }

    VertexData<int> vertexIndices(*Tb.mesh, -1);
    EdgeData<std::vector<size_t>> edgeVertices(*Tb.mesh);

    auto nEdgePoints = [&](Edge e) {
        return edgesOnFinalTriangulation[e].points.size() +
               edgesOnOriginalTriangulation[e].points.size();
    };
    std::vector<size_t> parentFaces;

    VertexData<size_t> vIdx = Tb.mesh->getVertexIndices();
    VertexData<Vertex> supportingToOriginal(*Tb.mesh);
    for (Vertex v : Tb.mesh->vertices())
        supportingToOriginal[v] = Ta.mesh->vertex(vIdx[v]);

    size_t nF = Tb.mesh->nFaces();
    size_t iF = 0;
    if (verbose) std::cout << "\t subdividing faces";
    for (Face f : Tb.mesh->faces()) {
        std::array<std::vector<double>, 3> originalBary, finalBary;

        std::array<std::vector<std::array<DirectSum<T1, double>, 2>>, 3>
            originalCornerData;
        std::array<std::vector<std::array<T2, 2>>, 3> finalCornerData;

        std::vector<DirectSum<T1, double>> longSideOppositeOriginalData;
        std::vector<T2> longSideOppositeCornerData;

        std::tie(finalBary, finalCornerData, longSideOppositeCornerData) =
            readEdgeCrossings(f, edgesOnFinalTriangulation, finalData, true,
                              Tc.logScaleFactors);

        std::tie(originalBary, originalCornerData,
                 longSideOppositeOriginalData) =
            readEdgeCrossings(f, edgesOnOriginalTriangulation, posAndParentFace,
                              false, Ta.logScaleFactors);

        size_t nNewVertices;
        std::vector<std::vector<size_t>> slicedFaces;
        std::vector<std::vector<DirectSum<T1, double>>> newCornerPositions;
        std::vector<std::vector<T2>> newCornerData;


        bool verboseSlicing = false;
        std::tie(nNewVertices, slicedFaces, newCornerPositions, newCornerData) =
            sliceTri(originalBary, finalBary, originalCornerData,
                     finalCornerData, longSideOppositeOriginalData,
                     longSideOppositeCornerData, verboseSlicing);

        std::vector<size_t> newVertexIndices;
        size_t firstNewVertex = vertices.size();
        size_t iV             = firstNewVertex;
        bool glueVertices     = true;

        if (glueVertices) {
            // Vertices are returned from sliceTri in this order
            // (ordered by the opposite halfedge)
            std::array<Vertex, 3> faceVertices{
                f.halfedge().next().next().vertex(), f.halfedge().vertex(),
                f.halfedge().next().vertex()};

            for (Vertex v : faceVertices) {
                if (vertexIndices[v] < 0) vertexIndices[v] = (int)iV++;
                newVertexIndices.push_back((size_t)vertexIndices[v]);
                refinedVertices[supportingToOriginal[v]] = vertexIndices[v];
            }

            size_t nBoundaryVertices = 3;
            for (Halfedge he : f.adjacentHalfedges()) {
                if (edgeVertices[he.edge()].size() > 0) {
                    const std::vector<size_t>& newEdgeVs =
                        edgeVertices[he.edge()];

                    // Use reversed list since it came from the twin
                    // halfedge
                    newVertexIndices.insert(std::end(newVertexIndices),
                                            std::rbegin(newEdgeVs),
                                            std::rend(newEdgeVs));
                } else {
                    std::vector<size_t> newEdgeVs;
                    for (size_t iPt = 0; iPt < nEdgePoints(he.edge()); iPt++) {
                        newEdgeVs.push_back(iV++);
                    }
                    newVertexIndices.insert(std::end(newVertexIndices),
                                            std::begin(newEdgeVs),
                                            std::end(newEdgeVs));
                    edgeVertices[he.edge()] = newEdgeVs;
                }
                nBoundaryVertices += nEdgePoints(he.edge());
            }

            for (size_t j = 0; j + nBoundaryVertices < nNewVertices; ++j) {
                newVertexIndices.push_back(iV++);
            }
        } else {

            // Vertices are returned from sliceTri in this order
            // (ordered by the opposite halfedge)
            std::array<Vertex, 3> faceVertices{
                f.halfedge().next().next().vertex(), f.halfedge().vertex(),
                f.halfedge().next().vertex()};

            for (Vertex v : faceVertices) {
                refinedVertices[supportingToOriginal[v]] = iV;
                newVertexIndices.push_back(iV++);
            }

            for (size_t j = 3; j < nNewVertices; ++j) {
                newVertexIndices.push_back(iV++);
            }
        }

        std::vector<T1> newVertexPositions(iV - firstNewVertex);

        // Extract list of new vertex positions from new corner
        // positions Also, re-index returned faces
        bool faceParentsOkay = true;
        for (size_t iP = 0; iP < slicedFaces.size(); ++iP) {
            int myFaceParentInt = std::round(newCornerPositions[iP][0].y);
            if (myFaceParentInt < 0) {
                cerr << "face parent is negative? " << myFaceParentInt << endl;
                cerr << "                         "
                     << newCornerPositions[iP][0].y << endl;
            }
            size_t myFaceParent = std::max(0, myFaceParentInt);
            if (isFrontFace(myFaceParent)) {

                bool printed = false;
                for (size_t iC = 0; iC < slicedFaces[iP].size(); ++iC) {
                    size_t iV = slicedFaces[iP][iC];
                    size_t cornerFaceParent =
                        std::round(newCornerPositions[iP][iC].y);
                    faceParentsOkay =
                        faceParentsOkay && (cornerFaceParent == myFaceParent);

                    if (newVertexIndices[iV] >= firstNewVertex) {
                        size_t idx = newVertexIndices[iV] - firstNewVertex;
                        // Record the position of this corner
                        newVertexPositions[idx] = newCornerPositions[iP][iC].x;
                    }

                    // Shift the indices of the vertices in the face
                    // by the number of existing vertices
                    slicedFaces[iP][iC] = newVertexIndices[iV];
                }

                parentFaces.push_back(myFaceParent);
                polygons.push_back(slicedFaces[iP]);
                interpolatedData.insert(std::end(interpolatedData),
                                        std::begin(newCornerData[iP]),
                                        std::end(newCornerData[iP]));
            }
        }
        // TODO: add checks back in?
        /*
        my_assert(faceParentsOkay, "face parents on refined mesh
        inconsistent " "between corners");

        my_assert(nNewVertices >= 3, "Not enough new vertices");
        my_assert(slicedFaces.size() >= 1, "Not enough new faces");
        */

        // add vertices, faces, and interpolated data to appropriate
        // lists
        vertices.insert(std::end(vertices), std::begin(newVertexPositions),
                        std::end(newVertexPositions));
        iF++;
        if (verbose)
            std::cout << "\r\t subdividing faces " << iF << " / " << nF;
    }
    if (verbose) std::cout << std::endl;

    std::vector<Vector3> dummyVertexPositions(vertices.size());
    SimplePolygonMesh soup(polygons, dummyVertexPositions);

    return std::make_tuple(soup, vertices, interpolatedData, refinedVertices);
}


template <typename T>
std::tuple<std::array<std::vector<double>, 3>,
           std::array<std::vector<std::array<T, 2>>, 3>, std::vector<T>>
readEdgeCrossings(Face f, const EdgeData<MeshPath>& paths,
                  const CornerData<T>& inputData, bool projectiveInterpolate,
                  const VertexData<double>& logScaleFactors, bool verbose) {

    std::array<std::vector<double>, 3> bary;
    std::array<std::vector<std::array<T, 2>>, 3> data;
    std::vector<T> oppositeCornerData;

    std::array<Halfedge, 3> fHalfedges{f.halfedge(), f.halfedge().next(),
                                       f.halfedge().next().next()};
    std::array<size_t, 3> nCrossings;
    for (size_t iE = 0; iE < 3; ++iE) {
        nCrossings[iE] = paths[fHalfedges[iE].edge()].points.size();
    }

    size_t fanEdge = getFanEdge(nCrossings);

    Vertex fanVertex;
    if (fanEdge < 3) {
        fanVertex = fHalfedges[fanEdge].next().next().vertex();
    }

    for (size_t iE = 0; iE < 3; ++iE) {
        Halfedge he = fHalfedges[iE];
        Edge e      = he.edge();

        MeshPath crossings = paths[e];

        // Make sure MeshPath points along halfedge orientation
        if (he != e.halfedge()) {
            crossings = twin(crossings);
        }

        if (crossings.points.empty()) {
            // Special cases for MeshPath's which don't cross any edges
            // of the original mesh. Note that such an edge cannot be a
            // fan edge

            Vertex vStart = crossings.start.corner.vertex();
            Vertex vEnd   = crossings.end.corner.vertex();

            Halfedge meshHe = crossings.start.corner.halfedge();
            if (meshHe.next().corner() != crossings.end.corner) {
                meshHe = crossings.end.corner.halfedge();
                verbose_assert(meshHe.next().corner() == crossings.start.corner,
                               "start and end corners are not compatible");
                meshHe = meshHe.twin();
            }

            verbose_assert(meshHe.vertex() == vStart, "Actually impossible");
            verbose_assert(meshHe.twin().vertex() == vEnd,
                           "Definitely shouldn't happen");

            double startBary = crossings.start.bary;
            double endBary   = crossings.end.bary;

            T startData = inputData[meshHe.corner()] * startBary;
            T endData   = inputData[meshHe.next().corner()] * endBary;

            std::array<T, 2> startPair{startData, startData};
            std::array<T, 2> endPair{endData, endData};

            bary[iE] = std::vector<double>{0, 1};
            data[iE] = std::vector<std::array<T, 2>>{startPair, endPair};

        } else {
            // Start point
            CornerPt start = crossings.start;
            T startData    = start.bary * inputData[start.corner];
            bary[iE].push_back(0);
            data[iE].push_back(std::array<T, 2>{startData, startData});

            // Intermediate HalfedgePts
            for (size_t i = 0; i < crossings.points.size(); ++i) {
                double t = crossings.baryCoords[i];
                if (t != t) {
                    cerr << "Tracing Error: t = " << t
                         << "\t|\t path length = " << crossings.points.size()
                         << endl;
                }
                if (t != t) {
                    t = 0.5;
                } else {
                    t = clamp(t, 0., 1.);
                }

                HalfedgePt hePt  = crossings.points[i];
                Halfedge frontHe = hePt.halfedge;
                Vector2 bc       = hePt.bary;

                T heSrcData = inputData[frontHe.corner()];
                T heDstData = inputData[frontHe.next().corner()];

                T frontVal = bc.x * heSrcData + bc.y * heDstData;

                // Src and Dst are still relative to frontHe's
                // orientation
                T heOppSrcData = inputData[frontHe.twin().next().corner()];
                T heOppDstData = inputData[frontHe.twin().corner()];

                T backVal = bc.x * heOppSrcData + bc.y * heOppDstData;

                bary[iE].push_back(t);
                data[iE].push_back(std::array<T, 2>{frontVal, backVal});
                if (iE == fanEdge) {
                    size_t fanStart = nCrossings[prev(iE)];
                    size_t fanEnd   = nCrossings[iE] - nCrossings[next(iE)] - 1;

                    double vScale =
                        projectiveInterpolate
                            ? exp(-logScaleFactors[frontHe.next().vertex()])
                            : 1;

                    if (i >= fanStart && i <= fanEnd) {
                        oppositeCornerData.push_back(
                            inputData[frontHe.next().corner()] * vScale);
                    }
                    if (i == fanEnd) {
                        oppositeCornerData.push_back(
                            inputData[frontHe.twin().corner()] * vScale);
                    }
                }
            }

            // end point
            CornerPt end = crossings.end;
            T endData    = end.bary * inputData[end.corner];
            bary[iE].push_back(1);
            data[iE].push_back(std::array<T, 2>{endData, endData});
        }
    }

    // Match up data on line endpoints properly
    for (size_t iE = 0; iE < 3; ++iE) {
        data[iE][0][0] = data[prev(iE)][data[prev(iE)].size() - 1][0];
        data[iE][data[iE].size() - 1][1] = data[next(iE)][0][1];
    }

    return std::make_tuple(bary, data, oppositeCornerData);
}

// Returns: the number of vertices in the refined triangle, the face-vertex
// adjacency list of the triangle, and the interpolated values of val1 and val2
// within each face
template <typename T1, typename T2>
std::tuple<size_t, std::vector<std::vector<size_t>>,
           std::vector<std::vector<T1>>, std::vector<std::vector<T2>>>
sliceTri(std::array<std::vector<double>, 3> bary1,
         std::array<std::vector<double>, 3> bary2,
         const std::array<std::vector<std::array<T1, 2>>, 3>& val1,
         const std::array<std::vector<std::array<T2, 2>>, 3>& val2,
         const std::vector<T1>& longSideCornerVal1,
         const std::vector<T2>& longSideCornerVal2, bool verbose) {

    if (verbose) {
        cout << "bary1: " << endl;
        for (size_t iE = 0; iE < 3; ++iE) {
            cout << "\t edge " << iE << " : ";
            for (double t : bary1[iE]) cout << " " << t;
            cout << endl;
        }
        cout << "bary2: " << endl;
        for (size_t iE = 0; iE < 3; ++iE) {
            cout << "\t edge " << iE << " : ";
            for (double t : bary2[iE]) cout << " " << t;
            cout << endl;
        }
    }

    //== Early exit if no slicing is required
    if (bary1[0].size() == 2 && bary1[1].size() == 2 && bary1[2].size() == 2 &&
        bary2[0].size() == 2 && bary2[1].size() == 2 && bary2[2].size() == 2) {
        size_t nVertices = 3;
        std::vector<std::vector<size_t>> faces{std::vector<size_t>{1, 2, 0}};
        std::vector<std::vector<T1>> val1Interp{
            {val1[0][0][0], val1[1][0][0], val1[2][0][0]}};
        std::vector<std::vector<T2>> val2Interp{
            {val2[0][0][0], val2[1][0][0], val2[2][0][0]}};

        return std::make_tuple(nVertices, faces, val1Interp, val2Interp);
    }

    // Check compatibility of values at corners (value at end of line should
    // match value at start of next line)
    for (size_t iE = 0; iE < 3; ++iE) {
        size_t N = bary1[iE].size();
        verbose_assert(bary1[iE][0] == 0,
                       "Edges must start with a barycentric coordiante of 0");
        verbose_assert(bary1[iE][N - 1] == 1,
                       "Edges must end with a barycentric coordiante of 1");
        verbose_assert(val1[iE].size() == N,
                       "value array must contain one pair per element of the "
                       "barycentric coordinate array");
        verbose_assert(
            val1[iE][N - 1] == val1[next(iE)][0],
            "Value pair at end of one line must match value pair at the "
            "start of the next line");

        N = bary2[iE].size();
        verbose_assert(bary2[iE][0] == 0,
                       "Edges must start with a barycentric coordiante of 0");
        verbose_assert(bary2[iE][N - 1] == 1,
                       "Edges must end with a barycentric coordiante of 1");
        verbose_assert(val2[iE].size() == N,
                       "value array must contain one pair per element of the "
                       "barycentric coordinate array");

        if (!(val2[iE][N - 1] == val2[next(iE)][0])) {
            cerr << "Bad values: " << val2[iE][N - 1][0] << ", "
                 << val2[iE][N - 1][1] << " vs " << val2[next(iE)][0][0] << ", "
                 << val2[next(iE)][0][1] << endl;
        }
        // verbose_assert(!(val2[iE][N - 1] != val2[next(iE)][0]),
        //           "Value pair at end of one line must match value pair at
        //           the" "start of the next line");
        verbose_assert(
            val2[iE][N - 1] == val2[next(iE)][0],
            "Value pair at end of one line must match value pair at the "
            "start of the next line");
    }

    // For convenience, lay out a concrete triangle. It doesn't matter what
    // this triangle is. It should be oriented counterclockwise
    std::array<Vector2, 3> oppVertices{Vector2{0, 0}, Vector2{1, 0},
                                       Vector2{0, 1}};

    std::vector<std::vector<size_t>> boundaryVertices(3);

    // Initialize vertex list with triangle vertices and add them
    // to the list of vertices on each boundary edge
    std::vector<Vector2> vertices{oppVertices[0], oppVertices[1],
                                  oppVertices[2]};
    // Add starting vertex to boundary list. We'll add the final vertex
    // after adding the intermediate ones
    for (size_t iE = 0; iE < 3; ++iE) {
        boundaryVertices[prev(iE)].push_back(iE);
    }

    // arrays of interpolated vertex data
    std::vector<std::array<T1, 2>> vertexValues1{
        val1[prev(0)][0], val1[prev(1)][0], val1[prev(2)][0]};
    std::vector<std::array<T2, 2>> vertexValues2{
        val2[prev(0)][0], val2[prev(1)][0], val2[prev(2)][0]};

    // arrays to tell if we've filled in data yet
    std::vector<char> hasValue1{true, true, true};
    std::vector<char> hasValue2{true, true, true};

    // Used for looping around boundary vertices in order to
    // interpolate values
    std::array<std::vector<std::pair<double, size_t>>, 3> boundaryVertexTimes;
    for (size_t iE = 0; iE < 3; ++iE) {
        // put first vertex at the beginning of the list. We'll add the last
        // vertex after populating the middle of the list
        boundaryVertexTimes[iE].push_back(std::make_pair(0, next(iE)));
        // boundaryVertexTimes[iE].push_back(std::make_pair(1, prev(iE)));
    }

    // List of which vertices lie on which edge - used to compute which
    // vertices are connected by lines
    std::array<std::vector<size_t>, 3> edgeVertices1, edgeVertices2;

    // Records what each boundary vertex's position is in its edge (counting
    // corner vertices)
    // Just set it to 0 for corners.
    std::vector<size_t> indexInEdge{0, 0, 0};

    // Read off data from triangle boundary edges
    for (size_t iE = 0; iE < 3; ++iE) {
        // function to compute a point from its barycentric
        // coordinates
        auto bc = [&](double t) {
            // next(iE) is the source of edge iE
            // prev(iE) is the destination of edge iE
            return (1 - t) * oppVertices[next(iE)] + t * oppVertices[prev(iE)];
        };

        if (bary1[iE].size() > 2 || bary2[iE].size() > 2) {
            size_t iA = 1, iB = 1;
            while (iA < bary1[iE].size() - 1 || iB < bary2[iE].size() - 1) {
                if (iB + 1 >= bary2[iE].size() ||
                    bary1[iE][iA] < bary2[iE][iB]) {
                    double t = bary1[iE][iA];
                    boundaryVertices[iE].push_back(vertices.size());
                    edgeVertices1[iE].push_back(vertices.size());
                    boundaryVertexTimes[iE].push_back(
                        std::make_pair(t, vertices.size()));

                    vertices.push_back(bc(t));
                    vertexValues1.push_back(val1[iE][iA]);
                    vertexValues2.push_back(std::array<T2, 2>());
                    hasValue1.push_back(true);
                    hasValue2.push_back(false);
                    indexInEdge.push_back(iA - 1);

                    iA++;
                } else {
                    double t = bary2[iE][iB];
                    boundaryVertices[iE].push_back(vertices.size());
                    edgeVertices2[iE].push_back(vertices.size());
                    boundaryVertexTimes[iE].push_back(
                        std::make_pair(t, vertices.size()));

                    vertices.push_back(bc(t));
                    vertexValues1.push_back(std::array<T1, 2>());
                    vertexValues2.push_back(val2[iE][iB]);
                    hasValue1.push_back(false);
                    hasValue2.push_back(true);
                    indexInEdge.push_back(iB - 1);

                    iB++;
                }
            }
        }
    }

    // Add in the last vertex to each edge
    for (size_t iE = 0; iE < 3; ++iE) {
        boundaryVertices[next(iE)].push_back(iE);
        boundaryVertexTimes[iE].push_back(std::make_pair(1, prev(iE)));
    }

    // ===========================================================
    //                Identify Fan Edges
    // ===========================================================

    // Identify which, if any, of the triangle's vertices have edges coming
    // from them, in lines1 and lines2. Note that since we index edges by
    // their opposite vertex, this is equivalent to finding a "fan edge"
    // TODO: this computation is copied in getLines
    size_t lines1FanEdge  = getFanEdge(edgeVertices1);
    size_t lines2FanEdge  = getFanEdge(edgeVertices2);
    auto isLine1FanVertex = [&](size_t iV) -> bool {
        return lines1FanEdge < 3 && iV == lines1FanEdge;
    };
    auto isLine2FanVertex = [&](size_t iV) -> bool {
        return lines2FanEdge < 3 && iV == lines2FanEdge;
    };
    bool hasSharedFanVertex =
        (lines1FanEdge < 3) && (lines1FanEdge == lines2FanEdge);

    auto incidentOnSharedFanVertex = [&](Line l) -> bool {
        return hasSharedFanVertex &&
               (l[0] == lines1FanEdge || l[1] == lines1FanEdge);
    };

    size_t firstLine1FanCrossing = edgeVertices1[prev(lines1FanEdge)].size();
    size_t firstLine2FanCrossing = edgeVertices2[prev(lines2FanEdge)].size();

    // ===========================================================
    //                Interpolate Boundary Data
    // ===========================================================

    for (size_t iE = 0; iE < 3; ++iE) {
        // std::sort(std::begin(boundaryVertexTimes[iE]),
        //           std::end(boundaryVertexTimes[iE]));
        const std::vector<std::pair<double, size_t>>& vertexTimes =
            boundaryVertexTimes[iE];
        // std::sort(std::begin(vertexTimes), std::end(vertexTimes));

        auto vtx  = [&](size_t i) { return vertexTimes[i].second; };
        auto time = [&](size_t i) { return vertexTimes[i].first; };
        auto findNextVal1Vertex = [&](size_t start) {
            size_t next = start + 1;
            while (!hasValue1[vtx(next)]) {
                next++;
                if (next >= vertexTimes.size()) {
                    cerr << "Ran off the end of vertexTimes (1)" << endl;
                }
            }
            return next;
        };
        auto findNextVal2Vertex = [&](size_t start) {
            size_t next = start + 1;
            while (!hasValue2[vtx(next)]) {
                next++;
                if (next >= vertexTimes.size()) {
                    cerr << "Ran off the end of vertexTimes (2)" << endl;
                }
            }
            return next;
        };

        // interpolate value1
        size_t lastValueVertex = 0;
        double lastT           = time(lastValueVertex);

        while (lastValueVertex < vertexTimes.size() - 1) {
            size_t nextValueVertex = findNextVal1Vertex(lastValueVertex);
            T1 lastValue           = vertexValues1[vtx(lastValueVertex)][1];
            T1 nextValue           = vertexValues1[vtx(nextValueVertex)][0];
            double nextT           = time(nextValueVertex);
            double deltaT          = nextT - lastT;
            for (size_t i = lastValueVertex + 1; i < nextValueVertex; ++i) {
                double t = (time(i) - lastT) / deltaT;
                if (t != t) {
                    cerr << "Error: t is NaN. setting to 0.5" << endl;
                    cerr << "        lastT = " << lastT << endl;
                    cerr << "        time(i) = " << time(i) << endl;
                    cerr << "        deltaT = " << deltaT << endl;
                    t = 0.5;
                }
                T1 val = (1 - t) * lastValue + t * nextValue;
                if (val != val) {
                    cerr << "Error: value1 is NaN" << endl;
                    cerr << "        t = " << t << endl;
                    cerr << "lastValue = " << lastValue << endl;
                    cerr << "nextValue = " << nextValue << endl;
                }

                if (t < 0 || t > 1) {
                    cerr << "invalid t = " << t << endl;
                    t = 0.5;
                }

                vertexValues1[vtx(i)] = std::array<T1, 2>{val, val};
                hasValue1[vtx(i)]     = true;

                // if (norm(getX(val)) < 1e-12) {
                //     cerr << "val is 0" << endl;
                // }

                if (verbose) {
                    cout << "interpolated val1: " << val << endl;
                }
            }
            lastValueVertex = nextValueVertex;
            lastT           = nextT;
        }

        // interpolate value2
        lastValueVertex = 0;
        lastT           = time(lastValueVertex);
        while (lastValueVertex < vertexTimes.size() - 1) {
            size_t nextValueVertex = findNextVal2Vertex(lastValueVertex);
            T2 lastValue           = vertexValues2[vtx(lastValueVertex)][1];
            T2 nextValue           = vertexValues2[vtx(nextValueVertex)][0];
            double nextT           = time(nextValueVertex);
            double deltaT          = nextT - lastT;
            for (size_t i = lastValueVertex + 1; i < nextValueVertex; ++i) {
                double t = (time(i) - lastT) / deltaT;
                if (t != t) {
                    cerr << "Error: t is NaN. setting to 0.5" << endl;
                    cerr << "        lastT = " << lastT << endl;
                    cerr << "        time(i) = " << time(i) << endl;
                    cerr << "        deltaT = " << deltaT << endl;
                    t = 0.5;
                }
                T2 val                = (1 - t) * lastValue + t * nextValue;
                vertexValues2[vtx(i)] = std::array<T2, 2>{val, val};
                if (verbose) {
                    cout << "interpolated val2: " << val << endl;
                }
            }
            lastValueVertex = nextValueVertex;
            lastT           = nextT;
        }
    }

    // Check that we've interpolated data to all boundary vertices
    for (bool b : hasValue1) {
        if (!b) {
            cerr << "Error: missed a vertex when interpolating positions"
                 << endl;
        }
    }

    // ===========================================================
    //                      Compute Lines
    // ===========================================================

    // For now, I keep the lists of lines in two vectors because
    // I know that each set of lines doesn't intersect itself.
    // After I compute the intersections, I merge the lists
    // together
    std::vector<Line> lines1 = getLines(edgeVertices1);
    std::vector<Line> lines2 = getLines(edgeVertices2);
    std::vector<std::vector<size_t>> line1Vertices(lines1.size());
    std::vector<std::vector<size_t>> line2Vertices(lines2.size());

    // Record that lines contain their endpoints
    for (size_t iL = 0; iL < lines1.size(); ++iL) {
        line1Vertices[iL].push_back(lines1[iL][0]);
        line1Vertices[iL].push_back(lines1[iL][1]);
    }

    for (size_t iL = 0; iL < lines2.size(); ++iL) {
        line2Vertices[iL].push_back(lines2[iL][0]);
        line2Vertices[iL].push_back(lines2[iL][1]);
    }


    // ===========================================================
    //                Compute Intersection Vertices
    // ===========================================================
    // Record the first vertex along each boundary edge
    std::array<size_t, 3> edgeLen{3, 3, 3};
    edgeLen[0] = bary1[0].size() - 2 + bary2[0].size() - 2;
    edgeLen[1] = bary1[1].size() - 2 + bary2[1].size() - 2;
    edgeLen[2] = bary1[2].size() - 2 + bary2[2].size() - 2;

    auto cyclicIndex = [&](size_t iV) -> size_t {
        if (iV == 1) {
            return 0;
        } else if (iV == 2) {
            return 1 + edgeLen[0];
        } else if (iV == 0) {
            return 2 + edgeLen[0] + edgeLen[1];
        } else if (iV < 3 + edgeLen[0]) {
            return 1 + (iV - 3);
        } else if (iV < 3 + edgeLen[0] + edgeLen[1]) {
            return 2 + (iV - 3);
        } else {
            return 3 + (iV - 3);
        }
    };

    std::vector<size_t> cyclicOrdering;
    for (size_t iV = 0; iV < vertices.size(); ++iV) {
        cyclicOrdering.push_back(cyclicIndex(iV));
    }

    auto cyclicIntersect = [&](Line a, Line b,
                               const std::vector<Vector2>& position, double& tA,
                               double& tB) -> bool {
        size_t a0 = cyclicOrdering[a[0]];
        size_t a1 = cyclicOrdering[a[1]];
        if (a0 > a1) std::swap(a0, a1);
        size_t b0 = cyclicOrdering[b[0]];
        size_t b1 = cyclicOrdering[b[1]];
        if (b0 > b1) std::swap(b0, b1);

        // return intersect(a, b, position, tA, tB);
        if ((a0 <= b0 && b0 <= a1 && a1 <= b1) ||
            (b0 <= a0 && a0 <= b1 && b1 <= a1)) {
            // compute the intersection, using the knowledge that they must
            // intersect
            unsafeIntersect(a, b, position, tA, tB);
            return true;
        } else {
            return false;
        }
    };

    std::vector<std::vector<size_t>> line1Crosses(lines1.size()),
        line2Crosses(lines2.size());

    for (size_t iA = 0; iA < lines1.size(); ++iA) {
        Line a = lines1[iA];
        for (size_t iB = 0; iB < lines2.size(); ++iB) {
            Line b = lines2[iB];
            double tA, tB;

            // Ignore intersections at fan edges
            if (incidentOnSharedFanVertex(a) && incidentOnSharedFanVertex(b)) {
                continue;
            }

            if (cyclicIntersect(a, b, vertices, tA, tB)) {
                if (verbose) {
                    cout << "INTERSECTION" << endl;
                }
                Vector2 intersection =
                    (1 - tA) * vertices[a[0]] + tA * vertices[a[1]];
                line1Vertices[iA].push_back(vertices.size());
                line2Vertices[iB].push_back(vertices.size());
                line1Crosses[iA].push_back(iB);
                line2Crosses[iB].push_back(iA);
                vertices.push_back(intersection);

                std::array<T1, 2> startVal1 = vertexValues1[a[0]];
                std::array<T1, 2> endVal1   = vertexValues1[a[1]];
                // If a starts at the fan vertex, we have to do something a
                // bit more complicated to find the right values at the
                // vertex Note that only the starting vertex can be the fan
                // vertex due to how we construct these lines
                if (isLine1FanVertex(a[0])) {
                    // We need to find the position of vertex a[1] in its
                    // edge. Since line a goes through a triangle vertex, we
                    // know that its other endpoint must lie on
                    // lines1FanEdge. We build the vertex list by first
                    // adding in the triangle vertices, then adding the
                    // vertices from each edge in order. So to find the
                    // index of a[1] in its edge, we just take a[1]'s index,
                    // subtract 3, and subtract the sizes of the previous
                    // edges
                    // size_t iEnd = a[1] - 3;
                    // for (size_t iE = 0; iE < lines1FanEdge; ++iE) {
                    //     iEnd -= edgeVertices1[iE].size();
                    //     iEnd -= edgeVertices2[iE].size();
                    // }
                    // iEnd -= firstLine1FanCrossing;
                    size_t iEnd = indexInEdge[a[1]] - firstLine1FanCrossing;

                    // TODO: right now I recompute these values for the
                    // segment at the end of line a later on. Could I store
                    // them somehow?
                    startVal1 = std::array<T1, 2>{longSideCornerVal1[iEnd + 1],
                                                  longSideCornerVal1[iEnd]};
                }

                std::array<T2, 2> startVal2 = vertexValues2[b[0]];
                std::array<T2, 2> endVal2   = vertexValues2[b[1]];
                if (isLine2FanVertex(b[0])) {
                    // size_t iEnd = b[1] - 3;
                    // for (size_t iE = 0; iE < lines2FanEdge; ++iE) {
                    //     iEnd -= edgeVertices1[iE].size();
                    //     iEnd -= edgeVertices2[iE].size();
                    // }
                    // iEnd -= edgeVertices1[lines2FanEdge].size();
                    // iEnd -= firstLine2FanCrossing;

                    size_t iEnd = indexInEdge[b[1]] - firstLine2FanCrossing;

                    startVal2 = std::array<T2, 2>{longSideCornerVal2[iEnd + 1],
                                                  longSideCornerVal2[iEnd]};
                }

                // Interpolate values to intersection point
                T1 value1Left  = (1 - tA) * startVal1[0] + tA * endVal1[1];
                T1 value1Right = (1 - tA) * startVal1[1] + tA * endVal1[0];

                T2 value2Left  = (1 - tB) * startVal2[0] + tB * endVal2[1];
                T2 value2Right = (1 - tB) * startVal2[1] + tB * endVal2[0];

                if (value1Left != value1Left) {
                    cerr << "Error: value1 is NaN" << endl;
                    cerr << "   tA = " << tA << endl;
                    cerr << "start = " << startVal1[0] << endl;
                    cerr << "  end = " << endVal1[1] << endl;
                    exit(1);
                } else if (value1Right != value1Right) {
                    cerr << "Error: value1 is NaN" << endl;
                    cerr << "   tA = " << tA << endl;
                    cerr << "start = " << startVal1[1] << endl;
                    cerr << "  end = " << endVal1[0] << endl;
                    exit(1);
                }


                // if (norm(getX(value1Left)) < 1e-12) {
                //     cerr << "value1Left is 0" << endl;
                // } else if (norm(getX(value1Right)) < 1e-12) {
                //     cerr << "value1Right is 0" << endl;
                // }

                vertexValues1.push_back(
                    std::array<T1, 2>{value1Left, value1Right});
                vertexValues2.push_back(
                    std::array<T2, 2>{value2Left, value2Right});
            }
        }
    }

    // Re-order data at boundary vertices so that it obeys our
    // right-hand-rule (i.e. value[0] is to the left of the relevant line,
    // value[1] is to the right) This is true at the line's start, but not
    // at its end, since right now value[0] is always earlier on the
    // boundary edge. So we swap the line's endpoint values. Note that for
    // line1 lines we only have to swap val1 and likewise for line2
    for (Line a : lines1) {
        std::swap(vertexValues1[a[1]][0], vertexValues1[a[1]][1]);
    }
    for (Line b : lines2) {
        std::swap(vertexValues2[b[1]][0], vertexValues2[b[1]][1]);
    }

    auto orderIntersections = [](const Line& line,
                                 std::vector<size_t>& lineVertices,
                                 const std::vector<size_t> lineCrossings,
                                 const std::vector<Line>& oppLines,
                                 const std::vector<size_t>& cyclicOrdering) {
        if (lineCrossings.empty()) return;

        size_t N      = cyclicOrdering.size();
        size_t iStart = cyclicOrdering[line[0]];

        auto cwiseIdx = [=](size_t iV) -> size_t {
            return (cyclicOrdering[iV] + N - iStart) % N;
        };
        auto counterwiseIdx = [=](size_t iV) -> size_t {
            return (iStart + N - cyclicOrdering[iV]) % N;
        };

        verbose_assert(lineCrossings.size() + 2 == lineVertices.size(),
                       "size err");
        std::vector<std::pair<std::pair<size_t, size_t>, size_t>>
            intersectionTimes;
        for (size_t iI = 0; iI < lineCrossings.size(); ++iI) {
            size_t iB      = lineCrossings[iI];
            Line lineCross = oppLines[iB];

            size_t iC =
                std::min(cwiseIdx(lineCross[0]), cwiseIdx(lineCross[1]));
            size_t iCC = std::min(counterwiseIdx(lineCross[0]),
                                  counterwiseIdx(lineCross[1]));
            intersectionTimes.push_back(
                std::make_pair(std::make_pair(iC, iCC), lineVertices[iI + 2]));
            // +2 because lineVertices starts with the line's start and
            // endpoints, which I don't consider here
        }

        // std::sort sorts pairs lexicographically, so this
        // effectively sorts by cyclic ordering coordinates
        // We sort by both the clockwise and counterclockwise distance.
        // Either one works, except in the case of a pair of fan edges
        // which share an endpoint. Since two distinct lines must differ
        // in at least one endpoint, we sort lexicographically by both
        // distances
        std::sort(std::begin(intersectionTimes), std::end(intersectionTimes));

        // Move endpoint to end
        lineVertices[lineVertices.size() - 1] = lineVertices[1];
        // Copy over intermediate intersections
        for (size_t iV = 0; iV < intersectionTimes.size(); ++iV) {
            lineVertices[iV + 1] = intersectionTimes[iV].second;
        }
    };

    // Sort intersections combinatorially
    for (size_t iA = 0; iA < lines1.size(); ++iA) {
        orderIntersections(lines1[iA], line1Vertices[iA], line1Crosses[iA],
                           lines2, cyclicOrdering);
    }
    for (size_t iB = 0; iB < lines2.size(); ++iB) {
        orderIntersections(lines2[iB], line2Vertices[iB], line2Crosses[iB],
                           lines1, cyclicOrdering);
    }

    // ===========================================================
    //     Split lines up into line segments
    //           (cut at intersections)
    // ===========================================================
    // List of all line segments in triangle
    std::vector<Segment> segments;
    // List of segments incident on each vertex
    std::vector<std::vector<size_t>> vertexSegments(vertices.size());
    std::vector<std::vector<size_t>> vertexSegmentDirections(vertices.size());
    // Store whether each segment is part of a line from a, a line from
    // b, or a boundary edge
    std::vector<LineType> segmentTypes;

    // Split interior lines
    splitIntoSegments(lines1, line1Vertices, vertices, LineType::A,
                      vertexSegments, vertexSegmentDirections, segments,
                      segmentTypes, true);
    splitIntoSegments(lines2, line2Vertices, vertices, LineType::B,
                      vertexSegments, vertexSegmentDirections, segments,
                      segmentTypes, true);

    // check orientation of edges and segments incident on fan vertices
    if (lines1FanEdge < 3) {
        for (Line a : lines1) {
            if (isLine1FanVertex(a[0]) || isLine1FanVertex(a[1])) {
                verbose_assert(isLine1FanVertex(a[0]),
                               "edge orientation wrong");
            }
        }

        for (size_t iSeg = 0; iSeg < segments.size(); ++iSeg) {
            if (segmentTypes[iSeg] == LineType::A) {
                Segment a = segments[iSeg];
                if (isLine1FanVertex(a[0]) || isLine1FanVertex(a[1])) {
                    verbose_assert(isLine1FanVertex(a[0]),
                                   "segment orientation wrong");
                }
            }
        }
    }


    // Split boundary lines
    std::vector<Line> boundaryLines{Line{1, 2}, Line{2, 0}, Line{0, 1}};

    // Storing the index of the first boundary segments lets us
    // loop over interior segments later. Looping over interior
    // segments is nice because each interior segment is in
    // exactly two faces
    size_t firstBoundarySegment = segments.size();
    splitIntoSegments(boundaryLines, boundaryVertices, vertices,
                      LineType::Boundary, vertexSegments,
                      vertexSegmentDirections, segments, segmentTypes, true);

    // Do weird stuff for segments incident on fan vertex
    // Map storing the values on the left and right of segments incident
    // on the fan vertex for lines1. This is special because it could
    // have arbitrarily high degree, unlike all other vertices which are
    // incident on at most 2 segments from lines1 If lines1 has no fan
    // edge, this map is empty
    //
    // Note that because of how getLines and splitIntoSegments work,
    // vertexSegments[lines1FanEdge] consists of the segments incident
    // on the fan vertex, ordered by the barycentric coordinate of the
    // far end of the line (the end on the long edge) Luckily, this is
    // the order that longSideCornerVal1 gives its values in
    //
    // Unfortunately, all segments from lines1 come before all segments
    // from lines2. So if both sets have the same fan vertex, the
    // ordering gets confusing. In this case, we compute which value
    // from lines2 goes in the lines1 corners and vice versa ahead of
    // time using barycentric coordinates along the fan edge
    //
    // Furthermore, we use the fact that all these segments are oriented
    // pointing away from the fan vertex
    std::unordered_map<size_t, std::array<T1, 2>> fan1SegmentValues;
    if (verbose) {
        cout << "LINES 1 FAN EDGE: " << lines1FanEdge << endl;
        for (size_t iE = 0; iE < 3; ++iE)
            cout << "\t" << edgeVertices1[iE].size() << endl;
    }
    if (lines1FanEdge < 3) { // check if lines1 has a fan edge
        // TODO: this computation is copied in getFanEdge
        // Note: lastFanLine1 is inclusive
        size_t firstFanLine1 = edgeVertices1[prev(lines1FanEdge)].size();
        size_t lastFanLine1  = edgeVertices1[lines1FanEdge].size() -
                              edgeVertices1[next(lines1FanEdge)].size() - 1;

        std::vector<T1> longSideCornerVal1ForLines2;
        if (lines1FanEdge == lines2FanEdge) {
            size_t fE            = lines1FanEdge;
            size_t firstFanLine2 = edgeVertices2[prev(lines2FanEdge)].size();
            size_t lastFanLine2  = edgeVertices2[lines2FanEdge].size() -
                                  edgeVertices2[next(lines2FanEdge)].size() - 1;
            size_t iL1 = firstFanLine1;
            for (size_t iL2 = firstFanLine2; iL2 <= lastFanLine2; ++iL2) {
                while (bary1[fE][iL1 + 1] < bary2[fE][iL2 + 1] &&
                       iL1 - firstFanLine1 < longSideCornerVal1.size() - 1) {
                    iL1++;
                }
                longSideCornerVal1ForLines2.push_back(
                    longSideCornerVal1[iL1 - firstFanLine1]);
                if (verbose) {
                    cout << ">>Using iL1 = " << iL1 << " " << bary1[fE][iL1 + 1]
                         << ", " << bary2[fE][iL2 + 1] << endl;
                }
            }
        }

        if (longSideCornerVal1.size() == 0) {
            cout << "Error: fan edge, but no fan corner values" << endl;
            WATCH(lines1FanEdge);
            WATCH(edgeVertices1);
            WATCH(bary1);
            WATCH(val1);
        }

        verbose_assert(longSideCornerVal1.size() > 0, "huhb?");
        if (longSideCornerVal1[0] != val1[prev(lines1FanEdge)][0][1]) {
            cout << "firstFanLine: " << firstFanLine1 << endl;
            cout << "edge val: " << val1[prev(lines1FanEdge)][0][1] << endl;
            cout << "edge lens (inclusive): "
                 << val1[prev(lines1FanEdge)].size() << "\t"
                 << val1[lines1FanEdge].size() << "\t"
                 << val1[next(lines1FanEdge)].size() << endl;
            cout << "corner values: ";
            for (T1 val : longSideCornerVal1) {
                cout << "\t" << val;
            }
            cout << endl;
        }

        verbose_assert(longSideCornerVal1[0] == val1[prev(lines1FanEdge)][0][1],
                       "Boundary data incompatible with fan data");

        size_t iA = 0, iB = 0;
        // cout << "first fan line: " << firstFanLine << endl;
        T1 lastVal = longSideCornerVal1[iA];
        for (size_t iSeg : vertexSegments[lines1FanEdge]) {
            // Only pay attention to segments of type A incident on this
            // vertex
            if (segmentTypes[iSeg] == LineType::A) {
                if (verbose) {
                    cout << "A" << endl;
                }
                iA++;
                T1 nextVal              = longSideCornerVal1[iA];
                T1 leftVal              = nextVal;
                T1 rightVal             = lastVal;
                fan1SegmentValues[iSeg] = std::array<T1, 2>{leftVal, rightVal};

                lastVal = nextVal;
            } else if (segmentTypes[iSeg] == LineType::B) {
                T1 currentVal = longSideCornerVal1ForLines2[iB];
                if (verbose) {
                    cout << "OTHER" << endl;
                    cout << "iB: " << iB << "\tval: " << currentVal << endl;
                }
                iB++;
                fan1SegmentValues[iSeg] =
                    std::array<T1, 2>{currentVal, currentVal};
            }
        }

        size_t nFanSectors = lastFanLine1 + 1 - firstFanLine1;
        if (iA != nFanSectors) {
            cout << edgeVertices2[lines1FanEdge].size() << "\t"
                 << edgeVertices2[prev(lines1FanEdge)].size() << "\t"
                 << edgeVertices2[next(lines1FanEdge)].size() << endl;
            cout << edgeVertices1[lines1FanEdge].size() << "\t"
                 << edgeVertices1[prev(lines1FanEdge)].size() << "\t"
                 << edgeVertices1[next(lines1FanEdge)].size() << "\t|\t"
                 << vertexSegments[lines1FanEdge].size() << endl;
        }
        verbose_assert(iA == nFanSectors,
                       "wrong number of corners on val1 line. iCorner = " +
                           std::to_string(iA) +
                           " but nFanSectors = " + std::to_string(nFanSectors));

        if (longSideCornerVal1[longSideCornerVal1.size() - 1] !=
            val1[prev(lines1FanEdge)][0][0]) {
            cout << "firstFanLine: " << firstFanLine1 << endl;
            cout << "edge val: " << val1[prev(lines1FanEdge)][0][0] << "\t"
                 << val1[prev(lines1FanEdge)][0][1] << endl;
            cout << "edge lens: " << val1[prev(lines1FanEdge)].size() << "\t"
                 << val1[lines1FanEdge].size() << "\t"
                 << val1[next(lines1FanEdge)].size() << endl;
            cout << "corner values: ";
            for (T1 val : longSideCornerVal1) {
                cout << "\t" << val;
            }
            cout << endl;
        }

        verbose_assert(longSideCornerVal1[longSideCornerVal1.size() - 1] ==
                           val1[prev(lines1FanEdge)][0][0],
                       "Boundary data incompatible with fan data");
        // verbose_assert(longSideCornerVal1[lastFanLine + 1] ==
        //           val1[next(lines1FanEdge)]
        //           [val1[next(lines1FanEdge)].size() - 1][0],
        //           "Boundary data incompatible with fan data");
    }

    std::unordered_map<size_t, std::array<T2, 2>> fan2SegmentValues;
    if (lines2FanEdge < 3) { // check if lines1 has a fan edge
        // TODO: this computation is copied in getFanEdge
        // Note: lastFanLine1 is inclusive
        size_t firstFanLine2 = edgeVertices2[prev(lines2FanEdge)].size();
        size_t lastFanLine2  = edgeVertices2[lines2FanEdge].size() -
                              edgeVertices2[next(lines2FanEdge)].size() - 1;

        std::vector<T2> longSideCornerVal2ForLines1;
        if (lines2FanEdge == lines1FanEdge) {
            size_t fE            = lines2FanEdge;
            size_t firstFanLine1 = edgeVertices1[prev(lines1FanEdge)].size();
            size_t lastFanLine1  = edgeVertices1[lines1FanEdge].size() -
                                  edgeVertices1[next(lines1FanEdge)].size() - 1;
            size_t iL2 = firstFanLine2;
            for (size_t iL1 = firstFanLine1; iL1 <= lastFanLine1; ++iL1) {
                while (bary2[fE][iL2 + 1] < bary1[fE][iL1 + 1] &&
                       iL2 - firstFanLine2 < longSideCornerVal2.size() - 1) {
                    iL2++;
                }
                longSideCornerVal2ForLines1.push_back(
                    longSideCornerVal2[iL2 - firstFanLine2]);
                if (verbose) {
                    cout << "Using iL2 = " << iL2 << " " << bary2[fE][iL2 + 1]
                         << ", " << bary1[fE][iL1 + 1] << endl;
                }
            }
        }

        size_t iB = 0, iA = 0;
        // cout << "first fan line: " << firstFanLine << endl;
        T2 lastVal = longSideCornerVal2[iB];
        for (size_t iSeg : vertexSegments[lines2FanEdge]) {
            // Only pay attention to segments of type B incident on this
            // vertex
            if (segmentTypes[iSeg] == LineType::B) {
                iB++;
                T2 nextVal              = longSideCornerVal2[iB];
                T2 leftVal              = nextVal;
                T2 rightVal             = lastVal;
                fan2SegmentValues[iSeg] = std::array<T2, 2>{leftVal, rightVal};

                lastVal = nextVal;
            } else if (segmentTypes[iSeg] == LineType::A) {
                T2 currentVal = longSideCornerVal2ForLines1[iA];
                iA++;
                fan2SegmentValues[iSeg] =
                    std::array<T2, 2>{currentVal, currentVal};
                if (verbose) {
                    cout << iSeg << "\tinterp: " << currentVal << endl;
                }
            }
        }

        size_t lastFanLine = edgeVertices2[lines2FanEdge].size() -
                             edgeVertices2[next(lines2FanEdge)].size() - 1;

        size_t nFanSectors = lastFanLine2 + 1 - firstFanLine2;
        if (iB != nFanSectors) {
            cout << edgeVertices2[lines1FanEdge].size() << "\t"
                 << edgeVertices2[prev(lines1FanEdge)].size() << "\t"
                 << edgeVertices2[next(lines1FanEdge)].size() << endl;
            cout << edgeVertices1[lines1FanEdge].size() << "\t"
                 << edgeVertices1[prev(lines1FanEdge)].size() << "\t"
                 << edgeVertices1[next(lines1FanEdge)].size() << "\t|\t"
                 << vertexSegments[lines1FanEdge].size() << endl;
        }
        verbose_assert(iB == nFanSectors,
                       "wrong number of corners on val1 line. iCorner = " +
                           std::to_string(iB) +
                           " but nFanSectors = " + std::to_string(nFanSectors));
    }

    // ===========================================================
    //                      Extract polygons
    // ===========================================================

    // Assuming v is one vertex in edge e, returns the other one
    auto opp = [](size_t v, Segment e) { return e[0] + e[1] - v; };

    // Sort vertex segments to go clockwise around each vertex. This
    // makes the polygon extraction easier

    // Handle the triangle's corners separately (the corners are the
    // first 3 vertices)

    /*
    // Split boundary lines
    std::vector<Line> boundaryLines{Line{1, 2}, Line{2, 0}, Line{0, 1}};

    // Storing the index of the first boundary segments lets us
    // loop over interior segments later. Looping over interior
    // segments is nice because each interior segment is in
    // exactly two faces
    size_t firstBoundarySegment = segments.size();
    splitIntoSegments(boundaryLines, boundaryVertices, vertices,
                      LineType::Boundary, vertexSegments, segments,
                      segmentTypes);
    */
    /*
    for (size_t iV = 0; iV < 3; ++iV) {
          // Boundary segments come after interior segments
          std::vector<size_t> interiorSegments = vertexSegments[iV];
          // Remove last two elements (boundary segments);
          size_t iB1 = interiorSegments.back();
          interiorSegments.pop_back();
          size_t iB2 = interiorSegments.back();
          interiorSegments.pop_back();

          Vector2 v     = vertices[iV];
          auto cmpAngle = [&](const size_t& iE, const size_t& iF) {
              Vector2 vE = vertices[opp(iV, segments[iE])] - v;
              Vector2 vF = vertices[opp(iV, segments[iF])] - v;

              double thetaE = atan2(vE.y, vE.x);
              double thetaF = atan2(vF.y, vF.x);

              return thetaE > thetaF;
          };
          std::sort(std::begin(interiorSegments),
    std::end(interiorSegments), cmpAngle);

          if (!cmpAngle(iB1, iB2)) std::swap(iB1, iB2);

          interiorSegments.insert(interiorSegments.begin(), iB1);
          interiorSegments.push_back(iB2);

          vertexSegments[iV] = interiorSegments;
    }

    // For general vertices, just sort incident edges by angle
    // TODO: we could use the cyclic ordering of line endpoints around
    the
    // triangle boundary. That sounds like a lot of work
    for (size_t iV = 3; iV < vertices.size(); ++iV) {
        Vector2 v     = vertices[iV];
        auto cmpAngle = [&](const size_t& iE, const size_t& iF) {
            Vector2 vE = vertices[opp(iV, segments[iE])] - v;
            Vector2 vF = vertices[opp(iV, segments[iF])] - v;

            double thetaE = atan2(vE.y, vE.x);
            double thetaF = atan2(vF.y, vF.x);

            return thetaE > thetaF;
        };
        std::sort(std::begin(vertexSegments[iV]),
    std::end(vertexSegments[iV]), cmpAngle);
    }
    */
    for (size_t iV = 0; iV < vertices.size(); ++iV) {
        std::vector<std::pair<size_t, size_t>> segmentsAndDirections;
        for (size_t iS = 0; iS < vertexSegments[iV].size(); ++iS) {
            // TODO: negating the index here underflows. I just want to
            // reverse the order, so I guess that's fine, but maybe I
            // should fix it
            segmentsAndDirections.push_back(
                {-cyclicIndex(vertexSegmentDirections[iV][iS]),
                 vertexSegments[iV][iS]});
        }

        std::sort(std::begin(segmentsAndDirections),
                  std::end(segmentsAndDirections));

        for (size_t iS = 0; iS < vertexSegments[iV].size(); ++iS) {
            vertexSegments[iV][iS] = segmentsAndDirections[iS].second;
        }
    }

    auto indexOf = [](size_t elem, std::vector<size_t> vec) {
        for (size_t i = 0; i < vec.size(); ++i) {
            if (vec[i] == elem) return i;
        }
        cout << "failed to find elem " << elem << " in vec {";
        for (size_t i : vec) {
            cout << i << ", ";
        }
        cout << " }" << endl;
        exit(1);
    };


    // returns which side of iSeg is the left side of the oriented edge
    // pointing in towards the vertex. If iSeg points towards the
    // vertex, then it is the left (i.e. 0) of iSeg. Otherwise it is the
    // right
    auto incomingLeft = [&](size_t vertex, size_t iSeg) {
        return (segments[iSeg][1] == vertex) ? 0 : 1;
    };

    // returns which side of iSeg is the left side of the oriented edge
    // pointing in towards the vertex. If iSeg points towards the
    // vertex, then it is the left (i.e. 0) of iSeg. Otherwise it is the
    // right
    auto outgoingLeft = [&](size_t vertex, size_t iSeg) {
        return (segments[iSeg][0] == vertex) ? 0 : 1;
    };

    auto getVal1 = [&](size_t vertex, double iSegIn, double iSegOut) -> T1 {
        T1 result;
        if (isLine1FanVertex(vertex)) {
            // Since we move counterclockwise, we always care about the
            // corner to the left of the segment. Grab Val1 from
            // whichever segment has type A
            if (segmentTypes[iSegIn] == LineType::A) {
                size_t leftSide = incomingLeft(vertex, iSegIn);

                if (fan1SegmentValues.find(iSegIn) == fan1SegmentValues.end()) {
                    cerr << "Error: key not present in map (1)" << endl;
                }

                result = fan1SegmentValues[iSegIn][leftSide];
                // return fan1SegmentValues[iSegIn][leftSide];

            } else if (segmentTypes[iSegOut] == LineType::A) {
                size_t leftSide = outgoingLeft(vertex, iSegOut);

                if (fan1SegmentValues.find(iSegOut) ==
                    fan1SegmentValues.end()) {
                    cerr << "Error: key not present in map (2)" << endl;
                }


                result = fan1SegmentValues[iSegOut][leftSide];
                // return fan1SegmentValues[iSegOut][leftSide];
            } else {
                // If neither segment has type A, then this is the
                // intersection of an B line with another B line or the
                // boundary. In this case, both possible values of val1
                // are the same on the B line So we just take the
                // segment of type B. Note that both segments cannot be
                // boundary segments because this is a fan vertex

                // TODO: fix this for real. It shouldn't be possible to
                // have 2 boundary segments meeting here
                if (segmentTypes[iSegIn] == LineType::Boundary &&
                    segmentTypes[iSegOut] == LineType::Boundary) {
                    return vertexValues1[vertex][0];
                }

                size_t iSeg =
                    (segmentTypes[iSegIn] == LineType::B) ? iSegIn : iSegOut;

                if (fan1SegmentValues.find(iSeg) == fan1SegmentValues.end()) {
                    cerr << "Error: key not present in map (3)" << endl;
                }


                // if (verbose && vertex == 0) {
                //     cout << "Segment type: " << segmentTypes[iSeg] <<
                //     endl; cout << "\t options were: " <<
                //     segmentTypes[iSegIn] << ",
                //     "
                //          << segmentTypes[iSegOut] << endl;
                //     cout << "\t\t\t" << fan1SegmentValues[iSeg][0] <<
                //     "\t"
                //          << fan1SegmentValues[iSeg][1] << endl;
                // }
                result = fan1SegmentValues[iSeg][0];
                // return fan1SegmentValues[iSeg][0];
            }
        } else {
            // Since we move counterclockwise, we always care about the
            // corner to the left of the segment. Grab Val1 from
            // whichever segment has type A
            if (segmentTypes[iSegIn] == LineType::A) {
                size_t leftSide = incomingLeft(vertex, iSegIn);
                result          = vertexValues1[vertex][leftSide];
                // return vertexValues1[vertex][leftSide];

            } else if (segmentTypes[iSegOut] == LineType::A) {
                size_t leftSide = outgoingLeft(vertex, iSegOut);
                result          = vertexValues1[vertex][leftSide];
                // return vertexValues1[vertex][leftSide];
            } else {
                // If neither segment has type A, then this is the
                // intersection of a B line with the boundary. In this
                // case, both possible values of val1 are the same
                result = vertexValues1[vertex][0];
                // return vertexValues1[vertex][0];
            }
        }

        // if (norm(getX(result)) < 1e-8 && false) {
        //     cerr << "val is 0 [getVal1]" << endl;
        //     cerr << "\t isLine1FanVertex(vertex) = "
        //          << (isLine1FanVertex(vertex) ? "true" : "false") <<
        //          endl;
        //     cerr << "\tsegmentTypes[iSegIn] = " << segmentTypes[iSegIn]
        //     << endl; cerr << "\tsegmentTypes[iSegOut] = " <<
        //     segmentTypes[iSegOut]
        //          << endl;
        //     cerr << "\tiSegIn: " << iSegIn << endl;
        //     cerr << "\tiSegOut: " << iSegOut << endl;
        //     cerr << "\tfirstLine1FanCrossing: " << firstLine1FanCrossing
        //          << endl;
        //     cerr << "\tlastLine1FanCrossing: ???" << endl;
        //     cerr << "\tvertexSegments.size(): " <<
        //     vertexSegments[vertex].size()
        //          << endl;
        // }

        return result;
    };

    auto getVal2 = [&](size_t vertex, double iSegIn, double iSegOut) {
        if (isLine2FanVertex(vertex)) {
            // Since we move counterclockwise, we always care about the
            // corner to the left of the segment. Grab Val2 from
            // whichever segment has type B
            if (segmentTypes[iSegIn] == LineType::B) {
                size_t leftSide = incomingLeft(vertex, iSegIn);
                if (vertex == 1 && verbose) {
                    cout << "HERE" << endl;
                }
                return fan2SegmentValues[iSegIn][leftSide];

            } else if (segmentTypes[iSegOut] == LineType::B) {
                size_t leftSide = outgoingLeft(vertex, iSegOut);
                if (vertex == 1 && verbose) {
                    cout << "OR HERE" << endl;
                }
                return fan2SegmentValues[iSegOut][leftSide];
            } else {
                // If neither segment has type B, then this is the
                // intersection of an A line with another A line or the
                // boundary. In this case, both possible values of val1
                // are the same on the A line So we just take the
                // segment of type A. Note that both segments cannot be
                // boundary segments because this is a fan vertex

                size_t iSeg =
                    (segmentTypes[iSegIn] == LineType::A) ? iSegIn : iSegOut;
                if (vertex == 1 && verbose) {
                    cout << "Really HERE? " << iSeg << endl;
                    cout << (segmentTypes[iSeg] == LineType::A ? "A"
                                                               : "BOUNDARY")
                         << endl;
                    cout << "\t" << fan2SegmentValues[iSeg][0] << "\t"
                         << fan2SegmentValues[iSeg][1] << endl;
                }
                return fan2SegmentValues[iSeg][0];
            }
        } else {
            // Since we move counterclockwise, we always care about the
            // corner to the left of the segment. Grab Val2 from
            // whichever segment has type B
            if (segmentTypes[iSegIn] == LineType::B) {
                size_t leftSide = incomingLeft(vertex, iSegIn);
                return vertexValues2[vertex][leftSide];

            } else if (segmentTypes[iSegOut] == LineType::B) {
                size_t leftSide = outgoingLeft(vertex, iSegOut);
                return vertexValues2[vertex][leftSide];
            } else {
                // If neither segment has type A, then this is the
                // intersection of a B line with the boundary. In this
                // case, both possible values of val1 are the same
                return vertexValues2[vertex][0];
            }
        }
    };

    if (verbose) {
        cout << "input Data:" << endl;
        cout << "difference: " << bary1[0][1] - bary2[0][1] << endl;
        for (size_t iE = 0; iE < 3; ++iE) {
            for (auto ptData : val1[iE]) {
                cout << "\t" << ptData[0] << ", " << ptData[1];
            }
            cout << endl;
            cout << "\t";
            for (double t : bary1[iE]) {
                cout << "\t" << t;
            }
            cout << endl;
        }
        for (auto ptData : longSideCornerVal1) {
            cout << "\t" << ptData;
        }
        cout << endl;
        cout << "input Data:" << endl;
        for (size_t iE = 0; iE < 3; ++iE) {
            for (auto ptData : val2[iE]) {
                cout << "\t" << ptData[0] << ", " << ptData[1];
            }
            cout << endl;
            cout << "\t";
            for (double t : bary2[iE]) {
                cout << "\t" << t;
            }
            cout << endl;
        }
        for (auto ptData : longSideCornerVal2) {
            cout << "\t" << ptData;
        }
        cout << endl;
    }

    std::vector<char> inFaceForwards(segments.size());
    std::vector<char> inFaceBackwards(segments.size());
    auto extractClockwisePolygon = [&](size_t iSeg, bool forwards = true) {
        std::vector<size_t> poly;
        std::vector<T1> polyVal1;
        std::vector<T2> polyVal2;

        // Don't construct faces for edges which are already in
        // some face
        if ((forwards && inFaceForwards[iSeg]) ||
            (!forwards && inFaceBackwards[iSeg])) {
            return std::make_tuple(false, poly, polyVal1, polyVal2);
        }

        Segment seg   = segments[iSeg];
        size_t vStart = (forwards) ? seg[0] : seg[1];
        size_t vCurr  = vStart;
        size_t vNext  = opp(vCurr, seg);

        bool first = true;
        double faceParent;

        size_t iter = 0;

        if (verbose) {
            cout << "FACE: vStart = " << vStart << ", vNext = " << vNext
                 << ", forwards = " << forwards << endl;
        }
        do {
            size_t localSegIdx = indexOf(iSeg, vertexSegments[vNext]);

            size_t N           = vertexSegments[vNext].size();
            size_t newSegIndex = vertexSegments[vNext][(localSegIdx + 1) % N];

            poly.push_back(vNext);
            polyVal1.push_back(getVal1(vNext, iSeg, newSegIndex));
            polyVal2.push_back(getVal2(vNext, iSeg, newSegIndex));

            if (verbose) {
                // if (vNext == lines1FanEdge) {
                //     cout << "\t\t(I'm at the fan) v" << endl;
                // }
                // cout << "\t\t\t isLine1FanVertex(v): "
                //      << isLine1FanVertex(vNext) << endl;
                T1 c1 = getVal1(vNext, iSeg, newSegIndex);
                T2 c2 = getVal2(vNext, iSeg, newSegIndex);
                cout << "\t" << c1 << "\t" << c2 << "\t vNext:  " << vNext
                     << endl;
                // if (vNext == 1) {
                //     cout << "vNext(1) values: " <<
                //     vertexValues2[vNext][0]
                //          << ", " << vertexValues2[vNext][1] << endl;
                //     cout << "isFan: " << isLine2FanVertex(vNext) <<
                //     endl; cout << "getVal2: " << getVal2(vNext, iSeg,
                //     newSegIndex)
                //          << endl;
                // }
            }

            iSeg  = newSegIndex;
            seg   = segments[iSeg];
            vCurr = vNext;
            vNext = opp(vCurr, seg);

            if (seg[0] == vCurr) {
                inFaceForwards[iSeg] = true;
            } else {
                inFaceBackwards[iSeg] = true;
            }
        } while (vCurr != vStart);

        return std::make_tuple(true, poly, polyVal1, polyVal2);
    };

    std::vector<std::vector<size_t>> polygons;
    std::vector<std::vector<T1>> cornerValues1;
    std::vector<std::vector<T2>> cornerValues2;

    bool foundForwardsPolygon, foundBackwardsPolygon;
    std::vector<size_t> polygon;
    std::vector<T1> polygonCorners1;
    std::vector<T2> polygonCorners2;

    if (firstBoundarySegment > 0) {
        // If there's at least one interior segment, then we can
        // restrict ourselves to interior segments when looking for
        // faces
        for (size_t iSeg = 0; iSeg < firstBoundarySegment; ++iSeg) {
            std::tie(foundForwardsPolygon, polygon, polygonCorners1,
                     polygonCorners2) = extractClockwisePolygon(iSeg);
            if (foundForwardsPolygon) {
                polygons.push_back(polygon);
                cornerValues1.push_back(polygonCorners1);
                cornerValues2.push_back(polygonCorners2);
            }

            std::tie(foundBackwardsPolygon, polygon, polygonCorners1,
                     polygonCorners2) = extractClockwisePolygon(iSeg, false);
            if (foundBackwardsPolygon) {
                polygons.push_back(polygon);
                cornerValues1.push_back(polygonCorners1);
                cornerValues2.push_back(polygonCorners2);
            }
        }
    } else {
        // If there are no interior segments, then the triangle does not
        // get split
        // TODO: exit early in this special case

        for (size_t iE = 0; iE < 3; ++iE) {
            if (val1[iE].size() > 2 || val2[iE].size() > 2) {
                cout << "bary1: " << endl;
                for (size_t jE = 0; jE < 3; ++jE) {
                    cout << "\t";
                    for (double t : bary1[jE]) cout << " " << t;
                    cout << endl;
                }
                cout << "bary2: " << endl;
                for (size_t jE = 0; jE < 3; ++jE) {
                    cout << "\t";
                    for (double t : bary2[jE]) cout << " " << t;
                    cout << endl;
                }
            }
            verbose_assert(val1[iE].size() == 2,
                           "a triangle with edge crossings should have "
                           "interior segments");
            verbose_assert(val2[iE].size() == 2,
                           "a triangle with edge crossings should have "
                           "interior segments");
        }

        // Vertices are in a weird order because we order them by
        // opposite halfedge
        polygons.push_back(std::vector<size_t>{1, 2, 0});
        std::vector<T1> triCornerValues1;
        std::vector<T2> triCornerValues2;
        for (size_t i = 0; i < 3; ++i) {
            triCornerValues1.push_back(val1[i][0][0]);
            triCornerValues2.push_back(val2[i][0][0]);
        }
        cornerValues1.push_back(triCornerValues1);
        cornerValues2.push_back(triCornerValues2);
    }

    return std::make_tuple(vertices.size(), polygons, cornerValues1,
                           cornerValues2);
}
} // namespace ImplementationDetails
} // namespace CEPS
