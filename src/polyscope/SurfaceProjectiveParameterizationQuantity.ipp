namespace polyscope {

template <typename T>
SurfaceCornerProjectiveParameterizationQuantity*
addProjectiveParameterizationQuantity(SurfaceMesh& mesh, std::string name,
                                      const T& data, ParamCoordsType type) {
    SurfaceCornerProjectiveParameterizationQuantity* q =
        new SurfaceCornerProjectiveParameterizationQuantity(
            name,
            applyPermutation(standardizeVectorArray<glm::vec3, 3>(data),
                             mesh.cornerPerm),
            type, ParamVizStyle::CHECKER, mesh);
    mesh.addQuantity(q);
    return q;
}
} // namespace polyscope
