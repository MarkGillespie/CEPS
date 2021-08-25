namespace polyscope {

template <typename T>
SurfaceVertexProjectiveSphericalParameterizationQuantity*
addProjectiveSphericalParameterizationQuantity(SurfaceMesh& mesh,
                                               std::string name,
                                               const T& data) {

    SurfaceVertexProjectiveSphericalParameterizationQuantity* q =
        new SurfaceVertexProjectiveSphericalParameterizationQuantity(
            name,
            applyPermutation(standardizeVectorArray<glm::vec4, 4>(data),
                             mesh.vertexPerm),
            mesh);
    mesh.addQuantity(q);
    return q;
}
} // namespace polyscope
