// Implement member functions for SpectralConformalParameterization class.
#include "spectral-conformal-parameterization.h"

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
SpectralConformalParameterization::SpectralConformalParameterization(ManifoldSurfaceMesh* inputMesh,
                                                                     VertexPositionGeometry* inputGeo) {

    this->mesh = inputMesh;
    this->geometry = inputGeo;
}

/*
 * Builds the complex conformal energy matrix EC = ED - A.
 *
 * Input:
 * Returns: A complex sparse matrix representing the conformal energy
 */
SparseMatrix<std::complex<double>> SpectralConformalParameterization::buildConformalEnergy() const {

    // TODO
    // The conformal energy matrix is given by EC = ED - A, where ED is the discrete Laplace operator and A is the
    // diagonal mass matrix. The diagonal mass matrix is given by the barycentric dual area of each vertex.
    SparseMatrix<std::complex<double>> ED = geometry->complexLaplaceMatrix();
    size_t numBoundaryLoops = this->mesh->nBoundaryLoops();
    size_t numVertices = this->mesh->nVertices();

    std::vector<Eigen::Triplet<std::complex<double>>> triplets;
    //for (size_t bloop = 0; bloop < numBoundaryLoops; ++bloop) {
    //    BoundaryLoop bl = this->mesh->boundaryLoop(bloop);        
    //    for (Halfedge he : bl.adjacentHalfedges()) {
    //        Vertex zi = he.tailVertex();
    //        Vertex zj = he.tipVertex();
    //        size_t i = zi.getIndex();
    //        size_t j = zj.getIndex();
    //        triplets.push_back(Eigen::Triplet<std::complex<double>>(i, j, std::complex<double>(0, -0.25)));
    //        triplets.push_back(Eigen::Triplet<std::complex<double>>(j, i, std::complex<double>(0, 0.25)));
    //    }       
    //}

    for (Face f : mesh->faces()) {
        for (Halfedge he : f.adjacentHalfedges()) {
            size_t i = he.tailVertex().getIndex();
            size_t j = he.tipVertex().getIndex();
            triplets.push_back(Eigen::Triplet<std::complex<double>>(i, j, std::complex<double>(0, -0.25)));
            triplets.push_back(Eigen::Triplet<std::complex<double>>(j, i, std::complex<double>(0, 0.25)));
        }
    }

    SparseMatrix<std::complex<double>> A(numVertices, numVertices);
    A.setFromTriplets(triplets.begin(), triplets.end());
    return ED/2 - A;
}


/*
 * Flattens the input surface mesh with 1 or more boundaries conformally.
 *
 * Input:
 * Returns: A MeshData container mapping each vertex to a vector of planar coordinates.
 */
VertexData<Vector2> SpectralConformalParameterization::flatten() const {

    // TODO
//    return VertexData<Vector2>(*mesh); // placeholder
    SparseMatrix<std::complex<double>> EC = buildConformalEnergy();
    Vector<std::complex<double>> x = solveInversePowerMethod(EC);
    VertexData<Vector2> flattened = VertexData<Vector2>(*mesh);
    for (Vertex v : mesh->vertices()) {
        flattened[v] = Vector2::fromComplex(x[v.getIndex()]);
    }
    return flattened;
}