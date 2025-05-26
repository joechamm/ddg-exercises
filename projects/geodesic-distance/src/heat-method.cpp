// Implement member functions HeatMethod class.
#include "heat-method.h"
#include "geometrycentral/numerical/linear_solvers.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
HeatMethod::HeatMethod(ManifoldSurfaceMesh* surfaceMesh, VertexPositionGeometry* geo) {

    this->mesh = surfaceMesh;
    this->geometry = geo;

    // TODO: Build Laplace and flow matrices.
    // Note: core/geometry.cpp has meanEdgeLength() function
    //this->A = identityMatrix<double>(1); // placeholder
    //this->F = identityMatrix<double>(1); // placeholder
    double meanEdgeLength = geometry->meanEdgeLength();
    double timeStep = meanEdgeLength * meanEdgeLength;
    size_t nVertices = mesh->nVertices();
    SparseMatrix<double> M = geometry->massMatrix();
    this->A = geometry->laplaceMatrix();
    this->F = M + timeStep * A; // flow matrix    
}

/*
 * Computes the vector field X = -∇u / |∇u|.
 *
 * Input: <u>, a dense vector representing the heat that is allowed to diffuse on the input mesh for a brief period of
 * time.
 * Returns: A MeshData container that stores a Vector3 per face.
 */
FaceData<Vector3> HeatMethod::computeVectorField(const Vector<double>& u) const {
    FaceData<Vector3> vecF(*mesh, {0, 0, 0}); // placeholder
    for (Face f : mesh->faces()) {
        Vector3 gradient = geometrycentral::Vector3::zero();    
        Vector3 faceNormal = geometry->faceNormal(f);

        for (Halfedge he : f.adjacentHalfedges()) {
            Vector3 edgeVector = geometry->inputVertexPositions[he.next().tipVertex()] -
                                 geometry->inputVertexPositions[he.next().tailVertex()];
            Vector3 edgePerp = edgeVector.rotateAround(faceNormal, M_PI / 2.0);
            gradient += edgePerp * (u[he.vertex().getIndex()]);
        }
        vecF[f.getIndex()] = -gradient.normalizeCutoff();
    }
    // TODO
    return vecF;
}

/*
 * Computes the integrated divergence ∇.X.
 *
 * Input: <X>, the vector field -∇u / |∇u| represented as a FaceData container
 * Returns: A dense vector
 */
Vector<double> HeatMethod::computeDivergence(const FaceData<Vector3>& X) const {
    Vector<double> divX = Vector<double>::Zero(mesh->nVertices());
    for (Face f : mesh->faces()) {
        Vector3 Xj = X[f.getIndex()];
        for (Halfedge he : f.adjacentHalfedges()) {
            Vector3 edgeVector =
                geometry->inputVertexPositions[he.tipVertex()] - geometry->inputVertexPositions[he.tailVertex()];
            double cot = geometry->cotan(he);
            double div = 0.5 * cot * dot(Xj, edgeVector);
            divX[he.tailVertex().getIndex()] += div;
            divX[he.tipVertex().getIndex()] -= div;
        }
    }
    // TODO
    return divX;
}

/*
 * Computes the geodesic distances φ using the heat method.
 *
 * Input: <delta>, a dense vector representing the heat sources, i.e., u0 = δ(x). Returns: A dense vector containing the
 * geodesic distances per vertex.
 */
Vector<double> HeatMethod::compute(const Vector<double>& delta) const {
    Eigen::SimplicialLLT<SparseMatrix<double>> llt(F);
    Vector<double> u = llt.solve(delta);
    FaceData<Vector3> X = computeVectorField(u);
    Vector<double> deltaPhi = computeDivergence(X);
    SparseMatrix<double> L = this->A;
    geometrycentral::PositiveDefiniteSolver<double> solver(L);
    
    Vector<double> phi = solver.solve(-deltaPhi);

    // Since φ is unique up to an additive constant, it should be shifted such that the smallest distance is zero
    this->subtractMinimumDistance(phi);

    return phi;
}