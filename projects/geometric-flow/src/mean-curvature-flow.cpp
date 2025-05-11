// Implement member functions for MeanCurvatureFlow class.
#include "mean-curvature-flow.h"
#include "geometrycentral/numerical/linear_solvers.h"

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
MeanCurvatureFlow::MeanCurvatureFlow(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    // Build member variables: mesh, geometry
    mesh = inputMesh;
    geometry = inputGeo;
}

/*
 * Build the mean curvature flow operator.
 *
 * Input: The mass matrix <M> of the mesh, and the timestep <h>.
 * Returns: A sparse matrix representing the mean curvature flow operator.
 */
SparseMatrix<double> MeanCurvatureFlow::buildFlowOperator(const SparseMatrix<double>& M, double h) const {

    // TODO
    SparseMatrix<double> L = geometry->laplaceMatrix();
    SparseMatrix<double> H = M + h * L;
    return H;
   
}

/*
 * Performs mean curvature flow.
 *
 * Input: The timestep <h>.
 * Returns:
 */
void MeanCurvatureFlow::integrate(double h) {

    // TODO
    // Note: Geometry Central has linear solvers: https://geometry-central.net/numerical/linear_solvers/
    // Note: Update positions via geometry->inputVertexPositions
    SparseMatrix<double> M = geometry->massMatrix();
    SparseMatrix<double> H = buildFlowOperator(M, h);

    Eigen::MatrixXd currentVpos(mesh->nVertices(), 3);

    for (Vertex v : mesh->vertices()) {
        Vector3 vPos = geometry->inputVertexPositions[v];
        currentVpos.row(v.getIndex()) = Eigen::Vector3d(vPos[0], vPos[1], vPos[2]);
    }

    geometrycentral::PositiveDefiniteSolver<double> solver(H);
    Eigen::MatrixXd nextVpos(mesh->nVertices(), 3);

    for (Eigen::Index c = 0; c < 3; ++c) {
        Vector<double> rhs = M * currentVpos.col(c);
        nextVpos.col(c) = solver.solve(rhs);
    }

    for (Vertex v : mesh->vertices()) {
        Eigen::Vector3d newPos = nextVpos.row(v.getIndex());
        geometry->inputVertexPositions[v] = {newPos.x(), newPos.y(), newPos.z()};
    }
}