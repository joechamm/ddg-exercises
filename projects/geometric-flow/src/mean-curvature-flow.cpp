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


    // Note: Geometry Central has linear solvers: https://geometry-central.net/numerical/linear_solvers/
    // Note: Update positions via geometry->inputVertexPositions
    SparseMatrix<double> M = geometry->massMatrix();
    SparseMatrix<double> H = buildFlowOperator(M, h);

    Eigen::MatrixXd currentPositions(mesh->nVertices(), 3);

    for(Vertex v : mesh->vertices()) {
        Vector3 pos = geometry->inputVertexPositions[v];
        currentPositions.row(v.getIndex()) = Eigen::Vector3d(pos.x, pos.y, pos.z);
    }

    geometrycentral::PositiveDefiniteSolver<double> solver(H);
    Eigen::MatrixXd newPositions(mesh->nVertices(), 3);

    for(size_t i = 0; i < 3; ++i) {
        Vector<double> b = M * currentPositions.col(i);
        newPositions.col(i) = solver.solve(b);
    }

    for(Vertex v : mesh->vertices()) {
        Eigen::Vector3d newPos = newPositions.row(v.getIndex());
        geometry->inputVertexPositions[v] = {newPos[0], newPos[1], newPos[2]};
    }
}