// Implement member functions for ScalarPoissonProblem class.
#include "scalar-poisson-problem.h"
#include "geometrycentral/numerical/linear_solvers.h"

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
ScalarPoissonProblem::ScalarPoissonProblem(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    mesh = inputMesh;
    geometry = inputGeo;

    // TODO: Build member variables A (Laplace matrix), M (mass matrix), total area
    //this->A = identityMatrix<double>(1); // placeholder
    //this->M = identityMatrix<double>(1); // placeholder
    //this->totalArea = 0;                 // placeholder
    this->A = geometry->laplaceMatrix(); // Laplace matrix
    this->M = geometry->massMatrix();    // Mass matrix
    this->totalArea = geometry->totalArea(); // Total area of the mesh
}

/*
 * Computes the solution of the poisson problem Ax = -M(rho - rhoBar), where A is the POSITIVE DEFINITE Laplace matrix
 * and M is the mass matrix.
 *
 * Input: <rho>, the density of vertices in the mesh.
 * Returns: The solution vector.
 */
Vector<double> ScalarPoissonProblem::solve(const Vector<double>& rho) const {

    // TODO
    // Note: Geometry Central has linear solvers: https://geometry-central.net/numerical/linear_solvers/
    // First we need to find rhoBar, which is the integral of rho over the mesh divided by the total area.
    SparseMatrix<double> L = this->A; // Laplace matrix

    // double rhoBarVal = ((M * rho).sum()) / totalArea; // rhoBar = M * rho / totalArea
    // Vector<double> rhoBar =
    //     Vector<double>::Constant(rho.rows(), rhoBarVal); // rhoBar is a vector of the same size as rho

    // Vector<double> rhs = -M * (rho - rhoBar); // right-hand side of the equation
    // PositiveDefiniteSolver<double> solver(L); // create a solver for the Laplace matrix
    // return solver.solve(rho);                 // solve the equation Ax = -M(rho - rhoBar)
    double rhoBar = 0.0;
    for (Eigen::Index i = 0; i < rho.rows(); i++) {
        if (rho[i]) {
            rhoBar += rho[i] * this->M.coeff(i, i);
        }
    }
    rhoBar /= this->totalArea;

    geometrycentral::PositiveDefiniteSolver<double> solver(L);
    Vector<double> rhs = -this->M * (rho - rhoBar * Vector<double>::Ones(L.rows()));
    Vector<double> sol = solver.solve(rhs);
    return sol;
}