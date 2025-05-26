// Implement member functions for HodgeDecomposition class.
#include "hodge-decomposition.h"
#include "geometrycentral/numerical/linear_solvers.h"

/*
 * Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
HodgeDecomposition::HodgeDecomposition(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    mesh = inputMesh;
    geometry = inputGeo;

    // TODO: build DEC operators
    this->hodge1 = geometry->buildHodgeStar1Form();
    this->hodge2 = geometry->buildHodgeStar2Form(); // placeholder
    this->d0 = geometry->buildExteriorDerivative0Form();     // placeholder
    this->d1 = geometry->buildExteriorDerivative1Form();     // placeholder

    // TODO: Build operator inverses.
    // Hint: Use the sparseInverseDiagonal() in utils/src/solvers.cpp to invert sparse diagonal matrices.
    this->hodge1Inv = sparseInverseDiagonal(this->hodge1);
    this->hodge2Inv = sparseInverseDiagonal(this->hodge2);
    this->d0T = this->d0.transpose();      // placeholder
    this->d1T = this->d1.transpose(); // placeholder

    // TODO: Construct 0-form Laplace matrix.
    // Shift matrix by a small constant (1e-8) to make it positive definite.
    this->A = d0T * hodge1 * d0 + identityMatrix<double>(mesh->nVertices()) * 1e-8; // placeholder

    // TODO: Construct 2-form matrix.
    this->B = d1 * hodge1Inv * d1T;
}

/*
 * Compute the 0-form potential α by solving the system 𝛿dα = 𝛿ω.
 *
 * Input: A primal 1-form on the edges of the input mesh.
 * Returns: The exact component dα of ω.
 */
Vector<double> HodgeDecomposition::computeExactComponent(const Vector<double>& omega) const {

    Vector<double> rhs = d0T * hodge1 * omega; // Compute the right-hand side of the equation 𝛿dα = 𝛿ω
    SparseMatrix<double> L = this->A;          // Use the Laplace matrix A
    geometrycentral::PositiveDefiniteSolver<double> solver(L); // Create a solver for the linear system
    Vector<double> alpha = solver.solve(rhs);                  // Solve the linear system to find the potential α
    return d0 * alpha;
}

/*
 * Compute the 2-form potential β by solving the system d𝛿β = dω.
 *
 * Input: A primal 1-form on the edges of the input mesh.
 * Returns: The coexact component 𝛿β of ω.
 */
Vector<double> HodgeDecomposition::computeCoExactComponent(const Vector<double>& omega) const {

    Vector<double> rhs = d1 * omega; // Compute the right-hand side of the equation d𝛿β = dω
    SparseMatrix<double> L = this->B; // Use the 2-form matrix B
    geometrycentral::SquareSolver<double> solver(L); // Create a solver for the linear system
    Vector<double> betaTilda = solver.solve(rhs);         // Solve the linear system to find the potential β
    Vector<double> deltaBeta = d1T * betaTilda;     // Compute the coexact component 𝛿β of ω
    return hodge1Inv * deltaBeta;                   // Return the coexact component 𝛿β of ω
}

/*
 * Compute the harmonic component γ = ω - dα - 𝛿β of ω.
 *
 * Input: A primal 1-form <omega> on the edges of the input mesh, the exact component <dAlpha> of ω, and the coexact
 * component <deltaBeta> of ω.
 * Returns: The coexact component 𝛿β of ω.
 */
Vector<double> HodgeDecomposition::computeHarmonicComponent(const Vector<double>& omega, const Vector<double>& dAlpha,
                                                            const Vector<double>& deltaBeta) const {

    return omega - dAlpha - deltaBeta; // placeholder
}