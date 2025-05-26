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
    this->hodge1 = identityMatrix<double>(1); // placeholder
    this->hodge2 = identityMatrix<double>(1); // placeholder
    this->d0 = identityMatrix<double>(1);     // placeholder
    this->d1 = identityMatrix<double>(1);     // placeholder

    this->hodge1 = geometry->buildHodgeStar1Form();
    this->hodge2 = geometry->buildHodgeStar2Form();
    this->d0 = geometry->buildExteriorDerivative0Form();
    this->d1 = geometry->buildExteriorDerivative1Form();
    // this->hodge1Inv = geometry->buildInverseHodgeStar1Form();
    // this->hodge2Inv = geometry->buildInverseHodgeStar2Form();
    this->hodge1Inv = sparseInverseDiagonal(hodge1);
    this->hodge2Inv = sparseInverseDiagonal(hodge2);
    this->d0T = d0.transpose();
    this->d1T = d1.transpose();

    // TODO: Build operator inverses.
    // Hint: Use the sparseInverseDiagonal() in utils/src/solvers.cpp to invert sparse diagonal matrices.
    // this->hodge1Inv = identityMatrix<double>(1); // placeholder
    // this->hodge2Inv = identityMatrix<double>(1); // placeholder
    // this->d0T = identityMatrix<double>(1);       // placeholder
    // this->d1T = identityMatrix<double>(1);       // placeholder

    // TODO: Construct 0-form Laplace matrix.
    // Shift matrix by a small constant (1e-8) to make it positive definite.
    // this->A = identityMatrix<double>(1); // placeholder
    this->A = d0T * hodge1 * d0 + (1e-8 * identityMatrix<double>(mesh->nVertices()));

    // TODO: Construct 2-form matrix.
    // this->B = identityMatrix<double>(1); // placeholder
    this->B = d1 * hodge1Inv * d1T + (1e-8 * identityMatrix<double>(mesh->nFaces()));
}

/*
 * Compute the 0-form potential α by solving the system 𝛿dα = 𝛿ω.
 *
 * Input: A primal 1-form on the edges of the input mesh.
 * Returns: The exact component dα of ω.
 */
Vector<double> HodgeDecomposition::computeExactComponent(const Vector<double>& omega) const {

    // TODO
 //   return Vector<double>::Zero(1); // placeholder
    Vector<double> rhs = d0T * hodge1 * omega; // right-hand side of the system
    SparseMatrix<double> L = this->A; // Laplace matrix for the 0-form
    geometrycentral::PositiveDefiniteSolver<double> solver(L); // create a solver for the Laplace matrix
    Vector<double> alpha = solver.solve(rhs); // solve the system
    return d0 * alpha; // return the exact component dα of ω
}

/*
 * Compute the 2-form potential β by solving the system d𝛿β = dω.
 *
 * Input: A primal 1-form on the edges of the input mesh.
 * Returns: The coexact component 𝛿β of ω.
 */
Vector<double> HodgeDecomposition::computeCoExactComponent(const Vector<double>& omega) const {

    // TODO
   // return Vector<double>::Zero(1); // placeholder
    Vector<double> rhs = d1 * omega; // right-hand side of the system
    SparseMatrix<double> L = this->B; // Laplace matrix for the 2-form
    geometrycentral::SquareSolver<double> solver(L); // create a solver for the Laplace matrix
    Vector<double> beta = solver.solve(rhs); // solve the system
    return d1 * beta; // return the coexact component 𝛿β of ω
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



    // // TODO
    // return Vector<double>::Zero(1); // placeholder
    return omega - dAlpha - deltaBeta;
}