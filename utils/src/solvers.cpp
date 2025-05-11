#include "solvers.h"
#include "geometrycentral/numerical/linear_solvers.h"

/*
 * Compute the inverse of a sparse diagonal matrix.
 *
 * Input: A sparse diagonal matrix <M>.
 * Returns: The inverse of M, which is also a sparse diagonal matrix.
 */
SparseMatrix<double> sparseInverseDiagonal(SparseMatrix<double>& M) {

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    SparseMatrix<double> inv(M.rows(), M.cols());
    for (int i = 0; i < M.rows(); i++) {
        tripletList.push_back(T(i, i, 1.0 / M.coeffRef(i, i)));
    }
    inv.setFromTriplets(tripletList.begin(), tripletList.end());
    return inv;
}

/*
 * Computes the residual of Ax - λx, where x has unit norm and λ = x.Ax.
 *
 * Input: <A>, the complex sparse matrix whose eigendecomposition is being computed; and <x>, the current guess for the
 * smallest eigenvector
 * Returns: The residual
 */
double residual(const SparseMatrix<std::complex<double>>& A, const Vector<std::complex<double>>& x) {

    // TODO
    // Compute the residual of Ax - λx, where x has unit norm and λ = x.Ax.
    // Vector<std::complex<double>> y = x.normalized(); // get unit x
    // Vector<std::complex<double>> v = A * y - (y.dot(A * y)) * y;
    // return v.norm();
    std::complex<double> lambda = x.adjoint() * A * x;

    return (A * x - lambda * x).norm();
}

/*
 * Solves Ax = λx, where λ is the smallest nonzero eigenvalue of A, and x is the corresponding eigenvector.
 *
 * Input: <A>, the complex positive definite sparse matrix whose eigendecomposition is being computed.
 * Returns: The smallest eigenvector of A.
 */
Vector<std::complex<double>> solveInversePowerMethod(const SparseMatrix<std::complex<double>>& A) {
    // double eps = 1e-10;
    // Vector<std::complex<double>> x(A.cols());
    // x.setOnes().normalize(); // initial guess
    // double res = residual(A, x);
    // while (res > eps) {
    //
    //     // TODO
    //     // Compute the inverse of A and multiply it by x to get y.
    //     Vector<std::complex<double>> y = A * x;
    //     std::complex<double> mean = y.mean();
    // }

    //// TODO
    // return Vector<std::complex<double>>::Zero(1);
    Vector<std::complex<double>> y = Vector<std::complex<double>>::Random(A.rows()).normalized();
    SparseMatrix<std::complex<double>> cpA = A;
    geometrycentral::PositiveDefiniteSolver<std::complex<double>> solver(cpA);
    Vector<std::complex<double>> ynext = y;
    double eps = 1e-10;
    double res = residual(cpA, y);
    while (res > eps) {
        solver.solve(ynext, y);
        y = ynext - Vector<std::complex<double>>::Constant(ynext.size(), ynext.mean());
        y.normalize();
        res = residual(cpA, y);
        ynext = y;       
    }

    return y;
}