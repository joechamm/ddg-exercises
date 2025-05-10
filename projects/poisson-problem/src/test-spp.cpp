#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/surface_mesh_factories.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "scalar-poisson-problem.h"
#include "gtest/gtest.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

namespace {

class ScalarPoissonProblemTest : public ::testing::Test {

  protected:
    ScalarPoissonProblem SPP;
    Vector<double> rho;
    Vector<double> phi_soln;
    SparseMatrix<double> laplace_soln;
    std::vector<double> mass_soln;

    ScalarPoissonProblemTest() {

        std::vector<double> r;
        std::vector<double> phi;
        typedef Eigen::Triplet<double> T;
        std::vector<T> tripletList;

        std::vector<Vector3> v;
        std::vector<std::array<int, 3>> f;

        std::string filepath = "../include/test-spp-soln.txt";
        std::ifstream input_file(filepath);
        std::string line;
        std::string X;

        if (input_file.is_open()) {
            while (!input_file.eof()) {
                getline(input_file, line);
                std::istringstream iss(line);
                iss >> X;
                if (X == "v") {
                    double x, y, z;
                    iss >> x >> y >> z;
                    v.push_back({x, y, z});

                } else if (X == "f") {
                    int a, b, c;
                    iss >> a >> b >> c;
                    f.push_back({a, b, c});

                } else if (X == "rho") {
                    double val;
                    iss >> val;
                    r.push_back(val);
                } else if (X == "phi") {
                    double val;
                    iss >> val;
                    phi.push_back(val);
                } else if (X == "T") {
                    double val, row, col;
                    iss >> val >> row >> col;
                    tripletList.push_back(T(row, col, val));
                } else if (X == "mass") {
                    double val;
                    iss >> val;
                    mass_soln.push_back(val);
                }
            }
        } else {
            std::cerr << "Oops, could not open input file <" << filepath
                      << "> for unit testing. Make sure filepath is correct." << std::endl;
            std::runtime_error("");
        }

        size_t nVertices = v.size();
        size_t nFaces = f.size();
        Eigen::MatrixXd vMat(nVertices, 3);
        Eigen::MatrixXi fMat(nFaces, 3);
        for (size_t i = 0; i < nVertices; i++) {
            vMat(i, 0) = v[i][0];
            vMat(i, 1) = v[i][1];
            vMat(i, 2) = v[i][2];
        }
        for (size_t i = 0; i < nFaces; i++) {
            fMat(i, 0) = f[i][0] - 1; // Geometry Central takes 0-indexed vertices
            fMat(i, 1) = f[i][1] - 1;
            fMat(i, 2) = f[i][2] - 1;
        }
        std::unique_ptr<ManifoldSurfaceMesh> mesh;
        std::unique_ptr<VertexPositionGeometry> geometry;
        std::tie(mesh, geometry) = makeManifoldSurfaceMeshAndGeometry(vMat, fMat);
        SPP = ScalarPoissonProblem(mesh.release(), geometry.release());

        size_t N = SPP.mesh->nVertices();
        rho = Vector<double>::Zero(N);
        phi_soln = Vector<double>::Zero(N);
        laplace_soln = SparseMatrix<double>(N, N);
        laplace_soln.setFromTriplets(tripletList.begin(), tripletList.end());
        for (size_t i = 0; i < N; i++) {
            rho[i] = r[i];
            phi_soln[i] = phi[i];
        }
    }

    virtual ~ScalarPoissonProblemTest() {
        delete SPP.mesh;
        delete SPP.geometry;
    }
};

// DEBUGGING 
void printCompareMatrix(const SparseMatrix<double>& sol, const SparseMatrix<double>& actual) {

    size_t solRows = sol.rows();
    size_t solCols = sol.cols();
    size_t actualRows = actual.rows();
    size_t actualCols = actual.cols();
    if (solRows == actualRows && solCols == actualCols) {
        std::cout << "Matrix sizes match. Dimensions: " << solRows << " x " << solCols << std::endl;    
    } else {
        std::cout << "Matrix sizes do not match. Dimensions: " << solRows << " x " << solCols << " vs. " << actualRows
                  << " x " << actualCols << std::endl
                  << "Ending Comparison..." << std::endl;
        return;
    }

    SparseMatrix<double> diff = sol - actual;
    std::cout << "Frobenius Norm of difference: " << diff.norm() << std::endl;
    std::cout << "Number of non-zero entries in the difference matrix: " << diff.nonZeros() << std::endl;
    std::cout << "Sum of non-zero entries in the difference matrix: " << diff.sum() << std::endl;
//    std::cout << "Difference Matrix: " << std::endl << diff << std::endl;
}

void dumpMatrix(const SparseMatrix<double>& mat, const std::string& filepath) {
 
    std::ofstream output_file(filepath);
    if (output_file.is_open()) {
        for (size_t i = 0; i < mat.outerSize(); ++i) {
            for (SparseMatrix<double>::InnerIterator it(mat, i); it; ++it) {
                output_file << it.value() << " " << it.row() << " " << it.col() << std::endl;
            }
        }
    } else {
        std::cerr << "Oops, could not open output file <" << filepath
                  << "> for unit testing. Make sure filepath is correct." << std::endl;
        std::runtime_error("");
    }
    output_file.close();
    std::cout << "Dumped matrix to file <" << filepath << ">" << std::endl;
}

void dumpVector(const Vector<double>& rho, const std::string& filepath) {
    std::ofstream output_file(filepath);
    if (output_file.is_open()) {
        output_file << "Vector size: " << rho.size() << std::endl;
        output_file << "L2 norm: " << rho.norm() << std::endl;
        output_file << rho << std::endl;
    } else {
        std::cerr << "Oops, could not open output file <" << filepath
                  << "> for unit testing. Make sure filepath is correct." << std::endl;
        std::runtime_error("");
    }
    output_file.close();
    std::cout << "Dumped vector to file <" << filepath << ">" << std::endl;
}

// Some sanity checks for the Laplace matrix and mass matrix.
TEST_F(ScalarPoissonProblemTest, laplaceMatrixSymmetric) {
    bool success = true;
    double eps = 1e-5;
    SparseMatrix<double>& A = SPP.A;
    size_t V = A.rows();
    for (size_t i = 0; i < V; i++) {
        for (SparseMatrix<double>::InnerIterator it(A, i); it; ++it) {
            if (abs(it.value() - A.coeffRef(it.col(), it.row())) > eps) {
                success = false;
                break;
            }
        }
    }
    EXPECT_TRUE(success) << "Laplace matrix is not symmetric";
}

TEST_F(ScalarPoissonProblemTest, laplaceMatrixRowSums) {
    bool success = true;
    double eps = 1e-5;
    SparseMatrix<double>& A = SPP.A; // A should equal A^T
    size_t V = A.rows();
    for (size_t i = 0; i < V; i++) {
        double rowSum = 0.;
        for (SparseMatrix<double>::InnerIterator it(A, i); it; ++it) {
            rowSum += it.value();
        }
        if (abs(rowSum) > eps) success = false;
    }

    std::cout << "Dumping Laplace matrix..." << std::endl;
    dumpMatrix(A, "laplace_matrix.txt");
    std::cout << "Dumping Laplace solution..." << std::endl;
    dumpMatrix(laplace_soln, "laplace_soln.txt");
    // DEBUGGING
    //if (!success) {
    //    std::cout << "Row sum of actual Laplace Matrix do not sum to 0. Dumping difference..." << std::endl;
    //    dumpMatrix(A, "laplace_matrix.txt");
    //    dumpMatrix(laplace_soln, "laplace_soln.txt");

    //}

    EXPECT_TRUE(success) << "Each row of Laplace matrix does not sum to zero";
}
TEST_F(ScalarPoissonProblemTest, laplaceMatrixDiagonalSums) {
    bool success = true;
    double eps = 1e-5;
    SparseMatrix<double>& A = SPP.A; // A should equal A^T
    size_t V = A.rows();
    for (size_t i = 0; i < V; i++) {
        double diag = 0.;
        double sum = 0.;
        for (SparseMatrix<double>::InnerIterator it(A, i); it; ++it) {
            if (it.row() != it.col()) {
                sum += it.value();
            } else {
                diag = it.value();
            }
        }
        if (abs(diag + sum) > eps) success = false;
    }
    EXPECT_TRUE(success)
        << "Off-diagonal entries of Laplace matrix do not sum to (negative) the diagonal entry in each row";
}

TEST_F(ScalarPoissonProblemTest, massMatrixDiagonal) {

    bool success = true;
    SparseMatrix<double>& M = SPP.M;
    size_t V = M.cols();
    for (size_t i = 0; i < V; i++) {
        for (SparseMatrix<double>::InnerIterator it(M, i); it; ++it) {
            if (it.row() != it.col()) {
                success = false;
                break;
            }
        }
    }
    EXPECT_TRUE(success) << "Mass matrix is not diagonal";
}

// Laplace/mass matrices also needed for geometric-flow
TEST_F(ScalarPoissonProblemTest, laplaceMatrix) {

    EXPECT_TRUE((laplace_soln - SPP.A).norm() < 1e-5) << "Laplace matrix is incorrect";
}

TEST_F(ScalarPoissonProblemTest, massMatrix) {

    bool success = true;
    for (Vertex v : SPP.mesh->vertices()) {
        size_t i = v.getIndex();
        if (abs(mass_soln[i] - SPP.M.coeffRef(i, i)) > 1e-5) {
            success = false;
            break;
        }
    }
    EXPECT_TRUE(success) << "Mass matrix is incorrect";
}

TEST_F(ScalarPoissonProblemTest, solve) {

    Vector<double> phi = SPP.solve(rho);
    // DEBUGGING
    std::cout << "Dumping phi..." << std::endl;
    dumpVector(phi, "phi.txt");
    std::cout << "Dumping phi_soln..." << std::endl;
    dumpVector(phi_soln, "phi_soln.txt");
    std::cout << "Dumping difference..." << std::endl;
    dumpVector(phi - phi_soln, "phi_diff.txt");
    
    EXPECT_TRUE((phi - phi_soln).norm() < 1e-6) << "Poisson problem solved incorrectly";
}
} // namespace

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}