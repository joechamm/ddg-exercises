// Implement member functions for TrivialConnections class.
#include "trivial-connections.h"
#include "geometrycentral/numerical/linear_solvers.h"
#include "tree-cotree.h"
#include "harmonic-bases.h"
#include "hodge-decomposition.h"

typedef Eigen::Triplet<double> Triplet;
/*
 * Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
TrivialConnections::TrivialConnections(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    mesh = inputMesh;
    geometry = inputGeo;

    // Initialize the Hodge Decomposition
    HodgeDecomposition hodgeDecomp(mesh, geometry);

    // Initialize the Tree-Cotree structure and build the generators.
    TreeCotree treeCotree(mesh, geometry);
    treeCotree.buildGenerators();

    // Initialize the Harmonic Bases.
    HarmonicBases harmonicBases(mesh, geometry);

    // Build the Harmonic Bases using the generators and Hodge Decomposition.
    this->bases = harmonicBases.compute(treeCotree.generators, hodgeDecomp);

    // homology generators
    this->generators = treeCotree.generators;

    // Build the period matrix P.
    this->P = this->buildPeriodMatrix();
    // 0-form Laplacian
    this->A = hodgeDecomp.A;
    // 1-form Hodge star
    this->hodge1 = hodgeDecomp.hodge1;
    // 0-form exterior derivative
    this->d0 = hodgeDecomp.d0;

    // TODO: Build harmonic bases
    //    this->bases; // placeholder;
    //
    //    // Build period matrix.
    //    this->P = this->buildPeriodMatrix();
    //
    //    // TODO: Store DEC operators
    //    this->A = identityMatrix<double>(1);      // placeholder
    //    this->hodge1 = identityMatrix<double>(1); // placeholder
    //    this->d0 = identityMatrix<double>(1);     // placeholder
    //}
}

/*
 * Builds the period matrix Pij = ∑_{ek ∈ li} (ξj)k, where li is the ith homology generator, ek is a dual edge in li and
 * ξj is the jth harmonic 1-form basis.
 *
 * Input:
 * Returns: A sparse matrix represending the period matrix.
 */
SparseMatrix<double> TrivialConnections::buildPeriodMatrix() const {
    size_t nBases = this->bases.size();

    std::vector<Triplet> triplets;

    for (size_t i = 0; i < nBases; i++) {
        std::vector<Halfedge> generator = this->generators[i];

        for (size_t j = 0; j < nBases; j++) {
            Vector<double> basis = this->bases[j];
            double sum = 0.0;

            for (Halfedge he : generator) {
                size_t k = he.edge().getIndex();
                double orientation = (he == he.edge().halfedge()) ? 1.0 : -1.0;

                sum += orientation * basis[k];
            }

            triplets.push_back(Triplet(i, j, sum));
        }
    }

    SparseMatrix<double> periodMatrix(nBases, nBases);
    periodMatrix.setFromTriplets(triplets.begin(), triplets.end());
    return periodMatrix;
}

//SparseMatrix<double> TrivialConnections::buildPeriodMatrix() const {
//    // TODO
//    size_t nBases = this->bases.size();
//    SparseMatrix<double> P(nBases, nBases);
//
//    for (size_t i = 0; i < this->generators.size(); i++) {
//        for (size_t j = 0; j < this->bases.size(); j++) {
//            double sum = 0.0;
//            for (Halfedge he : this->generators[i]) {
//                double orientation = (he == he.edge().halfedge()) ? 1.0 : -1.0;
//                sum += orientation * this->bases[j][he.edge().getIndex()];
//            }
//
//            P.coeffRef(i, j) = sum;
//        }
//    
//    }
//
//    return P;
//    
//    //    return identityMatrix<double>(1); // placeholder
//}

/*
 * Determine if a mesh satisfies Gauss-Bonnet.
 *
 * Input: A vector where the ith entry is the the index of the singularity at the ith vertex.
 * Returns: True if mesh satisfies Gauss-Bonnet, false otherwise.
 */
bool TrivialConnections::satsifyGaussBonnet(const Vector<double>& singularity) const {

    return (abs(singularity.sum() - geometry->eulerCharacteristic()) < 1e-8);
}

/*
 * Compute the dual 0-form potential β by solving the system d𝛿β = -K + 2π * singularity.
 *
 * Input: A vector where the ith entry is the the index of the singularity at the ith vertex.
 * Returns: The coexact component 𝛿β.
 */
Vector<double> TrivialConnections::computeCoExactComponent(const Vector<double>& singularity) const {
    size_t nVerts = mesh->nVertices();
    Vector<double> u = Vector<double>::Zero(nVerts);

    for (Vertex v : mesh->vertices()) {
        size_t i = v.getIndex();
        double angleDefect = geometry->angleDefect(v);
        u[i] = -angleDefect + 2.0 * PI * singularity[i];
    }

    SparseMatrix<double> L = this->A;
    geometrycentral::PositiveDefiniteSolver<double> solver(L);
    Vector<double> betaTilde = solver.solve(u);
    
    return hodge1 * d0 * betaTilde; // Return the coexact component 𝛿β
}
//Vector<double> TrivialConnections::computeCoExactComponent(const Vector<double>& singularity) const {
//    Vector<double> u = Vector<double>::Zero(mesh->nVertices());
//    /*for (Vertex v : mesh->vertices()) {
//        size_t i = v.getIndex();
//        double angleDefect = geometry->angleDefect(v);
//        u[i] = 2.0 * M_PI * singularity[i] - angleDefect;
//    }
//
//    SparseMatrix<double> L = this->A;
//    Vector<double> beta = solvePositiveDefinite(L, u);
//    return hodge1 * d0 * beta;*/
//
//    // u = - K + 2π * singularity
//    for (Vertex v : mesh->vertices()) {
//        size_t i = v.getIndex();
//        double K_i = geometry->vertexGaussianCurvature(v);
//        u[i] = 2.0 * PI * singularity[i] - K_i;
//    }
//
//    SparseMatrix<double> d0T = this->d0.transpose();
//    geometrycentral::PositiveDefiniteSolver<double> solver(d0T);
//    Vector<double> deltaBeta = solver.solve(u);
//    return deltaBeta;
//    // TODO
// //   return Vector<double>::Zero(1); // placeholder
//}


/*
 * Given an initial angle αi in face i, this function computes the new angle αj in the neighboring face j as
 * αj = αi - θij + θji, where θij and θji are the angles between the shared edge e and an arbitrary but fixed reference
 * direction in faces i and j. Repeating this procedure for n consecutive dual edges in a generator gives a sequence of
 * angles α0, . . . , αn with a resulting total angle defect equal to αn - α0. This corresponds to transporting a vector
 * around a generator by unfolding, sliding and refolding it across neighboring faces without any extra in plane
 * rotation.
 *
 * Input: A halfedge lying on the shared edge between face i and j, and the initial angle αi.
 * Returns: The new angle αj.
 */
double TrivialConnections::transportNoRotation(Halfedge he, double alphaI) const {

    Vector3 u = geometry->halfedgeVector(he);

    // Compute two orthonormal tangent vectors for each face.
    Face fi = he.face();
    Face fj = he.twin().face();
    Vector3 e1 = geometry->halfedgeVector(fi.halfedge()).normalize();
    Vector3 e2 = cross(geometry->faceNormal(fi), e1);
    Vector3 f1 = geometry->halfedgeVector(fj.halfedge()).normalize();
    Vector3 f2 = cross(geometry->faceNormal(fj), f1);
    double thetaIJ = atan2(dot(u, e2), dot(u, e1));
    double thetaJI = atan2(dot(u, f2), dot(u, f1));

    return alphaI - thetaIJ + thetaJI;
}

/*
 * Compute the harmonic component γ = ∑_{i = 1, ..., 2g} zi ξi by solving the system Pz = v - ∑𝛿β.
 * v - ∑𝛿β should be normalized to lie between -π and π.
 *
 * Input: The coexact component 𝛿β.
 * Returns: The harmonic component γ.
 */
Vector<double> TrivialConnections::computeHarmonicComponent(const Vector<double>& deltaBeta) const {
    size_t nBases = this->bases.size();
    size_t nEdges = mesh->nEdges();

    Vector<double> gamma = Vector<double>::Zero(nEdges);

    if (nBases > 0) {
        SparseMatrix<double> P = this->P;
        geometrycentral::SquareSolver<double> solver(P);

        // construct right hand side
        Vector<double> rhs = Vector<double>::Zero(nBases);
        for (size_t i = 0; i < nBases; i++) {
            std::vector<Halfedge> generator = this->generators[i];
            double sum = 0.0;

            for (Halfedge he : generator) {
                size_t k = he.edge().getIndex();
                double orientation = (he == he.edge().halfedge()) ? 1.0 : -1.0;

                sum += transportNoRotation(he, 0.0);
                sum -= orientation * deltaBeta[k];
            }

            // Normalize the sum to lie between -π and π
            while (sum < -PI) {
                sum += 2.0 * PI;
            }
            while (sum >= PI) {
                sum -= 2.0 * PI;
            }

            rhs[i] = sum;
        }

        // Solve the system Pz = rhs
        Vector<double> z = solver.solve(rhs);

        // Compute gamma
        for (size_t i = 0; i < nBases; i++) {
            Vector<double> basis = this->bases[i];
            gamma += z[i] * basis;        
        }
    }

    return gamma;
}

//Vector<double> TrivialConnections::computeHarmonicComponent(const Vector<double>& deltaBeta) const {
//    Vector<double> gamma = Vector<double>::Zero(mesh->nEdges());
//
//    if (this->bases.size() > 0) {
//    
//        Vector<double> v(this->generators.size());
//        for (size_t i = 0; i < this->generators.size(); i++) {
//            double sum = 0.0;
//            for (Halfedge he : this->generators[i]) {
//                double orientation = (he == he.edge().halfedge() ? 1.0 : -1.0);
//                sum += transportNoRotation(he, 0.0);
//                sum -= orientation * deltaBeta[he.edge().getIndex()];                
//            }
//            v[i] = sum - 2.0 * M_PI * std::floor(sum / (2.0 * M_PI));
//        }
//
//        SparseMatrix<double> PeriodMatrix = this->P;
//
//        geometrycentral::SquareSolver<double> solver(PeriodMatrix);
//        Vector<double> z = solver.solve(v);
//
//        SparseMatrix<double> d1 = geometry->buildExteriorDerivative1Form();
//        SparseMatrix<double> hodge1 = geometry->buildHodgeStar1Form();
//        SparseMatrix<double> d0T = geometry->buildExteriorDerivative0Form().transpose();
//        for (size_t i = 0; i < this->bases.size(); i++) {
//            gamma += z[i] * bases[i];
//
//            Vector<double> dGamma = d1 * bases[i];
//            Vector<double> deltaGamma = d0T * hodge1 * bases[i];
//            if (dGamma.norm() > 1e-5) {
//                std::cout << "dGamma.norm() = " << dGamma.norm() << std::endl;
//            }
//
//            if (deltaGamma.norm() > 1e-5) {
//                std::cout << "deltaGamma.norm() = " << deltaGamma.norm() << std::endl;
//            }
//        }
//    }
//    // TODO
//    //return Vector<double>::Zero(1); // placeholder
//    return gamma;
//}

/*
 * Compute the dual 1-form connections φ = 𝛿β + γ.
 *
 * Input: A vector where the ith entry is the the index of the singularity at the ith vertex.
 * Returns: A vector representing the connections.
 */
Vector<double> TrivialConnections::computeConnections(const Vector<double>& singularity) const {

    if (!this->satsifyGaussBonnet(singularity)) {
        std::cerr << "Singularities do not add up to the Euler characteristic of the mesh" << std::endl;
        return Vector<double>::Zero(mesh->nEdges());
    }
    // TODO: Compute connections on topological spheres
    //return Vector<double>::Zero(1); // placeholder

    Vector<double> deltaBeta = this->computeCoExactComponent(singularity);
    Vector<double> gamma = this->computeHarmonicComponent(deltaBeta);
    return deltaBeta + gamma;
}