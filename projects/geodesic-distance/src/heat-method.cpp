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
    // Build Laplace Matrix, A
    this->A = geometry->laplaceMatrix();

    // Build Mass Matrix, M
    SparseMatrix<double> M = geometry->massMatrix();

    // Timestep should be t = m * h^2, where m is a constant (according to video, m = 1 is good choice) and 
    // h is average spacing between nodes (average edge length).

    double h = geometry->meanEdgeLength();
    double t = h * h; // m = 1

    // Build Flow Matrix, F
    // Trying F = M + tA
    this->F = M + t * A;
}

/*
 * Computes the vector field X = -∇u / |∇u|.
 *
 * Input: <u>, a dense vector representing the heat that is allowed to diffuse on the input mesh for a brief period of
 * time.
 * Returns: A MeshData container that stores a Vector3 per face.
 */
FaceData<Vector3> HeatMethod::computeVectorField(const Vector<double>& u) const {
    ManifoldSurfaceMesh& meshRef = *mesh;
    FaceData<Vector3> X(meshRef, {0,0,0});
    
    // Per Heat Method Paper, section 3.2, the gradient of u for triangle face f is given by
    //      ∇u = (1 / 2A_f) * Sum over i of u_i * (N x e_i)
    // where A_f is the area of triangle face f, u_i is the value of u at vertex i, N is the face normal,
    // and e_i is the (counter-clockwise oriented) edge opposite vertex i, and N x e_i is the cross product.
    for(Face f : meshRef.faces()) {
        Vector3 gradU = Vector3::zero();
        Vector3 N = geometry->faceNormal(f);
        double oneOver2A_f = (1.0 / (2.0 * geometry->faceArea(f)));
        for(Halfedge he : f.adjacentHalfedges()) {
            double u_i = u[he.vertex().getIndex()];
            Vector3 e_i = geometry->inputVertexPositions[he.next().tipVertex()] - geometry->inputVertexPositions[he.tipVertex()];
            gradU += (u_i * cross(N, e_i));
        }

        X[f] = - unit(gradU);
    }


    return X;
}

/*
 * Computes the integrated divergence ∇.X.
 *
 * Input: <X>, the vector field -∇u / |∇u| represented as a FaceData container
 * Returns: A dense vector
 */
Vector<double> HeatMethod::computeDivergence(const FaceData<Vector3>& X) const {
    // The integrated divergence associated with vertex i is given by
    //      ∇.X = (1/2) * Sum over f of [ cot(θ_1) * (dot(e_1, X_f) + cot(θ_2) * (dot(e_2, X_f)) ]
    // where θ_1 and θ_2 are the angles opposite to vertex i, e_1 is the edge vector opposite θ_1 and pointing 
    // away from vertex i, e_2 the edge vector opposite θ_2 and pointing away from vertex i, and X_f is the value
    // of the vector field X = -∇u / |∇u| at face f.

    Vector<double> divX = Vector<double>::Zero(mesh->nVertices());
    // for(Vertex v : mesh->vertices()) {
    //     for(Face f : v.adjacentFaces()) { // face f is oriented ijk
    //         Vector3 X_f = X[f];
    //         // e_1 is the edge vector from vertex i to vertex j, and e_2 is the edge vector from vertex i to vertex k.
    //         // We need to find the halfedges he_1 and he_2 that are oriented from vertex i to vertex j and from vertex i to vertex k. 
    //         // Halfedge he = f.halfedge();
    //         // Halfedge he_1, he_2;  // he_1 is oriented CCW, he_2 is oriented CW, and both have vertex v as their tail.
    //         // try {
    //         //     if(he.vertex() == v) { // he is e_ij, so he_1 = he. he.next() is e_jk, and he.next().next() is e_ki, so he_2 = he_next().next().twin() is e_ik.
    //         //         he_1 = he;
    //         //         he_2 = he.next().next().twin();
    //         //     } else if(he.tipVertex() == v) {    // he is e_ki, so he_1 = he.next() is e_ij, and he_2 = he.twin() is e_ik.
    //         //         he_1 = he.next();
    //         //         he_2 = he.twin();
    //         //     } else if(he.next().tipVertex() == v) { // he is e_jk since he.next() is e_ki, so he_1 = he.next().next() is e_ij and he_2 = he.next().twin() is e_ik. 
    //         //         he_1 = he.next().next();
    //         //         he_2 = he.next().twin();                
    //         //     } else {
    //         //         throw std::runtime_error("Error: Halfedge does not have vertex as tail or tip.");
    //         //     }
    //         // } catch (const std::exception& e) {
    //         //     // Handle or log the error as needed
    //         //     std::cerr << "Exception: " << e.what() << std::endl;               
    //         // }

    //         // Vector3 e_1 = geometry->inputVertexPositions[he_1.tipVertex()] - geometry->inputVertexPositions[v];
    //         // Vector3 e_2 = geometry->inputVertexPositions[he_2.tipVertex()] - geometry->inputVertexPositions[v];
    //         // double cotan_1 = geometry->cotan(he_1);
    //         // double cotan_2 = geometry->cotan(he_2);
    //         // divX[v.getIndex()] += (0.5 * (cotan_1 * dot(e_1, X_f) + cotan_2 * dot(e_2, X_f)));            

    //         Halfedge he_1 = v.halfedge();
    //         Halfedge he_2 = he_1.next().next().twin();
    //         Vector3 e_1 = geometry->halfedgeVector(he_1);
    //         Vector3 e_2 = geometry->halfedgeVector(he_2);
    //         double cotan_1 = geometry->cotan(he_1);
    //         double cotan_2 = geometry->cotan(he_2);
    //         divX[v.getIndex()] += (0.5 * (cotan_1 * dot(e_1, X_f) + cotan_2 * dot(e_2, X_f)));
    //     }
    // }

    for(Face f : mesh->faces()) {
        Vector3 X_f = X[f];
        for(Halfedge he : f.adjacentHalfedges()) {
            Vector3 e = geometry->inputVertexPositions[he.tipVertex()] - geometry->inputVertexPositions[he.tailVertex()];
            double val = 0.5 * dot(e, X_f) * (geometry->cotan(he));
            divX[he.tailVertex().getIndex()] += val;
            divX[he.tipVertex().getIndex()] -= val;
        }
    }

    return divX;

    // // Try using discrete exterior calculus to compute divergence
    // SparseMatrix<double> H1 = geometry->buildHodgeStar1Form();
    // SparseMatrix<double> D0 = geometry->buildExteriorDerivative0Form();
    // //SparseMatrix<double> H0 = geometry->buildHodgeStar0Form();
    // //const auto diagonalH0 = geometry->buildHodgeStar0Form().diagonal();
    // //const auto diagonalH0Matrix = diagonalH0.asDiagonal();
    // //const auto H0_inv = diagonalH0Matrix.inverse();
    // SparseMatrix<double> H0_inv = geometry->buildInverseHodgeStar0Form();
    // // SparseMatrix<double> H0_inv = H0.inverse();
    // // SparseMatrix<double> H0_inv = H0.transpose().inverse();
    // SparseMatrix<double> divOperator = H0_inv * D0.transpose() * H1;
    // Vector<double> divX = divOperator * X.toVector();
 //   return divX;

}

/*
 * Computes the geodesic distances φ using the heat method.
 *
 * Input: <delta>, a dense vector representing the heat sources, i.e., u0 = δ(x). Returns: A dense vector containing the
 * geodesic distances per vertex.
 */
Vector<double> HeatMethod::compute(const Vector<double>& delta) const {

    // We first need to find u, the diffused heat at time t, which is the solution to the equation Fu = u0, where F is the
    // flow matrix and u0 is the initial heat distribution. We can use the Eigen library to solve this linear system using the fact
    // that F is symmetric positive definite. 
    // SparseMatrix<double> M = this->F;
    // geometrycentral::PositiveDefiniteSolver<double> solver(M);
    // Vector<double> u = solver.solve(delta);

    // // Now that we have u, we can compute the vector field X = -∇u / |∇u|.
    // FaceData<Vector3> X = computeVectorField(u);
    // // Get the divergence of X now.
    // Vector<double> divX = computeDivergence(X);

    // // The geodesic distance φ is the solution to the Poisson equation Δφ = ∇.X, where Δ is the Laplace-Beltrami
    // // operator and ∇.X is the divergence of the vector field X. We can solve this linear system using the Eigen library.
    // // We need to build the Laplace matrix again, since we need to solve the equation Aφ = divX, where A is the Laplace
    // // matrix and φ is the geodesic distance.
    // SparseMatrix<double> L = this->A;
    // geometrycentral::PositiveDefiniteSolver<double> solver2(L);
    // Vector<double> phi = solver2.solve(divX);

    // // Since φ is unique up to an additive constant, it should be shifted such that the smallest distance is zero
    // this->subtractMinimumDistance(phi);

    // return phi;

     Eigen::SimplicialLLT<SparseMatrix<double>> llt(F);
     Vector<double> u = llt.solve(delta);
//    SparseMatrix<double> M = F;
//    Vector<double> u = geometrycentral::solve(M, delta);
    FaceData<Vector3> X = computeVectorField(u);
    Vector<double> divX = computeDivergence(X);

    SparseMatrix<double> L = A;
    geometrycentral::PositiveDefiniteSolver<double> solver(L);
    Vector<double> phi = solver.solve(-divX);
 //   Vector<double> phi = geometrycentral::solve(L, - divX);

    this->subtractMinimumDistance(phi);
    return phi;
}