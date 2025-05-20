// Implement member functions HeatMethod class.
#include "heat-method.h"

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

    // TODO
    return Vector<double>::Zero(1); // placeholder
}

/*
 * Computes the geodesic distances φ using the heat method.
 *
 * Input: <delta>, a dense vector representing the heat sources, i.e., u0 = δ(x). Returns: A dense vector containing the
 * geodesic distances per vertex.
 */
Vector<double> HeatMethod::compute(const Vector<double>& delta) const {

    // TODO
    Vector<double> phi = Vector<double>::Zero(delta.rows());

    // Since φ is unique up to an additive constant, it should be shifted such that the smallest distance is zero
    this->subtractMinimumDistance(phi);

    return phi;
}