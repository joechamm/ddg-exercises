// Implement member functions for HarmonicBases class.
#include "harmonic-bases.h"
#include "geometrycentral/numerical/linear_solvers.h";

/*
 * Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
HarmonicBases::HarmonicBases(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    mesh = inputMesh;
    geometry = inputGeo;
}

Halfedge HarmonicBases::sharedHalfedge(Face f, Face g) const {
    for (Halfedge he : f.adjacentHalfedges()) {
        if (he.twin().face() == g) {
            return he;
        }
    }
    // If no shared halfedge is found, return an invalid Halfedge (or throw an exception if preferred).
    // This should not happen if the input faces are valid and share an edge.
    std::cerr << "Oops, HarmonicBases::sharedHalfedge: No shared halfedge found between faces " << f.getIndex()
              << " and " << g.getIndex() << std::endl;
    return f.halfedge(); // Return the first halfedge of face f as a fallback (not ideal, but prevents crash)
}

/*
 * Build a closed, but not exact, primal 1-form ω.
 *
 * Input: A std::vector of Halfedges representing a homology generator of the mesh.
 * Returns: A vector representing a closed primal 1-form.
 */
Vector<double> HarmonicBases::buildClosedPrimalOneForm(const std::vector<Halfedge>& generator) const {
    // Build teh 1-form omega from the generator. For every edge crossing from the "left" of the generator to the "right", set omega to +1, 
    // and for every edge crossing from the "right" to the "left", set omega to -1. All remaining edges should be set to 0. The resulting
    // 1-form is closed because the discrete exterior derivative on 1-forms is the (oriented) sum of edge values around each triangle. Since the generator is a cycle, this sum will
    // be zero for each triangle, thus making the 1-form closed.
    Vector<double> omega = Vector<double>::Zero(mesh->nEdges());

    // Iterate through the halfedges in the generator and assign values to the omega vector based on their orientation.
    // If the halfedge is oriented in the direction of the generator, assign a value of 1.0; otherwise, assign -1.0.
    // This ensures that the resulting 1-form is closed, as it will have a zero sum around each triangle in the mesh.
    // Note: The edge index is used to access the corresponding entry in the omega vector.
    // The orientation of the halfedge determines whether to assign a positive or negative value.

    Halfedge firstHe = generator.front();
    Halfedge lastHe = generator.back();
    double generatorOrientation =
        ((firstHe == firstHe.edge().halfedge())
             ? 1.0
             : -1.0); // Determine the orientation of the generator based on the first halfedge
       // Check if the generator is oriented in the same direction as the first halfedge
    //for (const Halfedge& he : generator) {
    //    size_t edgeIdx = he.edge().getIndex();
    //    
    //    Face face = he.face();
    //    face.
    //}

    size_t generatorSize = generator.size();

    for (size_t i = 0; i < generatorSize; ++i) {
        Halfedge he_i = generator[i];
        Halfedge he_j = generator[(i + 1) % generatorSize]; // Wrap around to the first halfedge after the last one
        Face face_i = he_i.face();
        Face face_j = he_j.face();
        Halfedge sharedHe = sharedHalfedge(face_i, face_j); // Get the shared halfedge between the two faces
        // Check if the shared halfedge is oriented in the same direction as the generator
        if (sharedHe.edge() == he_i.edge()) {
            // If the shared halfedge is oriented in the same direction as the generator, assign +1.0
            omega[he_i.edge().getIndex()] = generatorOrientation;
        } else {
            // If the shared halfedge is oriented in the opposite direction, assign -1.0
            omega[he_i.edge().getIndex()] = -generatorOrientation;
        }
    }



    //for (Halfedge he : generator) {
    //    omega[he.edge().getIndex()] = he.orientation() ? 1.0 : -1.0; // Assign a value of 1.0 to the halfedges in the generator
    //}
    // TODO
    return omega;
}

/*
 * Compute the harmonic bases [γ1, γ2 ... γn] of the input mesh.
 *
 * Input: A std::vector of homology generators of the mesh (which are in turn represented as std::vectors of halfedges),
 * and a HodgeDecomposition object. Returns:
 */
std::vector<Vector<double>> HarmonicBases::compute(const std::vector<std::vector<Halfedge>>& generators,
                                                   const HodgeDecomposition& hodgeDecomposition) const {

    SparseMatrix<double> d0 = hodgeDecomposition.d0;
    SparseMatrix<double> d0T = hodgeDecomposition.d0T;
    SparseMatrix<double> hodge1 = hodgeDecomposition.hodge1;
    SparseMatrix<double> codiff1 = d0T * hodge1;
    SparseMatrix<double> l0 = hodgeDecomposition.A;
    geometrycentral::PositiveDefiniteSolver<double> solver(l0);

    std::vector<Vector<double>> gammas;
    for (const std::vector<Halfedge>& generator : generators) {
        Vector<double> omega_i = buildClosedPrimalOneForm(generator);
        Vector<double> alpha_i = solver.solve(codiff1 * omega_i);
        Vector<double> dalpha_i = d0 * alpha_i;
        Vector<double> gamma_i = omega_i - dalpha_i;
        gammas.push_back(gamma_i);
    }



    //for (size_t i = 0; i < generators.size(); i++) {
    //    Vector<double> omega = buildClosedPrimalOneForm(generators[i]);
    //    Vector<double> dAlpha = hodgeDecomposition.computeExactComponent(omega);
    //    gammas.push_back(omega - dAlpha);
    //}


    return gammas; // placeholder
}