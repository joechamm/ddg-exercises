// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

typedef Eigen::Triplet<size_t> Triplet;

/*
 * Assign a unique index to each vertex, edge, and face of a mesh.
 * All elements are 0-indexed.
 *
 * Input: None. Access geometry via the member variable <geometry>, and pointer to the mesh via <mesh>.
 * Returns: None.
 */
void SimplicialComplexOperators::assignElementIndices() {

    // Needed to access geometry->vertexIndices, etc. as cached quantities.
    // Not needed if you're just using v->getIndex(), etc.
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    // You can set the index field of a vertex via geometry->vertexIndices[v], where v is a Vertex object (or an
    // integer). Similarly you can do edges and faces via geometry->edgeIndices, geometry->faceIndices, like so:
    size_t idx = 0;
    for (Vertex v : mesh->vertices()) {
        idx = geometry->vertexIndices[v];
    }

    for (Edge e : mesh->edges()) {
        idx = geometry->edgeIndices[e];
    }

    for (Face f : mesh->faces()) {
        idx = geometry->faceIndices[f];
    }

    // You can more easily get the indices of mesh elements using the function getIndex(), albeit less efficiently and
    // technically less safe (although you don't need to worry about it), like so:
    //
    //      v.getIndex()
    //
    // where v can be a Vertex, Edge, Face, Halfedge, etc. For example:

    for (Vertex v : mesh->vertices()) {
        idx = v.getIndex(); // == geometry->vertexIndices[v])
    }

    // Geometry Central already sets the indices for us, though, so this function is just here for demonstration.
    // You don't have to do anything :)
}

/*
 * Construct the unsigned vertex-edge adjacency matrix A0.
 *
 * Input:
 * Returns: The sparse vertex-edge adjacency matrix which gets stored in the global variable A0.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildVertexEdgeAdjacencyMatrix() const {

    // TODO
    // Note: You can build an Eigen sparse matrix from triplets, then return it as a Geometry Central SparseMatrix.
    // See <https://eigen.tuxfamily.org/dox/group__TutorialSparse.html> for documentation.
    size_t numVertices = mesh->nVertices();
    size_t numEdges = mesh->nEdges();
    std::vector<Triplet> triplets;
    for(Vertex v : mesh->vertices()) {
        for(Edge e : v.adjacentEdges()) {
            triplets.push_back(Triplet(e.getIndex(), v.getIndex(), 1));
        }
    }

    SparseMatrix<size_t> A0(numEdges, numVertices);
    A0.setFromTriplets(triplets.begin(), triplets.end());
    A0.makeCompressed();
    return A0;
}

/*
 * Construct the unsigned face-edge adjacency matrix A1.
 *
 * Input:
 * Returns: The sparse face-edge adjacency matrix which gets stored in the global variable A1.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildFaceEdgeAdjacencyMatrix() const {

    // TODO
    // Note: You can build an Eigen sparse matrix from triplets, then return it as a Geometry Central SparseMatrix.
    // See <https://eigen.tuxfamily.org/dox/group__TutorialSparse.html> for documentation.
    size_t numFaces = mesh->nFaces();
    size_t numEdges = mesh->nEdges();
    std::vector<Triplet> triplets;
    for(Face f : mesh->faces()) {
        for(Edge e : f.adjacentEdges()) {
            triplets.push_back(Triplet(f.getIndex(), e.getIndex(), 1));
        }
    }

    SparseMatrix<size_t> A1(numFaces, numEdges);
    A1.setFromTriplets(triplets.begin(), triplets.end());
    A1.makeCompressed();
    return A1;
}

/*
 * Construct a vector encoding the vertices in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |V|, where |V| = # of vertices in the mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildVertexVector(const MeshSubset& subset) const {

    // TODO

    // Note: You can build an Eigen sparse matrix from triplets, then return it as a Geometry Central SparseMatrix.
    // See <https://eigen.tuxfamily.org/dox/group__TutorialSparse.html> for documentation.
    // Note: You can build an Eigen sparse matrix from triplets, then return it as a Geometry Central SparseMatrix.
    // See <https://eigen.tuxfamily.org/dox/group__TutorialSparse.html> for documentation.
    size_t numVertices = mesh->nVertices();
    Vector<size_t> vertexVector(numVertices);
    vertexVector.setZero();
    for(size_t v : subset.vertices) {
        vertexVector[v] = 1;
    }
    return vertexVector;
}

/*
 * Construct a vector encoding the edges in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |E|, where |E| = # of edges in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildEdgeVector(const MeshSubset& subset) const {

    // TODO
    // Note: You can build an Eigen sparse matrix from triplets, then return it as a Geometry Central SparseMatrix.
    // See <https://eigen.tuxfamily.org/dox/group__TutorialSparse.html> for documentation.
    size_t numEdges = mesh->nEdges();
    Vector<size_t> edgeVector(numEdges);
    edgeVector.setZero();
    for(size_t e : subset.edges) {
        edgeVector[e] = 1;
    }
    return edgeVector;
}

/*
 * Construct a vector encoding the faces in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |F|, where |F| = # of faces in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildFaceVector(const MeshSubset& subset) const {

    // TODO
    // Note: You can build an Eigen sparse matrix from triplets, then return it as a Geometry Central SparseMatrix.
    // See <https://eigen.tuxfamily.org/dox/group__TutorialSparse.html> for documentation.
    size_t numFaces = mesh->nFaces();
    Vector<size_t> faceVector(numFaces);
    faceVector.setZero();
    for(size_t f : subset.faces) {
        faceVector[f] = 1;
    }
    return faceVector;
}

/*
 * Compute the simplicial star St(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The star of the given subset.
 */
MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {

    // TODO
    MeshSubset star = subset.deepCopy();
    Vector<size_t> vertexVector = buildVertexVector(subset);
    Vector<size_t> edgeVector = buildEdgeVector(subset);
    Vector<size_t> faceVector = buildFaceVector(subset);
    Vector<size_t> edgeIncidentVector = A0 * vertexVector;
    Vector<size_t> faceIncidentVector = A1 * (edgeIncidentVector + edgeVector);
    for(size_t i = 0; i < edgeIncidentVector.size(); i++) {
        if(edgeIncidentVector[i] > 0) {
            star.addEdge(i);
        }
    }
    for(size_t i = 0; i < faceIncidentVector.size(); i++) {
        if(faceIncidentVector[i] > 0) {
            star.addFace(i);
        }
    }

    return star;
}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {

    // TODO
    MeshSubset closure = subset.deepCopy();
    Vector<size_t> vertexVector = buildVertexVector(subset);
    Vector<size_t> edgeVector = buildEdgeVector(subset);
    Vector<size_t> faceVector = buildFaceVector(subset);
    Vector<size_t> edgeIncidentVector = A1.transpose() * faceVector + edgeVector;
    Vector<size_t> vertexIncidentVector = A0.transpose() * edgeIncidentVector + vertexVector;
    for(size_t i = 0; i < edgeIncidentVector.size(); i++) {
        if(edgeIncidentVector[i] > 0) {
            closure.addEdge(i);
        }
    }
    for(size_t i = 0; i < vertexIncidentVector.size(); i++) {
        if(vertexIncidentVector[i] > 0) {
            closure.addVertex(i);
        }
    }
    for(size_t i = 0; i < faceVector.size(); i++) {
        if(faceVector[i] > 0) {
            closure.addFace(i);
        }
    }
    return closure;
}

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {

    // TODO
    MeshSubset link = closure(star(subset));
    MeshSubset starClosure = star(closure(subset));
    link.deleteSubset(starClosure);
    return link;
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */
bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {

    // TODO
    return subset.equals(closure(subset)); // placeholder
}

/*
 * Check if the given subset S is a pure simplicial complex. If so, return the degree of the complex. Otherwise, return
 * -1.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: int representing the degree of the given complex (-1 if not pure)
 */
int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {

    // TODO
    if(!isComplex(subset)) {
        return -1;
    }
    Vector<size_t> faceVector = buildFaceVector(subset);
    Vector<size_t> edgeVector = buildEdgeVector(subset);
    Vector<size_t> vertexVector = buildVertexVector(subset);
    Vector<size_t> edgeCountVector = edgeVector.transpose() * A0;
    Vector<size_t> faceCountVector = faceVector.transpose() * A1;
    size_t pureVertices = (vertexVector.array() * edgeCountVector.array() > 0).count();
    size_t pureEdges = (edgeVector.array() * faceCountVector.array() > 0).count();
    bool allVerticesPure = subset.vertices.size() == pureVertices;
    bool allEdgesPure = subset.edges.size() == pureEdges;
    if(subset.faces.size()) return allEdgesPure && allVerticesPure ? 2 : -1;
    if(subset.edges.size()) return allVerticesPure ? 1 : -1;
    return 0;
}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The boundary of the given subset.
 */
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {

    // TODO
    MeshSubset boundary;
    int degree = isPureComplex(subset);
    if(degree == -1) {
        return boundary;
    }

    size_t numEdges = mesh->nEdges();
    size_t numVertices = mesh->nVertices();

    Vector<size_t> faceVector = buildFaceVector(subset);
    Vector<size_t> edgeVector = buildEdgeVector(subset);
    Vector<size_t> vertexVector = buildVertexVector(subset);
    Vector<size_t> countFaces = faceVector.transpose() * A1;
    Vector<size_t> countEdges = edgeVector.transpose() * A0;
    auto onceVertices = vertexVector.array() * countEdges.array() == 1;
    auto onceEdges = edgeVector.array() * countFaces.array() == 1;

    if(degree == 1) {
        for(size_t i = 0; i < numVertices; ++i) {
            if(onceVertices[i]) {
                boundary.addVertex(i);
            }
        }
    }

    if(degree == 2) {
        for(size_t i = 0; i < numEdges; ++i) {
            if(onceEdges[i]) {
                boundary.addEdge(i);
            }
        }
    }


    return closure(boundary);
}