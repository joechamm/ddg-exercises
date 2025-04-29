// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

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
        std::cout << "Vertex index: " << idx << std::endl;
        // Alternatively, you can use the following syntax to get the index of a vertex:
    }

    for (Edge e : mesh->edges()) {
        idx = geometry->edgeIndices[e];
        std::cout << "Edge index: " << idx << std::endl;
    }

    for (Face f : mesh->faces()) {
        idx = geometry->faceIndices[f];
        std::cout << "Face index: " << idx << std::endl;
    }

    // You can more easily get the indices of mesh elements using the function getIndex(), albeit less efficiently and
    // technically less safe (although you don't need to worry about it), like so:
    //
    //      v.getIndex()
    //
    // where v can be a Vertex, Edge, Face, Halfedge, etc. For example:

    for (Vertex v : mesh->vertices()) {
        idx = v.getIndex(); // == geometry->vertexIndices[v])
        std::cout << "Vertex index: " << idx << std::endl;
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
    size_t numVerts = mesh->nVertices();
    size_t numEdges = mesh->nEdges();
    std::vector<Eigen::Triplet<size_t>> triplets;
    for (Vertex v : mesh->vertices()) {
        size_t vIndex = v.getIndex();
        for (Edge e : v.adjacentEdges()) {
            size_t eIndex = e.getIndex();
            // The vertex index is the column and the edge index is the row
            triplets.push_back(Eigen::Triplet<size_t>(vIndex, eIndex, 1)); // Per Eigen documentation, the triplet constructor is (row, col, value), where value is 1 for non-zero entry in sparse matrix
        }
    }

    // SparseMatrix constructor is (rows, cols) 
    SparseMatrix<size_t> vertEdgeAdjacencyMatrix(numEdges, numVerts);
    vertEdgeAdjacencyMatrix.setFromTriplets(triplets.begin(), triplets.end());
    return vertEdgeAdjacencyMatrix;
}

/*
 * Construct the unsigned face-edge adjacency matrix A1.
 *
 * Input:
 * Returns: The sparse face-edge adjacency matrix which gets stored in the global variable A1.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildFaceEdgeAdjacencyMatrix() const {

    // TODO
    size_t numEdges = mesh->nEdges();
    size_t numFaces = mesh->nFaces();
    std::vector<Eigen::Triplet<size_t>> triplets;
    for (Face f : mesh->faces()) {
        size_t fIndex = f.getIndex();
        for (Edge e : f.adjacentEdges()) {
            size_t eIndex = e.getIndex();
            // The triangle face index is the row and edge index is the column
            triplets.push_back(Eigen::Triplet<size_t>(fIndex, eIndex, 1));  // Per Eigen documentation, the triplet constructor is (row, col, value), where value is 1 for non-zero entry in sparse matrix
        }
    }

    // SparseMatrix constructor is (rows, cols)
    SparseMatrix<size_t> faceEdgeAdjacencyMatrix(numFaces, numEdges);
    faceEdgeAdjacencyMatrix.setFromTriplets(triplets.begin(), triplets.end());
    return faceEdgeAdjacencyMatrix;
}

/*
 * Construct a vector encoding the vertices in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |V|, where |V| = # of vertices in the mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildVertexVector(const MeshSubset& subset) const {

    // TODO
    Vector<size_t> vertexVector(mesh->nVertices());
    vertexVector.setZero();
    for (Vertex v : mesh->vertices()) {

        if (subset.vertices.find(v.getIndex()) != subset.vertices.end()) {
            vertexVector[v.getIndex()] = 1; // Mark the vertex as selected
        }
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
    Vector<size_t> edgeVector(mesh->nEdges());
    edgeVector.setZero();
    for (Edge e : mesh->edges()) {
        if (subset.edges.find(e.getIndex()) != subset.edges.end()) {
            edgeVector[e.getIndex()] = 1; // Mark the edge as selected
        }
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
    Vector<size_t> faceVector(mesh->nFaces());
    faceVector.setZero();
    for (Face f : mesh->faces()) {
        if (subset.faces.find(f.getIndex()) != subset.faces.end()) {
            faceVector[f.getIndex()] = 1; // Mark the face as selected
        }
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
    // Note: The star of a vertex is the set of all edges and faces that are incident to it.
    // The star of an edge is the set of all faces that are incident to it.
    // The star of a face is the set of all edges and vertices that are incident to it.
    // The star of a subset of simplices is the union of the stars of each simplex in the subset.
    // Note: The star of a vertex is the set of all edges and faces that are incident to it.
    // Check if A0 or A1 has non-zero entries, then check if the subset is empty
    try {
        if (A0.nonZeros() == 0) {
            throw std::runtime_error("Error: Adjacency matrix A0 is empty.");
        } else if (A1.nonZeros() == 0) {
            throw std::runtime_error("Error: Adjacency matrix A1 is empty.");
        }
        if (subset.vertices.empty() && subset.edges.empty() && subset.faces.empty()) {
            throw std::runtime_error("Error: Subset is empty.");
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return subset; // placeholder
    }

    // first gather the star of the vertices in the subset
    // Use the adjacency matrix A0 to get the edges associated with the vertices in the subset
    // A0 is Column-major, so we'll to traverse over the columns of A0
    // The cth column of A0 corresponds to the cth vertex, and the rth row of A0 corresponds to the rth edge
    // The non-zero entries in the cth column of A0 correspond to the edges incident to the cth vertex

    // The star of the subset is the union of the stars of each simplex in the subset
    std::set<size_t> starVertices;
    std::set<size_t> starEdges;
    std::set<size_t> starFaces;

    // First we'll iterate over all the vertices in the subset, and use the adjacency matrices A0 and A1 to get the
    // edges and faces.
    for (size_t vIndex : subset.vertices) {
        starVertices.insert(vIndex);
        // iterate over the vIndex column of A0. The non-zero entries in this column correspond to the edges incident to
        // the vertex
        for (SparseMatrix<size_t>::InnerIterator it0(A0, vIndex); it0; ++it0) {
            // the row of our non-zero element is the index of an edge containing the vIndex vertex
            size_t eIndex = it0.row();
            starEdges.insert(eIndex);
            // now add the faces that contain this edge
            for (SparseMatrix<size_t>::InnerIterator it1(A1, eIndex); it1; ++it1) {
                // the row of our non-zero element is the index of a face containing the edge
                size_t fIndex = it1.row();
                starFaces.insert(fIndex);
            }
        }        
    }

    // Next we'll iterate over all the edges in the subset, and use the adjacency matrix A1 to get the faces
    for (size_t eIndex : subset.edges) {
        starEdges.insert(eIndex);
        for (SparseMatrix<size_t>::InnerIterator it1(A1, eIndex); it1; ++it1) {
        
        // the row of our non-zero element is the index of a face containing the edge
            size_t fIndex = it1.row();
            starFaces.insert(fIndex);
        }
    }

    // Finally, we'll just iterate over all the faces in the subset and add them to the star
    for (size_t fIndex : subset.faces) {
        starFaces.insert(fIndex);
    }

    MeshSubset starSubset(starVertices, starEdges, starFaces);

    return starSubset;
}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {

    // TODO
    // The closure of S is the smallest simplicial subcomplex that contains each simplex in S.
    // The closure of a face is the is the face, its edges, and vertices.
    // The closure of an edge is the edge and its vertices.
    // The closure of a vertex is the vertex itself.
    
    // we'll start by just copying the subset's vertices, edges, and faces
    MeshSubset closureSubset = subset.deepCopy();
    
    // Then we'll iterate over all the edges in the subset and add their vertices
    for (size_t eIndex : subset.edges) {
        Edge e = mesh->edge(eIndex);
        Vertex v0 = e.firstVertex();
        Vertex v1 = e.secondVertex();
        closureSubset.vertices.insert(v0.getIndex());
        closureSubset.vertices.insert(v1.getIndex());
    }

    // Next we'll iterate over all the faces in the subset and add their edges and vertices
    for (size_t fIndex : subset.faces) {
        Face f = mesh->face(fIndex);
        for (Edge e : f.adjacentEdges()) {
            closureSubset.edges.insert(e.getIndex());
            Vertex v0 = e.firstVertex();
            Vertex v1 = e.secondVertex();
            closureSubset.vertices.insert(v0.getIndex());
            closureSubset.vertices.insert(v1.getIndex());
        }
    }
    
    return closureSubset;
}

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {

    // TODO
    return subset; // placeholder
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */
bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {

    // TODO
    return false; // placeholder
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
    return -1; // placeholder
}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The boundary of the given subset.
 */
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {

    // TODO
    return subset; // placeholder
}