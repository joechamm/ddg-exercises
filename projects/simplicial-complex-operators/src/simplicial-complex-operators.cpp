// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

typedef Eigen::Triplet<size_t> T;

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
    std::vector<T> triplets;
    for (Vertex v : mesh->vertices()) {
        size_t vIndex = v.getIndex();
        Eigen::Index col = (Eigen::Index)vIndex;
        for (Edge e : v.adjacentEdges()) {
            size_t eIndex = e.getIndex();

            Eigen::Index row = (Eigen::Index)eIndex;
            // The vertex index is the column and the edge index is the row
            triplets.push_back(T(row, col, (size_t)1)); // Per Eigen documentation, the triplet constructor is (row, col, value), where value is 1 for non-zero entry in sparse matrix
        }
    }

    // SparseMatrix constructor is (rows, cols)
    Eigen::Index rows = (Eigen::Index)numEdges;
    Eigen::Index cols = (Eigen::Index)numVerts;
    SparseMatrix<size_t> vertEdgeAdjacencyMatrix(rows, cols);
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
    std::vector<T> triplets;
    for (Face f : mesh->faces()) {
        size_t fIndex = f.getIndex();
        Eigen::Index row = (Eigen::Index)fIndex;
        for (Edge e : f.adjacentEdges()) {
            size_t eIndex = e.getIndex();
            Eigen::Index col = (Eigen::Index)eIndex;
            // The triangle face index is the row and edge index is the column
            triplets.push_back(T(row, col, (size_t)1)); // Per Eigen documentation, the triplet constructor is (row, col, value), where value is 1 for non-zero entry in sparse matrix
        }
    }

    // SparseMatrix constructor is (rows, cols)
    Eigen::Index rows = (Eigen::Index)numFaces;
    Eigen::Index cols = (Eigen::Index)numEdges;
    SparseMatrix<size_t> faceEdgeAdjacencyMatrix(rows, cols);
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
    size_t numVerts = mesh->nVertices();
    Vector<size_t> vertexVector(numVerts); // Create a vector of size |V|, where |V| = # of vertices in mesh   
        
    vertexVector.setZero();
    for (size_t v : subset.vertices) {
        vertexVector[v] = 1; // Mark the vertex as selected
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
    size_t numEdges = mesh->nEdges();
    Vector<size_t> edgeVector(numEdges); // Create a vector of size |E|, where |E| = # of edges in mesh

    edgeVector.setZero();
    for (size_t e : subset.edges) {
        edgeVector[e] = 1; // Mark the edge as selected
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
    size_t numFaces = mesh->nFaces();
    Vector<size_t> faceVector(numFaces);
    faceVector.setZero();

    for (size_t f : subset.faces) {
        faceVector[f] = 1; // Mark the face as selected
    }

    return faceVector;
}

/*
 * Compute the simplicial star St(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The star of the given subset.
 */
//MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {
//
//    // TODO
//    // Note: The star of a vertex is the set of all edges and faces that are incident to it.
//    // The star of an edge is the set of all faces that are incident to it.
//    // The star of a face is the set of all edges and vertices that are incident to it.
//    // The star of a subset of simplices is the union of the stars of each simplex in the subset.
//    // Note: The star of a vertex is the set of all edges and faces that are incident to it.
//    // Check if A0 or A1 has non-zero entries, then check if the subset is empty
//    try {
//        if (A0.nonZeros() == 0) {
//            throw std::runtime_error("Error: Adjacency matrix A0 is empty.");
//        } else if (A1.nonZeros() == 0) {
//            throw std::runtime_error("Error: Adjacency matrix A1 is empty.");
//        }
//        if (subset.vertices.empty() && subset.edges.empty() && subset.faces.empty()) {
//            throw std::runtime_error("Error: Subset is empty.");
//        }
//    } catch (const std::exception& e) {
//        std::cerr << e.what() << std::endl;
//        return subset; // placeholder
//    }
//
//    // first gather the star of the vertices in the subset
//    // Use the adjacency matrix A0 to get the edges associated with the vertices in the subset
//    // A0 is Column-major, so we'll to traverse over the columns of A0
//    // The cth column of A0 corresponds to the cth vertex, and the rth row of A0 corresponds to the rth edge
//    // The non-zero entries in the cth column of A0 correspond to the edges incident to the cth vertex
//
//    // The star of the subset is the union of the stars of each simplex in the subset
//    std::set<size_t> starVertices;
//    std::set<size_t> starEdges;
//    std::set<size_t> starFaces;
//
//    // First we'll iterate over all the vertices in the subset, and use the adjacency matrices A0 and A1 to get the
//    // edges and faces.
//    for (size_t vIndex : subset.vertices) {
//        starVertices.insert(vIndex);
//        // iterate over the vIndex column of A0. The non-zero entries in this column correspond to the edges incident to
//        // the vertex
//        for (SparseMatrix<size_t>::InnerIterator it0(A0, vIndex); it0; ++it0) {
//            // the row of our non-zero element is the index of an edge containing the vIndex vertex
//            size_t eIndex = it0.row();
//            starEdges.insert(eIndex);
//            // now add the faces that contain this edge
//            for (SparseMatrix<size_t>::InnerIterator it1(A1, eIndex); it1; ++it1) {
//                // the row of our non-zero element is the index of a face containing the edge
//                size_t fIndex = it1.row();
//                starFaces.insert(fIndex);
//            }
//        }
//    }
//
//    // Next we'll iterate over all the edges in the subset, and use the adjacency matrix A1 to get the faces
//    for (size_t eIndex : subset.edges) {
//        starEdges.insert(eIndex);
//        for (SparseMatrix<size_t>::InnerIterator it1(A1, eIndex); it1; ++it1) {
//
//            // the row of our non-zero element is the index of a face containing the edge
//            size_t fIndex = it1.row();
//            starFaces.insert(fIndex);
//        }
//    }
//
//    // Finally, we'll just iterate over all the faces in the subset and add them to the star
//    for (size_t fIndex : subset.faces) {
//        starFaces.insert(fIndex);
//    }
//
//    MeshSubset starSubset(starVertices, starEdges, starFaces);
//
//    return starSubset;
//}
MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {
    // Start by copying the subset. We'll add the star vertices, edges, and faces to this copy.
    MeshSubset star = subset.deepCopy();
    // Get some useful information about the mesh.
    size_t numVerts = mesh->nVertices();
    size_t numEdges = mesh->nEdges();
    size_t numFaces = mesh->nFaces();
    // Build the vertex, edge, and face vectors for use in computing adjacency.
    Vector<size_t> vertexVector = buildVertexVector(subset);
    Vector<size_t> edgeVector = buildEdgeVector(subset);
    Vector<size_t> faceVector = buildFaceVector(subset);
    // Convert the vertices in subset S to an index vector we can multiply by A0 to get the edges that contain vertices
    // in S. The resultant vector will be an n dimension vector, where n = |E|. The the ith component will be the number
    // of vertices in S that are incident to the ith edge. We can then add the non-zero entries in this vector to the
    // star edges.
    //for (size_t i = 0; i < numVerts; ++i) {
    //    star.addVertex(i); // Add the vertex to the star vertices
    //}
    //Vector<size_t> edgeIncidentVector = A0 * vertexVector + edgeVector; // A0 is the vertex-edge adjacency matrix
    Vector<size_t> edgeIncidentVector = A0 * vertexVector;
    
    // Iterate over the edgeIncidentVector and add the non-zero entries to the star edges
    for (size_t i = 0; i < numEdges; i++) {
        if (edgeIncidentVector[i] > 0) {
            star.addEdge(i); // Add the edge to the star edges
        }
    }
    // Now we can do the same for the faces. The resultant vector will be an n dimension vector, where n = |F|. The the
    // ith component will be the number of edges in S that are incident to the ith face. We can then add the non-zero
    // entries in this vector to the star faces.
 //   Vector<size_t> faceIncidentVector = A1 * edgeVector + faceVector; // A1 is the face-edge adjacency matrix
    Vector<size_t> faceIncidentVector = A1 * (edgeIncidentVector + edgeVector);

    // Iterate over the faceIncidentVector and add the non-zero entries to the star faces
    for (size_t i = 0; i < numFaces; i++) {
        if (faceIncidentVector[i] > 0) {
            star.addFace(i); // Add the face to the star faces
        }
    }

    // Debugging - print star
    star.printVertices();
    star.printEdges();
    star.printFaces();
 
    return star; // Return the star of the subset
}

/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
//MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {
//
//    // TODO
//    // The closure of S is the smallest simplicial subcomplex that contains each simplex in S.
//    // The closure of a face is the is the face, its edges, and vertices.
//    // The closure of an edge is the edge and its vertices.
//    // The closure of a vertex is the vertex itself.
//
//    // we'll start by just copying the subset's vertices, edges, and faces
//    MeshSubset closureSubset = subset.deepCopy();
//
//    // Then we'll iterate over all the edges in the subset and add their vertices
//    for (size_t eIndex : subset.edges) {
//        Edge e = mesh->edge(eIndex);
//        Vertex v0 = e.firstVertex();
//        Vertex v1 = e.secondVertex();
//        closureSubset.vertices.insert(v0.getIndex());
//        closureSubset.vertices.insert(v1.getIndex());
//    }
//
//    // Next we'll iterate over all the faces in the subset and add their edges and vertices
//    for (size_t fIndex : subset.faces) {
//        Face f = mesh->face(fIndex);
//        for (Edge e : f.adjacentEdges()) {
//            closureSubset.edges.insert(e.getIndex());
//            Vertex v0 = e.firstVertex();
//            Vertex v1 = e.secondVertex();
//            closureSubset.vertices.insert(v0.getIndex());
//            closureSubset.vertices.insert(v1.getIndex());
//        }
//    }
//
//    return closureSubset;
//}

MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {
    // Get some useful information about the mesh.
    size_t numVerts = mesh->nVertices();
    size_t numEdges = mesh->nEdges();
    size_t numFaces = mesh->nFaces();
    // Build the vertex, edge, and face vectors for use in computing adjacency.
    Vector<size_t> vertexVector = buildVertexVector(subset);
    Vector<size_t> edgeVector = buildEdgeVector(subset);
    Vector<size_t> faceVector = buildFaceVector(subset);

    // Multiplying the transpose of A1 by the face vector will give us a vector of size |E|, where
    // the ith component will be the number of faces in S that are incident to the ith edge. Adding this
    // to edgeVector will give us all the edges in the closure.
    Vector<size_t> edgeIncidentVector = A1.transpose() * faceVector + edgeVector; // A1 is the face-edge adjacency matrix
    // We can use the same approach for the vertices. Multiplying the transpose of A0 by the edge vector will give us a
    // vector of size |V|, where the ith component will be the number of edges in S that are incident to the ith vertex.
    // Adding this to vertexVector will give us all the vertices in the closure.
    Vector<size_t> vertexIncidentVector =
        A0.transpose() * edgeIncidentVector + vertexVector; // A0 is the vertex-edge adjacency matrix

    // Now we'll build the closure subset
    MeshSubset closureSubset;
    // Add the vertices in the closure
    for (size_t i = 0; i < numVerts; i++) {
        if (vertexIncidentVector[i] > 0) {
            closureSubset.addVertex(i); // Add the vertex to the closure vertices
        }
    }

    // Add the edges in the closure
    for (size_t i = 0; i < numEdges; i++) {
        if (edgeIncidentVector[i] > 0) {
            closureSubset.addEdge(i); // Add the edge to the closure edges
        }
    }

    // Add the faces in the closure
    for (size_t i = 0; i < numFaces; i++) {
        if (faceVector[i] > 0) {
            closureSubset.addFace(i); // Add the face to the closure faces
        }
    }
    return closureSubset; // Return the closure of the subset
}

/*
 * Compute the closed star of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the star of the given subset.
 */ 
MeshSubset SimplicialComplexOperators::closedStar(const MeshSubset& subset) const {
    return closure(star(subset));   
}


/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
//MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {
//
//    // TODO
//    // The link of S is the set theoretic difference of the closure of the star of S minus the star of the closure of S
//    // Lk(S) = Cl(St(S)) \ St(St(Cl(S)))
//    // Make sure that subset is not empty
//    try {
//        if (subset.vertices.empty() && subset.edges.empty() && subset.faces.empty()) {
//            throw std::runtime_error("Error: Subset is empty.");
//        }
//    } catch (const std::exception& e) {
//        std::cerr << e.what() << std::endl;
//        return subset;
//    }
//
//    // start by getting the star and closure of S
//    MeshSubset stS = star(subset);    // St(S)
//    MeshSubset clS = closure(subset); // Cl(S)
//    // then get the closure of the star
//    MeshSubset clStS = closure(stS); // Cl(St(S))
//    // and the star of the closure
//    MeshSubset stClS = star(clS); // St(Cl(S))
//
//    // declare the link subset
//    MeshSubset linkS;
//    // use the set difference of the closure of the star and the star of the closure only the vertices, edges, and faces
//    std::set_difference(clStS.vertices.begin(), clStS.vertices.end(), stClS.vertices.begin(), stClS.vertices.end(),
//                        std::inserter(linkS.vertices, linkS.vertices.end()));
//    std::set_difference(clStS.edges.begin(), clStS.edges.end(), stClS.edges.begin(), stClS.edges.end(),
//                        std::inserter(linkS.edges, linkS.edges.end()));
//    std::set_difference(clStS.faces.begin(), clStS.faces.end(), stClS.faces.begin(), stClS.faces.end(),
//                        std::inserter(linkS.faces, linkS.faces.end()));
//
//    return linkS;
//}
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {
    // The link of S is the set theoretic difference of the closure of the star of S minus the star of the closure of S.
    // We'll start with Cl(St(S)), then remove the star of the closure of S.
    MeshSubset link = closedStar(subset); // Cl(St(S))
    MeshSubset starClosure = star(closure(subset)); // St(Cl(S))
    // Now we'll remove the star of the closure from the link
    link.deleteSubset(starClosure);
    return link;
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */
//bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {
//
//    // TODO
//    // Check if subset is empty
//    if (subset.vertices.empty() && subset.edges.empty() && subset.faces.empty()) {
//        return false; // The subset is empty
//    }
//    // A simplicial complex is a set of simplices such that every face of a simplex in the set is also in the set.
//    // We'll use the closure operator to check if the subset is a simplicial complex.
//    // If the closure of the subset is equal to the subset, then it is a simplicial complex.
//    MeshSubset closureSubset = closure(subset);
//    // Check if the closure of the subset is equal to the subset
//    if (closureSubset.vertices == subset.vertices && closureSubset.edges == subset.edges &&
//        closureSubset.faces == subset.faces) {
//        return true; // The subset is a simplicial complex
//    }
//
//    return false; // placeholder
//}
bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {
    return subset.equals(closure(subset));
}


/*
 * Check if the given subset S is a pure simplicial complex. If so, return the degree of the complex. Otherwise, return
 * -1.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: int representing the degree of the given complex (-1 if not pure)
 */
//int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {
//
//    // TODO
//    // Check if subset is empty
//    if (subset.vertices.empty() && subset.edges.empty() && subset.faces.empty()) {
//        return -1; // The subset is empty
//    }
//    // A complex K is a pure k-simplicial complex if every simplex in K is contained in some simplex of degree k
//    // (possibly itself). Let's grab the closure first and check and if the subset if a complex.
//    //
//    MeshSubset closureSubset = closure(subset);
//    // Check if the closure of the subset is equal to the subset
//    if (closureSubset.vertices != subset.vertices || closureSubset.edges != subset.edges ||
//        closureSubset.faces != subset.faces) {
//        return -1; // The subset is not a simplicial complex
//    }
//
//    if (!subset.faces.empty()) {
//        // Check if every edge is contained in a face in the subset.
//        for (size_t eIndex : subset.edges) {
//            // iterate over the eIndex column of A1. The non-zero entries in this column correspond to the faces
//            // incident to the edge
//            for (SparseMatrix<size_t>::InnerIterator it1(A1, eIndex); it1; ++it1) {
//                // the row of our non-zero element is the index of a face containing the edge
//                size_t fIndex = it1.row();
//                // check if the face is contained in the subset
//                if (subset.faces.find(fIndex) == subset.faces.end()) {
//                    return -1; // The subset is not a pure complex
//                }
//            }
//        }
//
//        // Now check if every vertex is contained in an edge in the subset
//        for (size_t vIndex : subset.vertices) {
//            // iterate over the vIndex column of A0. The non-zero entries in this column correspond to the edges
//            // incident to the vertex
//            for (SparseMatrix<size_t>::InnerIterator it0(A0, vIndex); it0; ++it0) {
//                // the row of our non-zero element is the index of an edge containing the vertex
//                size_t eIndex = it0.row();
//                // check if the edge is contained in the subset
//                if (subset.edges.find(eIndex) == subset.edges.end()) {
//                    return -1; // The subset is not a pure complex
//                }
//            }
//        }
//
//        // Every edge is contained in a face in the subset, and every vertex is contained in an edge in the subset.
//        // Degree of the complex is 2
//        return 2;
//    }
//
//    if (!subset.edges.empty()) {
//        // Check if every vertex is contained in an edge in the subset
//        for (size_t vIndex : subset.vertices) {
//            // iterate over the vIndex column of A0. The non-zero entries in this column correspond to the edges
//            // incident to the vertex
//            for (SparseMatrix<size_t>::InnerIterator it0(A0, vIndex); it0; ++it0) {
//                // the row of our non-zero element is the index of an edge containing the vertex
//                size_t eIndex = it0.row();
//                // check if the edge is contained in the subset
//                if (subset.edges.find(eIndex) == subset.edges.end()) {
//                    return -1; // The subset is not a pure complex
//                }
//            }
//        }
//        // Every vertex is contained in an edge in the subset.
//        // Degree of the complex is 1
//        return 1;
//    }
//
//    // If we reach here, then the subset is a pure complex of degree 0
//    return 0;
//
//    //// now add the faces that contain this edge
//    // for (SparseMatrix<size_t>::InnerIterator it1(A1, eIndex); it1; ++it1) {
//    //     // the row of our non-zero element is the index of a face containing the edge
//    //     size_t fIndex = it1.row();
//    //     starFaces.insert(fIndex);
//    // }
//
//    ////if (!subset.faces.empty()) { // degree = 2
//    ////    // Check if every edge is contained in a face in the subset.
//    ////    // If we take a vector V of dimension |E|, with 1 in the ith component if edge i is in S and a 0 otherwise,
//    ///we can multiply V by A1 since A1 contains |E| columns. /    // The result will be a vector of dimension |F|, with
//    ///the jth component equal to the number of edges in S contained in the jth face. /    Vector<size_t> edgeVector =
//    ///buildEdgeVector(subset); /    Vector<size_t> faceIncidentVector = A1 * edgeVector; // A1 is the face-edge
//    ///adjacency matrix /    // For each face f with index j in our complex, check the jth component of the
//    ///faceIncidentVector.
//
//    ////}
//
//
//    // int degree = 2;
//    //// If subset faces is empty, then the degree is 1, and if edges is empty and vertices is non-empty, then the
//    ///degree / is 0.
//    // if (subset.faces.empty()) {
//    //     if (subset.edges.empty()) {
//    //         // degree is 0 since we already checked that the subset is not empty
//    //         degree = 0;
//    //     } else {
//    //         // degree is 1 since we already checked that the subset is not empty
//    //         degree = 1;
//    //     }
//    // }
//
//    //// Check if the subset is a pure complex.
//    //// If degree = 2, then we check if each vertex is contained in an edge in subset, and if each edge is contained in
//    ///a face / in subset.
//    // if (degree == 2) {
//    //     for (size_t vIndex : subset.vertices) {
//    //         // check if the vertex is contained in an edge in the subset
//    //         if (subset.edges.find(vIndex) == subset.edges.end()) {
//    //             return -1; // The subset is not a pure complex
//    //         }
//    //     }
//    //     for (size_t eIndex : subset.edges) {
//    //         // check if the edge is contained in a face in the subset
//    //         if (subset.faces.find(eIndex) == subset.faces.end()) {
//    //             return -1; // The subset is not a pure complex
//    //         }
//    //     }
//    // }
//
//    //// If degree = 1, we just need to check if each vertex is contained in an edge in subset.
//
//
//    // return degree;
//}

int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {
    if (!isComplex(subset)) {    
        return -1; // The subset is not a pure complex
    }

    // Try using andy1li's approach
    Vector<size_t> vertexVector = buildVertexVector(subset);
    Vector<size_t> edgeVector = buildEdgeVector(subset);
    Vector<size_t> faceVector = buildFaceVector(subset);
    Vector<size_t> countEdges = edgeVector.transpose() * A0;
    Vector<size_t> countFaces = faceVector.transpose() * A1;

    size_t pureVertices = (vertexVector.array() * countEdges.array() > 0).count();
    size_t pureEdges = (edgeVector.array() * countFaces.array() > 0).count();

    bool allVerticesPure = subset.vertices.size() == pureVertices;
    bool allEdgesPure = subset.edges.size() == pureEdges;

    if (subset.faces.size()) return allEdgesPure && allVerticesPure ? 2 : -1;
    if (subset.edges.size()) return allVerticesPure ? 1 : -1;
    return 0;



    //// Get some useful information about the mesh.
    //size_t numVerts = mesh->nVertices();
    //size_t numEdges = mesh->nEdges();
    //size_t numFaces = mesh->nFaces();
    //// Build the vertex, edge, and face vectors for use in computing adjacency.
    //Vector<size_t> vertexVector = buildVertexVector(subset);
    //Vector<size_t> edgeVector = buildEdgeVector(subset);
    //Vector<size_t> faceVector = buildFaceVector(subset);

    //// Multiplying the transpose of A1 by the face vector will give us a vector of size |E|, where
    //// the ith component will be the number of faces in S that are incident to the ith edge. Adding this
    //// to edgeVector will give us all the edges in the closure.
    //Vector<size_t> edgeIncidentVector =
    //    A1.transpose() * faceVector; // A1 is the face-edge adjacency matrix
    //// We can use the same approach for the vertices. Multiplying the transpose of A0 by the edge vector will give us a
    //// vector of size |V|, where the ith component will be the number of edges in S that are incident to the ith vertex.
    //// Adding this to vertexVector will give us all the vertices in the closure.
    //Vector<size_t> vertexIncidentVector =
    //    A0.transpose() * edgeIncidentVector + vertexVector; // A0 is the vertex-edge adjacency matrix

    //// Now we'll build the closure subset
    //MeshSubset closureSubset;
    //// Add the vertices in the closure
    //for (size_t i = 0; i < numVerts; i++) {
    //    if (vertexIncidentVector[i] > 0) {
    //        closureSubset.addVertex(i); // Add the vertex to the closure vertices
    //    }
    //}

    //// Add the edges in the closure
    //for (size_t i = 0; i < numEdges; i++) {
    //    if (edgeIncidentVector[i] > 0) {
    //        closureSubset.addEdge(i); // Add the edge to the closure edges
    //    }
    //}

    //// Add the faces in the closure
    //for (size_t i = 0; i < numFaces; i++) {
    //    if (faceVector[i] > 0) {
    //        closureSubset.addFace(i); // Add the face to the closure faces
    //    }
    //}
    //return closureSubset; // Return the closure of the subset

}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The boundary of the given subset.
 */
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {

    // TODO
    // The boundary, bd(K'), of a pure k-subcomplex K' of simplicial complex K is the closure of the set of all simplices, s, that are proper faces of exactly one simplex of K'. 
    MeshSubset boundary;
    int degree = isPureComplex(subset);
    if (degree == -1) {
        return boundary; // The subset is not a pure complex
    }

    size_t numVertices = mesh->nVertices();
    size_t numEdges = mesh->nEdges();
    size_t numFaces = mesh->nFaces();

    Vector<size_t> vertexVector = buildVertexVector(subset);
    Vector<size_t> edgeVector = buildEdgeVector(subset);
    Vector<size_t> faceVector = buildFaceVector(subset);
    Vector<size_t> countEdges = edgeVector.transpose() * A0;
    Vector<size_t> countFaces = faceVector.transpose() * A1;
    auto onceVertices = vertexVector.array() * countEdges.array() == 1;
    auto onceEdges = edgeVector.array() * countFaces.array() == 1;

    if (degree == 1) {
        for (size_t i = 0; i < numVertices; ++i) {
            if (onceVertices[i]) {
                boundary.addVertex(i);
            }
        }
    }

    if (degree == 2) {
        for (size_t i = 0; i < numEdges; ++i) {
            if (onceEdges[i]) {
                boundary.addEdge(i);
            }
        }
    }

    return closure(boundary);
}