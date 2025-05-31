// Implement member functions for TreeCotree class.
#include "tree-cotree.h"
#include <queue>
#include <set>
#include <algorithm>

/*
    Primal Spanning Tree(T) - A spanning tree of the graph formed by the vertices and edges of the triangulation.
    It connects all the vertices of the mesh without forming any cycles.
 
    Dual Mesh - A graph where each vertex corresponds to a face(triangle) of the primal mesh, and an edge 
    exists between two dual vertices if their corresponding primal faces share an edge.
 
    Dual Spanning Tree(T⋆) - A spanning tree of the dual graph. It connects all the faces of the primal mesh 
    through their adjacencies without forming any cycles in the dual graph.
         
    Dual Edge(e_ij*) - An edge in the dual graph connecting two dual vertices(faces i and j of the primal mesh). 
    This dual edge corresponds to a primal edge e_ij that is shared by faces i and j.
 
    Generators of the Fundamental Group - A minimal set of loops whose homotopy classes can generate all other 
    loops on the surface up to homotopy. For a closed orientable surface of genus g, the fundamental group has 
    2g generators.
*/

/*
 * Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
TreeCotree::TreeCotree(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    mesh = inputMesh;
    geometry = inputGeo;
}

/*
 * Build a primal spanning tree on a mesh without boundary. More specifically, populate the member variable
 * <vertexParent>, which is a std::map that maps each vertex of the input mesh to its parent in the primal spanning
 * tree.
 *
 * Input:
 * Returns:
 */
void TreeCotree::buildPrimalSpanningTree() {
    Vertex root = *(mesh->vertices().begin()); // Get the first vertex as the root

    // Mark each vertex as its own parent
    for (Vertex v : mesh->vertices()) {
        vertexParent[v] = v;
    }

    std::queue<Vertex> vertexQueue; // Queue for BFS
    vertexQueue.push(root);         // Start BFS from the root
    while (!vertexQueue.empty()) {
        Vertex currentVertex = vertexQueue.front(); // Get the current vertex
        vertexQueue.pop();                          // Remove it from the queue
        for (Vertex neighbor : currentVertex.adjacentVertices()) {
            // Check if the neighbor has been visited yet
            if (vertexParent[neighbor] == neighbor && neighbor != root) {
                vertexParent[neighbor] = currentVertex; // Set the parent of the neighbor
                vertexQueue.push(neighbor);             // Add the neighbor to the queue
            }
        }
    }
}

//void TreeCotree::buildPrimalSpanningTree() {
//    std::queue<Vertex> vertexQueue;
//    Vertex root = *(mesh->vertices().begin());
//
//    // Initialze the queue for a Breadth-First Search.
//    vertexQueue.push(root);
//
//    // Keep track of visited vertices to avoid cycles.
//    std::set<Vertex> visitedVertices;
//    visitedVertices.insert(root);
//    vertexParent[root] = root; // Set the root's parent to itself
//
//    // Perform BFS to build the primal spanning tree.
//    while (!vertexQueue.empty()) {
//        Vertex currentVertex = vertexQueue.front();
//        vertexQueue.pop();
//        // Iterate through the negihbors of the current vertex.
//        for (Halfedge he : currentVertex.outgoingHalfedges()) {
//            Vertex neighbor = he.tipVertex();
//            // Check if the neighbor has already been visited.
//            if (visitedVertices.find(neighbor) == visitedVertices.end()) {
//                // Check if the edge is in the dual tree and add it to the tree if not.
//                if (!inDualSpanningCotree(he)) {
//                    // Neighbor has not been visited yet, so add it to the tree.
//                    visitedVertices.insert(neighbor);
//                    vertexParent[neighbor] = currentVertex; // Set parent of neighbor to currentVertex
//                    vertexQueue.push(neighbor);             // Add neighbor to the queue
//                }
//
//            }
//        }        
//    }
//}
    //void TreeCotree::buildPrimalSpanningTree() {
//
//    // TODO
//    std::queue<Vertex> vertexQueue;
//    Vertex root = *(mesh->vertices().begin());
//
//    vertexQueue.push(root);
//    vertexParent[root] = root;
//    while (!vertexQueue.empty()) {
//        Vertex v = vertexQueue.front();
//        vertexQueue.pop();
//        for (Halfedge he : v.outgoingHalfedges()) {
//            Vertex w = he.tipVertex();
//            if (!vertexParent.count(w)) {
//                vertexParent[w] = v; // Set parent of w to v
//                vertexQueue.push(w); // Add w to the queue
//            }
//        }
//    }
//}
/*
 * Check whether a halfedge is in the primal spanning tree.
 *
 * Input: A halfedge <he>
 * Returns: True if <he> is in the primal spanning tree, false otherwise.
 */
bool TreeCotree::inPrimalSpanningTree(Halfedge he) {
    Vertex v = he.vertex();
    Vertex w = he.tipVertex();

    // Check if the halfedge is in the primal spanning tree by checking if either vertex is a child of the other.
    if (vertexParent.count(v)) {
        if (vertexParent[v] == w) {
            return true; // v is a child of w
        }
    }

    if (vertexParent.count(w)) {
        if (vertexParent[w] == v) {
            return true; // w is a child of v
        }
    }
    
    // Return false if neither vertex is a child of the other.
    return false;
}

/*
 * Build a dual spanning tree on a mesh without boundary. More specificially, populate the member variable <faceParent>,
 * which is a std::map that maps each face of the input mesh to its parent in the dual spanning tree.
 *
 * Input:
 * Returns:
 */
void TreeCotree::buildDualSpanningCoTree() {
    Face root = *(mesh->faces().begin()); // Get the first face as the root

    // Mark each face as its own parent
    for (Face f : mesh->faces()) {
        faceParent[f] = f;
    }

    // Initialize the queue for a Breadth-First Search.
    std::queue<Face> faceQueue;
    faceQueue.push(root);
    while (!faceQueue.empty()) {
        Face currentFace = faceQueue.front();
        faceQueue.pop();
        // Iterate through the adjacent halfedges
        for (Halfedge he : currentFace.adjacentHalfedges()) {
            // Make sure we don't cross an edge in the primal spanning tree
            if (!inPrimalSpanningTree(he)) {
                Face neighborFace = he.twin().face();

                // Check if the neighbor face has already been visited.
                if (faceParent[neighborFace] == neighborFace && neighborFace != root) {
                    faceParent[neighborFace] = currentFace;
                    faceQueue.push(neighborFace); // Add the neighbor face to the queue
                }
            }
        }
    }
}

//void TreeCotree::buildDualSpanningCoTree() {
//    // Initialize the faceParent map, and initialize the queue with the root face for a Breadth-First Search.
//    std::queue<Face> faceQueue;
//    Face root = *(mesh->faces().begin());
//    faceQueue.push(root);
//    faceParent[root] = root;
//
//    // Keep track of visited faces to avoid cycles.
//    std::set<Face> visitedFaces;
//    visitedFaces.insert(root);
//
//    // Perform BFS to build the dual spanning tree.
//    while (!faceQueue.empty()) {
//        Face currentFace = faceQueue.front();
//        faceQueue.pop();
//
//        // Iterate through the neighbors of the current face.
//        for (Face neighbor : currentFace.adjacentFaces()) {
//            // Check if the neighbor has already been visited.
//            if (visitedFaces.find(neighbor) == visitedFaces.end()) {
//                // Neighbor has not been visited yet, so add it to the tree.
//                visitedFaces.insert(neighbor);
//                faceParent[neighbor] = currentFace; // Set parent of neighbor to currentFace
//                faceQueue.push(neighbor);           // Add neighbor to the queue                
//            }
//        }
//    }
//}

// void TreeCotree::buildDualSpanningCoTree() {
//
//     // TODO
//     std::queue<Face> faceQueue;
//     // Initialize the faceParent map, and initialize the queue with the root face for a Breadth-First Search.
//     Face root = *(mesh->faces().begin());
//     faceQueue.push(root);
//     faceParent[root] = root;
//
//     while (!faceQueue.empty()) {
//         Face f = faceQueue.front();
//         faceQueue.pop();
//         for (Halfedge he : f.adjacentHalfedges()) {
//             if (inPrimalSpanningTree(he)) {
//                 continue;
//             }
//
//             Face g = he.twin().face();
//             if (!faceParent.count(g)) {
//                 faceParent[g] = f; // Set parent of g to f
//                 faceQueue.push(g); // Add g to the queue
//             }
//         }
//
//     }
//
// }

/*
 * Check whether a halfedge is in the dual spanning tree.
 *
 * Input: A halfedge <he>
 * Returns: True if <he> is in the dual spanning tree, false otherwise.
 */
bool TreeCotree::inDualSpanningCotree(Halfedge he) {

    // Check if the halfedge is in the dual spanning tree by checking if either face is a child of the other.
    Face f = he.face();
    Face g = he.twin().face();

    if (faceParent.count(f)) {
        if (faceParent[f] == g) {
            return true; // f is a child of g
        }
    }

    if (faceParent.count(g)) {
        if (faceParent[g] == f) {
            return true; // g is a child of f
        }
    }

    // Return false if neither face is a child of the other.
    return false; 
}

/*
 * Returns a halfedge lying on the shared edge between face f and g.
 *
 * Input: Two adjacent faces <f> and <g>.
 * Returns: A halfedge lying on the shared edge between face f and g.
 */
Halfedge TreeCotree::sharedHalfedge(Face f, Face g) const {

    for (Halfedge he : f.adjacentHalfedges()) {
        if (he.twin().face() == g) {
            return he;
        }
    }
    // Should never get here!
    std::cerr << "Oops, TreeCotree::sharedHalfedge() received bad input." << std::endl;
    return f.halfedge();
}

/*
 * Compute the homology generators of the mesh.
 *
 * Input:
 * Returns:
 */
void TreeCotree::buildGenerators() {
    // For a closed surface of genus g, there are 2g generators.
    int g = mesh->genus();
    generators.clear();
//    generators.reserve(2 * g);

    // Build primal spanning tree T of primal edges.
    buildPrimalSpanningTree();

    // Build dual spanning tree T* of dual edges.
    buildDualSpanningCoTree();

    // collect dual edges that are neight in primal spanning tree nor in dual spanning cotree
    for (Edge e : mesh->edges()) {
        Halfedge he = e.halfedge();
        // Only consider edges that are not in the primal spanning tree or the dual spanning cotree.
        if (!inPrimalSpanningTree(he) && !inDualSpanningCotree(he)) {
            // trace faces back to root

            std::vector<Halfedge> tempGenerator1; // Set to store faces for the first part of the generator
//            std::set<Halfedge> tempGenerator1; // Set to store faces for the first part of the generator
            Face currentFace = he.face();
            while (faceParent[currentFace] != currentFace) {
                Face parentFace = faceParent[currentFace]; // Get the parent face
                Halfedge sharedHe = sharedHalfedge(currentFace, parentFace); // Get the halfedge shared with the parent face
                tempGenerator1.push_back(sharedHe);          // Add the halfedge to the set
//                tempGenerator1.insert(sharedHe); // Add the halfedge to the set
                currentFace = parentFace;       // Move to the parent face
            }

            std::vector<Halfedge> tempGenerator2; // Set to store faces for the second part of the generator
//            std::set<Halfedge> tempGenerator2; // Set to store faces for the second part of the generator
            currentFace = he.twin().face();    // Start from the twin face
            while (faceParent[currentFace] != currentFace) {
                Face parentFace = faceParent[currentFace]; // Get the parent face
                Halfedge sharedHe = sharedHalfedge(currentFace, parentFace); // Get the halfedge shared with the parent face
                tempGenerator2.push_back(sharedHe);             // Add the halfedge to the set
//                tempGenerator2.insert(sharedHe); // Add the halfedge to the set
                currentFace = parentFace;                                    // Move to the parent face
            }

            //// Remove common halfedges
            size_t m = tempGenerator1.size() - 1;
            size_t n = tempGenerator2.size() - 1;
            while (tempGenerator1[m] == tempGenerator2[n]) {
                m--;
                n--;
            }

            std::vector<Halfedge> generator; // Final generator
            generator.push_back(he);         // Add the original halfedge to the generator
            for (int i = 0; i <= m; i++) {
                Halfedge heTwin = tempGenerator1[i].twin();    // Get the twin halfedge
                generator.push_back(heTwin); // Add the halfedges from the first part of the generator
            }

            for (int i = n; i >= 0; i--) {
                generator.push_back(tempGenerator2[i]); // Add the halfedges from the second part of the generator
            }

            // remove the common halfedges
            //std::vector<Halfedge> generator; // Final generator
            //std::set_symmetric_difference(tempGenerator1.begin(), tempGenerator1.end(), tempGenerator2.begin(),
            //                              tempGenerator2.end(),
            //                              generator.begin()); // Compute the symmetric difference of the two sets
            
           

            generators.push_back(generator); // Add the generator to the list of generators

        }

    }
}

//void TreeCotree::buildGenerators() {
//
//    // For a closed surface of genus g, there are 2g generators.
//    int g = mesh->genus();
//    generators.clear();
//    generators.reserve(2 * g);
//
//    // Build spanning tree T of primal edges.
//    buildPrimalSpanningTree();
//
//    // Build spanning tree T* of dual edges.
//    buildDualSpanningCoTree();
//
//    // For each dual edge e_ij* that is neight contained in T* nor crossed by T, follow both of its endpoints back to
//    // the root of T*.
//    for (Edge e : mesh->edges()) {
//        Halfedge he = e.halfedge();
//        // Only consider edges that are not in the primal spanning tree or dual spanning cotree.
//        bool inT = inPrimalSpanningTree(he);
//        bool inTstar = inDualSpanningCotree(he);
//        if (!inT && !inTstar) {
//            // Find the two faces f and g that share the edge e.
//            // If the edge is not shared by two faces, skip it.
//            Face f = he.face();
//            Face g = he.twin().face();
//
//            std::vector<Halfedge> generator;
//            // Start with the halfedge he.
//            generator.push_back(he);
//            // Follow the face parent chain for face f until reaching the root face.
//            Face fparent = faceParent[f];
//            while (fparent != f) {
//                Halfedge sharedHe = sharedHalfedge(fparent, f); // Get the halfedge shared with the parent face
//                generator.push_back(sharedHe);
//                f = fparent;             // Move to the parent face
//                fparent = faceParent[f]; // Update the parent face
//            }
//
//            // Now follow the face parent chain for face g until reaching the root face.
//            std::vector<Halfedge> backward;
//            Face gparent = faceParent[g];
//            while (gparent != g) {
//                Halfedge sharedHe = sharedHalfedge(g, gparent); // Get the halfedge shared with the parent face
//                auto it = std::find(generator.begin(), generator.end(), sharedHe);
//                if (it != generator.end()) {
//                    generator.erase(it, generator.end()); // Remove the halfedge if it exists in the generator
//                    break;                                // Stop if we found a common halfedge
//                }
//                backward.push_back(sharedHe); // Add to backward chain
//                g = gparent;                  // Move to the parent face
//                gparent = faceParent[g];      // Update the parent face
//            }
//
//            // Append the backward chain in reverse order to the generator.
//            generator.insert(generator.end(), backward.rbegin(), backward.rend()); // HOW DOES THIS WORK?
//
//            // Add the generator to the list of generators.
//            generators.push_back(generator);
//        }
//    }
//    
//}
//
//void TreeCotree::buildGenerators() {
//
//    // For a closed surface of genus g, there are 2g generators.
//    int g = mesh->genus();
//    generators.clear();
//    generators.reserve(2 * g);
//
//    // Build spanning tree T* of dual edges.
//    buildDualSpanningCoTree();
//
//    // Build spanning tree T of primal edges.
//    buildPrimalSpanningTree();
//
//    /*
//    
//    */
//
//    // For each dual edge e_ij* that is neight contained in T* nor crossed by T, follow both of its endpoints back to the root of T*.
//    for (Edge e : mesh->edges()) {
//        Halfedge he = e.halfedge();
//        // Only consider edges that are not in the primal spanning tree or dual spanning cotree.
//        bool inT = inPrimalSpanningTree(he);
//        bool inTstar = inDualSpanningCotree(he);
//        if (!inT && !inTstar) {
//            // Find the two faces f and g that share the edge e.
//            // If the edge is not shared by two faces, skip it.
//            Face f = he.face();
//            Face g = he.twin().face();
//
//            std::vector<Halfedge> generator;
//            // Start with the halfedge he.
//            generator.push_back(he);
//            // Follow the face parent chain for face f until reaching the root face.
//            Face fparent = faceParent[f];
//            while (fparent != f) {
//                Halfedge sharedHe = sharedHalfedge(fparent, f); // Get the halfedge shared with the parent face
//                generator.push_back(sharedHe);
//                f = fparent; // Move to the parent face
//                fparent = faceParent[f]; // Update the parent face
//            }
//
//            // Now follow the face parent chain for face g until reaching the root face.
//            std::vector<Halfedge> backward;
//            Face gparent = faceParent[g];
//            while (gparent != g) {
//                Halfedge sharedHe = sharedHalfedge(g, gparent); // Get the halfedge shared with the parent face
//                auto it = std::find(generator.begin(), generator.end(), sharedHe);
//                if (it != generator.end()) {
//                    generator.erase(it, generator.end()); // Remove the halfedge if it exists in the generator
//                    break; // Stop if we found a common halfedge
//                }
//                backward.push_back(sharedHe); // Add to backward chain
//                g = gparent; // Move to the parent face
//                gparent = faceParent[g]; // Update the parent face
//            }
//
//            // Append the backward chain in reverse order to the generator.
//            generator.insert(generator.end(), backward.rbegin(), backward.rend()); // HOW DOES THIS WORK?
//
//            // Add the generator to the list of generators.
//            generators.push_back(generator);
//
//
//
//
//
//        }
//        if (inPrimalSpanningTree(he)) {
//            continue;
//        }
//        if (inDualSpanningCotree(he)) {
//            continue;
//        }
//        // TODO
//        // Find the two faces f and g that share the edge e.
//        Face f = he.face();
//        Face g = he.twin().face();
//        
//        generators.push_back(std::vector<Halfedge>());
//        std::vector<Halfedge>& generator = generators.back();
//        generator.push_back(he);
//        Face fparent = faceParent.at(f);
//        while (fparent != f) {
//            generator.push_back(sharedHalfedge(fparent, f));
//            f = fparent;
//            fparent = faceParent.at(f);
//        }
//
//        std::vector<Halfedge> backward;
//        Face gparent = faceParent.at(g);
//        while (gparent != g) {
//            he = sharedHalfedge(gparent, g);
//            auto it = std::find(generator.begin(), generator.end(), he);
//            if (it != generator.end()) {
//                generator.erase(it, generator.end());
//                break;
//            }
//            backward.push_back(sharedHalfedge(g, gparent));
//            g = gparent;
//            gparent = faceParent.at(g);
//        }
//
//        generator.insert(generator.end(), backward.rbegin(), backward.rend());
//    
//    }
//    // TODO: Build generators and populate this->generators
//}