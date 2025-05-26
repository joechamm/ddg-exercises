// Implement member functions for TreeCotree class.
#include "tree-cotree.h"
#include <queue>

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

    // TODO
    std::queue<Vertex> vertexQueue;
    Vertex root = *(mesh->vertices().begin());

    vertexQueue.push(root);
    vertexParent[root] = root;
    while (!vertexQueue.empty()) {
        Vertex v = vertexQueue.front();
        vertexQueue.pop();
        for (Halfedge he : v.outgoingHalfedges()) {
            Vertex w = he.tipVertex();
            if (!vertexParent.count(w)) {
                vertexParent[w] = v; // Set parent of w to v
                vertexQueue.push(w); // Add w to the queue
            }
        }
    }
}
/*
 * Check whether a halfedge is in the primal spanning tree.
 *
 * Input: A halfedge <he>
 * Returns: True if <he> is in the primal spanning tree, false otherwise.
 */
bool TreeCotree::inPrimalSpanningTree(Halfedge he) {
    Vertex v = he.vertex();
    Vertex w = he.tipVertex();

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
    
    // TODO
    return false; // placeholder
}

/*
 * Build a dual spanning tree on a mesh without boundary. More specificially, populate the member variable <faceParent>,
 * which is a std::map that maps each face of the input mesh to its parent in the dual spanning tree.
 *
 * Input:
 * Returns:
 */
void TreeCotree::buildDualSpanningCoTree() {

    // TODO
    std::queue<Face> faceQueue;
    Face root = *(mesh->faces().begin());
    faceQueue.push(root);
    faceParent[root] = root;

    while (!faceQueue.empty()) {
        Face f = faceQueue.front();
        faceQueue.pop();
        for (Halfedge he : f.adjacentHalfedges()) {
            if (inPrimalSpanningTree(he)) {
                continue;
            }

            Face g = he.twin().face();
            if (!faceParent.count(g)) {
                faceParent[g] = f; // Set parent of g to f
                faceQueue.push(g); // Add g to the queue
            }
        }
    
    }

}

/*
 * Check whether a halfedge is in the dual spanning tree.
 *
 * Input: A halfedge <he>
 * Returns: True if <he> is in the dual spanning tree, false otherwise.
 */
bool TreeCotree::inDualSpanningCotree(Halfedge he) {

    // TODO
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
    return false; // placeholder
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

    // order doesn't matter in a mesh without boundary
    buildPrimalSpanningTree();
    buildDualSpanningCoTree();

    for (Edge e : mesh->edges()) {
        Halfedge he = e.halfedge();
        if (inPrimalSpanningTree(he)) {
            continue;
        }
        if (inDualSpanningCotree(he)) {
            continue;
        }
        // TODO
        // Find the two faces f and g that share the edge e.
        Face f = he.face();
        Face g = he.twin().face();
        
        generators.push_back(std::vector<Halfedge>());
        std::vector<Halfedge>& generator = generators.back();
        generator.push_back(he);
        Face fparent = faceParent.at(f);
        while (fparent != f) {
            generator.push_back(sharedHalfedge(fparent, f));
            f = fparent;
            fparent = faceParent.at(f);
        }

        std::vector<Halfedge> backward;
        Face gparent = faceParent.at(g);
        while (gparent != g) {
            he = sharedHalfedge(gparent, g);
            auto it = std::find(generator.begin(), generator.end(), he);
            if (it != generator.end()) {
                generator.erase(it, generator.end());
                break;
            }
            backward.push_back(sharedHalfedge(g, gparent));
            g = gparent;
            gparent = faceParent.at(g);
        }

        generator.insert(generator.end(), backward.rbegin(), backward.rend());
    
    }
    // TODO: Build generators and populate this->generators
}