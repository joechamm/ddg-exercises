// PLEASE READ:
//
// This file additional geometry routines for the VertexPositionGeometry class in Geometry Central. Because we are
// "inside" the class, we no longer have to call
//
//          geometry->inputVertexPositions[v], etc.
//
// We can just call
//
//          this->inputVertexPositions[v], etc.
//
// or simply
//
//          inputVertexPositions[v], etc.
//
// In addition, we no longer access the corresponding surface mesh via
//
//          mesh->vertices(), etc.
//
// but instead <mesh> is not a pointer anymore, so we use
//
//          mesh.vertices(), etc.
//
// Functions in this file can be called from other projects simply by using geometry->buildHodgeStar0Form(), etc. where
// "geometry" is a pointer to a VertexPositionGeometry. This avoids having to declare a GeometryRoutines object in every
// project, and also mimics the way that geometry routines are normally called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.

#include "geometrycentral/surface/vertex_position_geometry.h"

namespace geometrycentral {
namespace surface {


/*
 * Build Hodge operator on 0-forms.
 * By convention, the area of a vertex is 1.
 * 
 * The discrete Hodge Star 0-form operator takes discrete 0-forms (vertex values) to their dual 2-forms (face values).
 * It is an |V| x |V| diagonal matrix, where |V| is the number of vertices in the mesh and the diagonal elements are the
 * dual areas of the vertices.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar0Form() const {

    // TODO
    // Since the area of a vertex is 1, the diagonal elements are equal to the dual area of the vertex.
    size_t nVertices = mesh.nVertices();
    std::vector<Eigen::Triplet<double>> triplets(nVertices);
    for (size_t i = 0; i < nVertices; ++i) {
        Vertex v = mesh.vertex(i);
        double dualArea = barycentricDualArea(v);
        triplets[i] = Eigen::Triplet<double>(i, i, dualArea);
    }

    SparseMatrix<double> hodgeStar0Form(nVertices, nVertices);
    hodgeStar0Form.setFromTriplets(triplets.begin(), triplets.end());
    return hodgeStar0Form;
}

/*
 * Build Inverse Hodge operator on 0-forms.
 * 
 * Input:
 * Returns: A sparse diagonal matrix that is the inverse of the matrix representation H0 of the Hodge star operator on
 * 0-forms. The diagonal elements are 1 / dual area of the vertices. 
 */
SparseMatrix<double> VertexPositionGeometry::buildInverseHodgeStar0Form() const {
    // TODO
    // The inverse of the Hodge star operator on 0-forms is the same as the Hodge star operator on 2-forms.
    size_t nVertices = mesh.nVertices();
    std::vector<Eigen::Triplet<double>> triplets(nVertices);
    for (size_t i = 0; i < nVertices; ++i) {
        Vertex v = mesh.vertex(i);
        double dualArea = barycentricDualArea(v);
        triplets[i] = Eigen::Triplet<double>(i, i, 1.0 / dualArea);
    }
    SparseMatrix<double> inverseHodgeStar0Form(nVertices, nVertices);
    inverseHodgeStar0Form.setFromTriplets(triplets.begin(), triplets.end());
    return inverseHodgeStar0Form;
}

/*
 * Build Hodge operator on 1-forms.
 *
 * The discrete Hodge Star 1-form operator takes discrte 1-forms (edge values) to their dual 1-forms (edge values).
 * It is an |E| x |E| diagonal matrix, where |E| is the number of edges in the mesh and the diagonal elements are the
 * lengths of the dual edges.
 * 
 * 
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar1Form() const {

    // TODO
    // The Hodge star operator on 1-forms is the cotangent Laplace operator, (1/2)(cot a + cot b)
    // where a and b are the angles opposite the halfedge. We can use the cotan method on the halfedge
    // to get cot a, and the twin operator on our halfedge allows us to get cot b.
    /* size_t nEdges = mesh.nEdges();
     std::vector<Eigen::Triplet<double>> triplets(nEdges);
     for (size_t i = 0; i < nEdges; ++i) {
         Edge e = mesh.edge(i);
         double cotA = cotan(e.halfedge());
         double cotB = cotan(e.halfedge().twin());
         triplets[i] = Eigen::Triplet<double>(i, i, 0.5 * (cotA + cotB));
     }

     SparseMatrix<double> hodgeStar1Form(nEdges, nEdges);
     hodgeStar1Form.setFromTriplets(triplets.begin(), triplets.end());
     return hodgeStar1Form;*/
    std::vector<Eigen::Triplet<double>> triplets;
    for (Edge e : mesh.edges()) {
        Halfedge he = e.halfedge();
        double cotanHe = 0.0;
        double cotanTwin = 0.0;
        if (he.isInterior()) {
            cotanHe = cotan(he);
        }
        if (he.twin().isInterior()) {
            cotanTwin = cotan(he.twin());
        }

        double weight = 0.5 * (cotanHe + cotanTwin);


 //       double area = (cotan(e.halfedge()) + cotan(e.halfedge().twin())) / 2.0;
        triplets.push_back(Eigen::Triplet<double>(e.getIndex(), e.getIndex(), weight));
    }

    SparseMatrix<double> hodgeStar1Form(mesh.nEdges(), mesh.nEdges());
    hodgeStar1Form.setFromTriplets(triplets.begin(), triplets.end());
    return hodgeStar1Form;
}

/*
 * Build Hodge operator on 2-forms.
 *
 * The discrete Hodge Star 2-form operator takes discrete 2-forms (face values) to their dual 0-forms (vertex values).
 * It is an |F| x |F| diagonal matrix, where |F| is the number of faces in the mesh and the diagonal elements are the
 * dual areas of the faces.
 * 
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 2-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar2Form() const {

    // TODO
    // The Hodge star operator on 2-forms is 1 over the area of the face.
    size_t nFaces = mesh.nFaces();
    std::vector<Eigen::Triplet<double>> triplets(nFaces);
    for (size_t i = 0; i < nFaces; ++i) {
        Face f = mesh.face(i);
        double area = faceArea(f);
        triplets[i] = Eigen::Triplet<double>(i, i, 1.0 / area);
    }

    SparseMatrix<double> hodgeStar2Form(nFaces, nFaces);
    hodgeStar2Form.setFromTriplets(triplets.begin(), triplets.end());
    return hodgeStar2Form;
}

/*
 * Build exterior derivative on 0-forms.
 *
 * The discrete exterior derivative on 0-forms takes discrete 0-forms (vertex values) to their 1-forms (edge
 * values). It is an |E| x |V| sparse matrix, where |E| is the number of edges and |V| is the number of vertices.
 * For each row i, the entry in column j is 1 if the vertex j is the head of edge i, and -1 if it is the tail of edge i.
 * 
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative0Form() const {

    // TODO
    // The exterior derivative on 0-forms is an |E| x |V| matrix, where |E| is the number of edges and |V| is the number of
    // vertices. Each row corresponds to an edge, and each column corresponds to a vertex. The entries are 1 if the
    // vertex is the head of the edge, -1 if it is the tail of the edge, and 0 otherwise.
    size_t nEdges = mesh.nEdges();
    size_t nVertices = mesh.nVertices();
    std::vector<Eigen::Triplet<double>> triplets;
    for (size_t i = 0; i < nEdges; ++i) {
        Halfedge he = mesh.edge(i).halfedge();
        Vertex head = he.tipVertex();
        Vertex tail = he.tailVertex();
        triplets.push_back(Eigen::Triplet<double>(i, head.getIndex(), 1.0));
        triplets.push_back(Eigen::Triplet<double>(i, tail.getIndex(), -1.0));        
    }

    SparseMatrix<double> exteriorDerivative0Form(nEdges, nVertices);
    exteriorDerivative0Form.setFromTriplets(triplets.begin(), triplets.end());
    return exteriorDerivative0Form;
}

/*
 * Build exterior derivative on 1-forms.
 *
 * The exterior derivative on 1-forms takes discrete 1-forms (edge values) to 2-form (face values). It is an |F| x |E| sparse matrix, Row i corresponds to face i, and column j corresponds to edge j. The entries are 1 if the face is
 * oriented in the same direction as the edge, and -1 if it is oriented in the opposite direction.
 * 
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative1Form() const {

    // TODO
    // The exterior derivative on 1-forms is an m x n matrix, where m is the number of faces and n is the number of
    // edges. Each row corresponds to a face, and each column corresponds to an edge. The entries are 1 if the face is
    // oriented in the same direction as the edge, -1 if it is oriented in the opposite direction, and 0 otherwise.
    size_t nFaces = mesh.nFaces();
    size_t nEdges = mesh.nEdges();
    std::vector<Eigen::Triplet<double>> triplets;
    for (size_t i = 0; i < nFaces; ++i) {
        Face f = mesh.face(i);        
        for (Halfedge he : f.adjacentHalfedges()) {
            double orientation = he.getIndex() < he.twin().getIndex() ? 1.0 : -1.0;
            triplets.push_back(Eigen::Triplet<double>(i, he.edge().getIndex(), orientation));
        }
    }

    SparseMatrix<double> exteriorDerivative1Form(nFaces, nEdges);
    exteriorDerivative1Form.setFromTriplets(triplets.begin(), triplets.end());
    return exteriorDerivative1Form;
}

SparseMatrix<double> VertexPositionGeometry::buildVertexEdgeAdjacencyMatrix() const {
    size_t numVerts = mesh.nVertices();
    size_t numEdges = mesh.nEdges();
    std::vector<Eigen::Triplet<double>> triplets;
    for (Vertex v : mesh.vertices()) {
        size_t vIndex = v.getIndex();
        Eigen::Index col = (Eigen::Index)vIndex;
        for (Edge e : v.adjacentEdges()) {
            size_t eIndex = e.getIndex();
            Eigen::Index row = (Eigen::Index)eIndex;
            triplets.push_back(Eigen::Triplet<double>(row, col, 1.0));
        }
    }

    Eigen::Index rows = (Eigen::Index)numEdges;
    Eigen::Index cols = (Eigen::Index)numVerts;
    SparseMatrix<double> vertexEdgeAdjacency(rows, cols);
    vertexEdgeAdjacency.setFromTriplets(triplets.begin(), triplets.end());
    return vertexEdgeAdjacency;
}

} // namespace surface
} // namespace geometrycentral