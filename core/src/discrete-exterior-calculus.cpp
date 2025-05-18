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

typedef Eigen::Triplet<double> Triplet;

/*
 * Build Hodge operator on 0-forms.
 * By convention, the area of a vertex is 1.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar0Form() const {

     // Since the area of a vertex is 1, the diagonal elements are equal to the dual area of the vertex.
    size_t nVertices = mesh.nVertices();
    std::vector<Triplet> triplets(nVertices);
    for (size_t i = 0; i < nVertices; ++i) {
        Vertex v = mesh.vertex(i);
        double dualArea = barycentricDualArea(v);
        triplets[i] = Triplet(i, i, dualArea);
    }

    SparseMatrix<double> hodgeStar0Form(nVertices, nVertices);
    hodgeStar0Form.setFromTriplets(triplets.begin(), triplets.end());
    return hodgeStar0Form;
}

/*
 * Build Hodge operator on 1-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar1Form() const {

    std::vector<Triplet> triplets;
    for(Edge e : mesh.edges()) {
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
        triplets.push_back(Triplet(e.getIndex(), e.getIndex(), weight));
    }

    SparseMatrix<double> hodgeStar1Form(mesh.nEdges(), mesh.nEdges());
    hodgeStar1Form.setFromTriplets(triplets.begin(), triplets.end());
    return hodgeStar1Form;
}

/*
 * Build Hodge operator on 2-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 2-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar2Form() const {

    size_t nFaces = mesh.nFaces();
    std::vector<Triplet> triplets(nFaces);
    for (size_t i = 0; i < nFaces; ++i) {
        Face f = mesh.face(i);
        double area = faceArea(f);
        triplets[i] = Triplet(i, i, 1.0 / area);
    }

    SparseMatrix<double> hodgeStar2Form(nFaces, nFaces);
    hodgeStar2Form.setFromTriplets(triplets.begin(), triplets.end());
    return hodgeStar2Form;
}

/*
 * Build exterior derivative on 0-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative0Form() const {

    size_t nEdges = mesh.nEdges();
    size_t nVertices = mesh.nVertices();
    std::vector<Triplet> triplets;
    for(size_t i = 0; i < nEdges; ++i) {
        Halfedge he = mesh.edge(i).halfedge();
        Vertex head = he.tipVertex();
        Vertex tail = he.tailVertex();
        triplets.push_back(Triplet(i, head.getIndex(), 1.0));
        triplets.push_back(Triplet(i, tail.getIndex(), -1.0));
    }

    SparseMatrix<double> exteriorDerivative0Form(nEdges, nVertices);
    exteriorDerivative0Form.setFromTriplets(triplets.begin(), triplets.end());
    return exteriorDerivative0Form;
}

/*
 * Build exterior derivative on 1-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative1Form() const {

    size_t nFaces = mesh.nFaces();
    size_t nEdges = mesh.nEdges();
    std::vector<Triplet> triplets;
    for(size_t i = 0; i < nFaces; ++i) {
        Face f = mesh.face(i);
        for(Halfedge he : f.adjacentHalfedges()) {
            double orientation = he.getIndex() % 2 == 0 ? 1.0 : -1.0;
            triplets.push_back(Triplet(i, he.edge().getIndex(), orientation));
        }
    }

    SparseMatrix<double> exteriorDerivative1Form(nFaces, nEdges);
    exteriorDerivative1Form.setFromTriplets(triplets.begin(), triplets.end());
    return exteriorDerivative1Form;
}

} // namespace surface
} // namespace geometrycentral