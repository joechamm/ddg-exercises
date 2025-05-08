// PLEASE READ:
//
// This file implements additional geometry routines for the VertexPositionGeometry class in Geometry Central. Because
// we are "inside" the class, we no longer have to call
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
// Functions in this file can be called from other projects simply by using geometry->cotan(he),
// geometry->barycentricDualArea(v), etc. where "geometry" is a pointer to a VertexPositionGeometry. This avoids having
// to declare a GeometryRoutines object in every project, and also mimics the way that geometry routines are normally
// called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.


#include "geometrycentral/surface/vertex_position_geometry.h"
#include <complex>

namespace geometrycentral {
namespace surface {

/*
 * Compute the Euler characteristic of the mesh.
 */
int VertexPositionGeometry::eulerCharacteristic() const {
    return (int)mesh.nVertices() - (int)mesh.nEdges() + (int)mesh.nFaces();
}

/*
 * Compute the mean length of all the edges in the mesh.
 *
 * Input:
 * Returns: The mean edge length.
 */
double VertexPositionGeometry::meanEdgeLength() const {

    double total = 0.0;
    for (Edge e : mesh.edges()) {
        total += edgeLength(e);
    }
    return total / mesh.nEdges();
}

/*
 * Compute the total surface area of the mesh.
 *
 * Input:
 * Returns: The surface area of the mesh.
 */
double VertexPositionGeometry::totalArea() const {

    double total = 0.0;
    for (Face f : mesh.faces()) {
        total += faceArea(f);
    }
    return total;
}

/*
 * Computes the cotangent of the angle opposite to a halfedge. (Do NOT use built-in function for this)
 *
 * Input: The halfedge whose cotan weight is to be computed.
 * Returns: The cotan of the angle opposite the given halfedge.
 */
//double VertexPositionGeometry::cotan(Halfedge he) const {
//
//    // TODO
//    // 
//    // 
//    // 
//    //Vector3 a = inputVertexPositions[he.next().tipVertex()];
//    //Vector3 b = inputVertexPositions[he.tailVertex()];
//    //Vector3 c = inputVertexPositions[he.tipVertex()];
//    //Vector3 u = b - a;
//    //Vector3 v = c - a;
//    ////Vector3 u = inputVertexPositions[he.next().tipVertex()] -
//    ////            inputVertexPositions[he.tipVertex()]; // vector from tip of this halfedge to tip of next halfedge
//    ////Vector3 v = inputVertexPositions[he.tailVertex()] -
//    ////            inputVertexPositions[he.tipVertex()]; // vector from tip of this halfedge to tail of this halfedge
//
//    //return dot(u, v) / cross(u, v).norm(); // cos(theta) / sin(theta) = cot(theta)
//
//    // Let the triangle ijk be oriented with the halfedge ki, so that the tail of this halfedge is vertex k and the tip
//    // is vertex i. Using Heron's formula, the area of the triangle ijk is A_ijk = Sqrt[s * (s - l_ij) * (s - l_jk) * (s - l_ki)], 
//    // where s = (l_ij + l_jk + l_ki) / 2. We also have A_ijk = 1/2 * l_ki * l_jk * sin(theta_k_ij), where theta_k_ij is the angle at
//    // vertex k, and l_ki and l_jk are the lengths of the edges ki and jk, respectively. Combining these two equations with the
//    // law of cosines, (l_ij)^2 = (l_jk)^2 + (l_ki)^2 - 2 * l_jk * l_ki * cos(theta_k_ij), we can calculate cotan( theta_k_ij ) =
//    // cos( theta_k_ij ) / sin( theta_k_ij). Solving for cos and sin here we get:
//    
//    //      cos(theta_k_ij) = ( (l_jk)^2 + (l_ki)^2 - (l_ij)^2 ) / (2 * l_jk * l_ki)
//    //      sin(theta_k_ij) = 2 * A_ijk / (l_ki * l_jk)
//    //      cot(theta_k_ij) = cos(theta_k_ij) / sin(theta_k_ij) = ( (l_jk)^2 + (l_ki)^2 - (l_ij)^2 ) * ( l_ki * l_jk ) / ( (2 * l_jk * l_ki) * (2 * A_ijk) )
//    //      cot(theta_k_ij) = ( (l_jk)^2 + (l_ki)^2 - (l_ij)^2 ) / (4 * A_ijk)
//
//    // Let's compute the lengths of the edges l_ij, l_jk, and l_ki. We can use the halfedge to ensure we have the
//    // correct orientation. Since we're at k, the halfedge is ki, and we can traverse around to get ij then jk.
//    Halfedge he_ij = he.next(); // halfedge ij
//    double l_ki = norm(inputVertexPositions[he.tipVertex()] - inputVertexPositions[he.tailVertex()]); // length of edge ki
//    double l_ij = norm(inputVertexPositions[he_ij.tipVertex()] - inputVertexPositions[he.tipVertex()]); // length of edge ij
//    double l_jk = norm(inputVertexPositions[he.tailVertex()] - inputVertexPositions[he_ij.tipVertex()]); // length of edge jk
//    
//    // Use Heron's formula to compute the area of the triangle ijk.
//    double s = (l_ij + l_jk + l_ki) / 2.0; // semi-perimeter
//    double A_ijk = sqrt(s * (s - l_ij) * (s - l_jk) * (s - l_ki)); // area of triangle ijk
//    // Now we can compute the cotangent of the angle opposite the halfedge.
//    double cotan =
//        ((l_jk * l_jk) + (l_ki * l_ki) - (l_ij * l_ij)) / (4.0 * A_ijk); // cotangent of the angle opposite the halfedge
//    return cotan;
//}
double VertexPositionGeometry::cotan(Halfedge he) const {
    // TODO
    // The cotangent of the angle opposite the halfedge e_ij in triangle ijk is given by dot(e_kj, e_ki) / cross(e_kj, e_ki).norm().
    // We can get the halfedge e_ij from the halfedge he, and then we can get the other two edges using the halfedge's next and twin. 
    // The halfedge e_ij is the one that points from vertex i to vertex j, so we can get the other two edges by traversing the halfedge 
    // to get the other two vertices.
    Halfedge he_jk = he.next(); // halfedge jk
    Vector3 e_kj = inputVertexPositions[he_jk.tailVertex()] - inputVertexPositions[he_jk.tipVertex()]; // vector from k to j
    Vector3 e_ki = inputVertexPositions[he.tailVertex()] - inputVertexPositions[he_jk.tipVertex()]; // vector from k to i
    double dot_product = dot(e_kj, e_ki); // dot product of e_kj and e_ki = |e_kj| * |e_ki| * cos(theta_k_ij)
    double cross_productNorm = cross(e_kj, e_ki).norm(); // norm of cross product of e_kj and e_ki = |e_kj| * |e_ki| * sin(theta_k_ij)
    return (dot_product / cross_productNorm); // cotangent of the angle opposite the halfedge
}

/*
 * Computes the barycentric dual area of a vertex.
 *
 * Input: The vertex whose barycentric dual area is to be computed.
 * Returns: The barycentric dual area of the given vertex.
 */
double VertexPositionGeometry::barycentricDualArea(Vertex v) const {

    // TODO
    double total = 0.0;
    for (Face f : v.adjacentFaces()) {
        total += faceArea(f);
    }
    return (total / 3.0); 
}

/*
 * Computes the angle (in radians) at a given corner. (Do NOT use built-in function for this)
 *
 *
 * Input: The corner at which the angle needs to be computed.
 * Returns: The angle clamped between 0 and π.
 */
double VertexPositionGeometry::angle(Corner c) const {

    // TODO
    // Use law of cosines to get interior angle at corner i of triangle ijk. For vertices i, j, k, we first compute
    // l_ij, l_jk, and l_ki. We can use the halfedge to ensure we have the correct orientation. The interior angle theta_i_jk,
    // for the corner i in the ijk oriented triangle, will then be given by theta_i_jk = arccos( ( (l_ij)^2 + (l_ki)^2 - (l_jk)^2 ) / 2 * l_ij * l_ki ) ).
    // We're starting at vertex i, so the first halfedge is ij.
    //double l_ij = norm(inputVertexPositions[c.halfedge().tipVertex()] - inputVertexPositions[c.vertex()]); // length of edge ij
    //double l_jk = norm(inputVertexPositions[c.halfedge().next().tipVertex()] - inputVertexPositions[c.halfedge().tipVertex()]); // length of edge jk    
    //double l_ki = norm(inputVertexPositions[c.halfedge().next().next().tipVertex()] - inputVertexPositions[c.vertex()]); // length of edge ki
    //double ctheta_i_jk = ((l_ij * l_ij) + (l_ki * l_ki) - (l_jk * l_jk)) / (2.0 * l_ij * l_ki);
    //ctheta_i_jk = clamp(ctheta_i_jk, -1.0, 1.0); // clamp to avoid NaN
    //return acos(ctheta_i_jk);                    // angle at corner i
    Halfedge heA = c.halfedge();
    Halfedge heOpp = heA.next();
    Halfedge heB = heOpp.next();

    double lOpp = norm(inputVertexPositions[heOpp.tipVertex()] -
                       inputVertexPositions[heOpp.tailVertex()]); // length of edge opposite the corner
    double lA = norm(inputVertexPositions[heA.tipVertex()] - inputVertexPositions[heA.tailVertex()]);
    double lB = norm(inputVertexPositions[heB.tipVertex()] - inputVertexPositions[heB.tailVertex()]);

    double q = (lA * lA + lB * lB - lOpp * lOpp) / (2.0 * lA * lB);
    q = clamp(q, -1.0, 1.0);

    return std::acos(q);
}

/*
 * Computes the signed angle (in radians) between two adjacent faces. (Do NOT use built-in function for this)
 *
 * Input: The halfedge (shared by the two adjacent faces) on which the dihedral angle is computed.
 * Returns: The dihedral angle.
 */
double VertexPositionGeometry::dihedralAngle(Halfedge he) const {

    // TODO
    // The dihedral angle between the two normals, N_ijk and N_ijl, of two adjacent faces, ijk and ijl,
    // is given by atan2((e_ij / |e_ij|) . (N_ijk x N_jil), (N_ijk . N_jil))
    // where e_ij is the edge between vertices i and j, and N_ijk and N_ijl are the normals of the two faces.
    Face f1 = he.face();
    Face f2 = he.twin().face();
    Vector3 n1 = faceNormal(f1);
    Vector3 n2 = faceNormal(f2);
    Vector3 unitEdge = unit(halfedgeVector(he));
    Vector3 f1xf2 = cross(n1, n2);
    return atan2(dot(unitEdge, f1xf2), dot(n1, n2));
}

/*
 * Computes the normal at a vertex using the "equally weighted" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "equally weighted" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalEquallyWeighted(Vertex v) const {

    // TODO
    // we'll just add up the normals of all incident faces, then normalize
    Vector3 normalSum = {0.0, 0.0, 0.0};
    // for (Face f : v.adjacentFaces()) {
    //     // for each face, we need the normal of the face.
    //     // we can get the normal of the face using the halfedge of the face.
    //     Halfedge he = f.halfedge();
    //     Vector3 pA = inputVertexPositions[he.vertex()];
    //     he = he.next();
    //     Vector3 pB = inputVertexPositions[he.vertex()];
    //     he = he.next();
    //     Vector3 pC = inputVertexPositions[he.vertex()];
    //     normalSum += unit(cross(pB - pA, pC - pA));
    // }

    // return unit(normalSum);
    for (Face f : v.adjacentFaces()) {
        Halfedge fh1 = f.halfedge();
        Halfedge fh2 = fh1.next();
        Vector3 N = cross(inputVertexPositions[fh1.tipVertex()] - inputVertexPositions[fh1.tailVertex()],
                          inputVertexPositions[fh2.tipVertex()] - inputVertexPositions[fh2.tailVertex()]);
        normalSum += unit(N);
    }
    return unit(normalSum); // normalize the normal vector
}

/*
 * Computes the normal at a vertex using the "tip angle weights" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "tip angle weights" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalAngleWeighted(Vertex v) const {

    // TODO
    Vector3 normalSum = {0.0, 0.0, 0.0};
    for (Face f : v.adjacentFaces()) {
        // for each face, we need the corner angle at the vertex, and the normal of the face.
        // we can get the normal of the face using the halfedge of the face.
        // we can get the corner angle at the vertex using the halfedge of the face.
        Vertex vj, vk;
        Corner c;
        Halfedge he = f.halfedge();
        // For triangle face ijk, we'll assume that v is vertex i. There are three possibilities
        // for the halfedge.
        // 1. v is the tail vertex of the halfedge, in which case j will be the tip vertex and k will be the tip of
        // next. Our corner c will be the halfedge corner.
        // 2. v is the tip vertex of the halfedge, in which case k will be the tail vertex and j will be the tail of
        // next. Our corner will be the next halfedge's corner.
        // 3. v is neither tip nor tail of the halfedge, which means that the tail vertex is j and the tip vertex is k.
        // In this case, the corner will belong to the second halfedge.
        if (he.vertex() == v) {
            vj = he.tipVertex();
            vk = he.next().tipVertex();
            c = he.corner();
        } else if (he.tipVertex() == v) {
            vk = he.tailVertex();
            vj = he.next().tipVertex();
            c = he.next().corner();
        } else {
            vk = he.tipVertex();
            vj = he.tailVertex();
            c = he.next().next().corner();
        }

        Vector3 eij = inputVertexPositions[vj] - inputVertexPositions[v];
        Vector3 eik = inputVertexPositions[vk] - inputVertexPositions[v];
        Vector3 ejk = inputVertexPositions[vk] - inputVertexPositions[vj];
        double l_ij = norm(eij); // length of edge ij
        double l_ik = norm(eik); // length of edge ik
        double l_jk = norm(ejk); // length of edge jk


         // Use law of cosines to get interior angle at corner i of triangle ijk. For vertices i, j, k, we first compute
        // l_ij, l_ik, and l_jk. We can use the halfedge to ensure we have the correct orientation. The interior angle
        // theta_i_jk, for the corner i in the ijk oriented triangle, will then be given by theta_i_jk = arccos( (
        // (l_ij)^2 + (l_ki)^2 - (l_jk)^2 ) / 2 * l_ij * l_ki ) ). We're starting at vertex v = vi, so the first halfedge is
        // ij.

        double chteta_i_jk = clamp(((l_ij * l_ij) + (l_ik * l_ik) - (l_jk * l_jk)) / (2.0 * l_ij * l_ik), -1.0, 1.0); // clamp to avoid NaN
        double theta_i_jk = acos(chteta_i_jk); // angle at corner i
        // Now we can compute the normal of the face using the halfedge.
        Vector3 faceNormal = unit(cross(eij, eik));
        normalSum += theta_i_jk * faceNormal; // add the normal of the face weighted by the angle      
    }
    
    return unit(normalSum); // normalize the normal vector
}

/*
 * Computes the normal at a vertex using the "inscribed sphere" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "inscribed sphere" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalSphereInscribed(Vertex v) const {

    // TODO
    Vector3 normalSum = {0.0, 0.0, 0.0};
    for (Face f : v.adjacentFaces()) {
        // for each face, we need the corner angle at the vertex, and the normal of the face.
        // we can get the normal of the face using the halfedge of the face.
        // we can get the corner angle at the vertex using the halfedge of the face.
        Vertex vj, vk;
        Halfedge he = f.halfedge();
        if (he.vertex() == v) {
            vj = he.tipVertex();
            vk = he.next().tipVertex();
        } else if (he.tipVertex() == v) {
            vk = he.tailVertex();
            vj = he.next().tipVertex();
        } else {
            vk = he.tipVertex();
            vj = he.tailVertex();
        }

        Vector3 eij = inputVertexPositions[vj] - inputVertexPositions[v];
        Vector3 eik = inputVertexPositions[vk] - inputVertexPositions[v];
        double l_ijSq = dot(eij, eij); // squared length of eij
        double l_ikSq = dot(eik, eik); // squared length of eik
        double lsqInv = 1.0 / (l_ijSq * l_ikSq); // inverse of the squared length of the edge
        normalSum += cross(eij, eik) * lsqInv;   // add the normal of the face weighted by the angle       
    }

    return unit(normalSum); // normalize the normal vector
}

/*
 * Computes the normal at a vertex using the "face area weights" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "face area weighted" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalAreaWeighted(Vertex v) const {

    // TODO
    Vector3 normalSum = {0.0, 0.0, 0.0};
    for (Face f : v.adjacentFaces()) {
        // for each face, we need the corner angle at the vertex, and the normal of the face.
        // we can get the normal of the face using the halfedge of the face.
        // we can get the corner angle at the vertex using the halfedge of the face.
        Vertex vj, vk;
        Halfedge he = f.halfedge();
        if (he.vertex() == v) {
            vj = he.tipVertex();
            vk = he.next().tipVertex();
        } else if (he.tipVertex() == v) {
            vk = he.tailVertex();
            vj = he.next().tipVertex();
        } else {
            vk = he.tipVertex();
            vj = he.tailVertex();
        }
        Vector3 eij = inputVertexPositions[vj] - inputVertexPositions[v];
        Vector3 eik = inputVertexPositions[vk] - inputVertexPositions[v];
        Vector3 ejk = inputVertexPositions[vk] - inputVertexPositions[vj];
        // Use Heron's formula to find the face area.
        double l_ij = norm(eij);
        double l_ik = norm(eik);
        double l_jk = norm(ejk);
        double s = (l_ij + l_jk + l_ik) / 2.0;
        double A_ijk = sqrt(s * (s - l_ij) * (s - l_jk) * (s - l_ik)); // area of triangle ijki) * (s - l_jk) * (s - l_ik));
        Vector3 faceNormal = cross(eij, eik);
        normalSum += A_ijk * unit(faceNormal);
    }

    return unit(normalSum); // normalize the normal vector
}

/*
 * Computes the normal at a vertex using the "Gauss curvature" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "Gauss curvature" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalGaussianCurvature(Vertex v) const {

    // TODO
    Vector3 normalSum = {0.0, 0.0, 0.0};
    // for (Corner c : v.adjacentCorners()) {
    //     Halfedge he = c.halfedge();
    //     double dangle = dihedralAngle(he);
    //     Vector3 eij = (he.vertex() == v) ? inputVertexPositions[he.tipVertex()] - inputVertexPositions[v] :
    //     inputVertexPositions[v] - inputVertexPositions[he.vertex()]; normalSum += dangle * unit(eij);
    // }

    // return unit(normalSum); // normalize the normal vector

    for (Corner c : v.adjacentCorners()) {
        Halfedge he = c.halfedge();
        normalSum +=
            dihedralAngle(he) * unit(inputVertexPositions[he.tipVertex()] - inputVertexPositions[he.tailVertex()]);
    }
    return unit(normalSum); // normalize the normal vector
}
/*
 * Computes the normal at a vertex using the "mean curvature" method (equivalent to the "area gradient" method).
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "mean curvature" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalMeanCurvature(Vertex v) const {

    // TODO
    Vector3 normalSum = {0.0, 0.0, 0.0};
    //for (Edge e : v.adjacentEdges()) {
    //    Halfedge he_k = e.halfedge();
    //    double cotan_k = cotan(he_k);
    //    double cotan_l = cotan(he_k.twin());
    //    Vector3 eij = he_k.vertex() == v ? halfedgeVector(he_k) : halfedgeVector(he_k.twin());
    //    normalSum += (0.5 * (cotan_k + cotan_l)) * eij;
    //}

    //return unit(normalSum); // normalize the normal vector

    //for (Halfedge he : v.outgoingHalfedges()) {
    //    Vector3 e_ij = inputVertexPositions[he.tipVertex()] - inputVertexPositions[he.tailVertex()];
    //    normalSum += 0.5 * (cotan(he) + cotan(he.twin())) * e_ij;
    //
    //}
    //for (Halfedge he : v.outgoingHalfedges()) {
    //    Vector3 e_ij = inputVertexPositions[he.tipVertex()] - inputVertexPositions[v];
    //    normalSum += 0.5 * (cotan(he) + cotan(he.twin())) * e_ij;
    //}
    // The mean curvature normal is given by (1/2)(sum over all halfedges e_ij of (cotan(alpha_k_ij) + cotan(beta_l_ij) * e_ij), where alpha_k_ij is the angle opposite the halfedge e_ij in the triangle ijk, and beta_l_ij is the
    // angle opposite halfedge e_ij in triangle ilj.
    for (Halfedge he : v.outgoingHalfedges()) {
        Vector3 e_ij = inputVertexPositions[he.tipVertex()] - inputVertexPositions[v];
        double cotan_k = cotan(he);
        double cotan_l = cotan(he.twin());
        normalSum += (0.5 * (cotan_k + cotan_l)) * e_ij; // add the normal of the face weighted by the angle
    }

    return unit(normalSum); // normalize the normal vector
}

/*
 * Computes the angle defect at a vertex.
 *
 * Input: The vertex whose angle defect is to be computed.
 * Returns: The angle defect of the given vertex.
 */
double VertexPositionGeometry::angleDefect(Vertex v) const {

    // TODO
    double total = 0.0;
    for (Corner c : v.adjacentCorners()) {
        total += angle(c);
    }

    return (2.0 * PI) - total; // angle defect = 2 * pi - sum of angles
}

/*
 * Computes the total angle defect of the mesh.
 *
 * Input:
 * Returns: The total angle defect
 */
double VertexPositionGeometry::totalAngleDefect() const {

    double total = 0.0;
    for (Vertex v : mesh.vertices()) {
        total += angleDefect(v);
    }
    return total; // total angle defect = sum of angle defects of all vertices
}
/*
 * Computes the (integrated) scalar mean curvature at a vertex.
 *
 * Input: The vertex whose mean curvature is to be computed.
 * Returns: The mean curvature at the given vertex.
 */
double VertexPositionGeometry::scalarMeanCurvature(Vertex v) const {
    double curvature = 0.0;
    for (Halfedge he : v.outgoingHalfedges()) {
        Vector3 e_ij = inputVertexPositions[he.tipVertex()] - inputVertexPositions[v];
        double dangle = dihedralAngle(he);
        curvature += dangle * norm(e_ij);
    }
    return 0.5 * curvature;
}

/*
 * Computes the circumcentric dual area of a vertex.
 *
 * Input: The vertex whose circumcentric dual area is to be computed.
 * Returns: The circumcentric dual area of the given vertex.
 */
double VertexPositionGeometry::circumcentricDualArea(Vertex v) const {
    //double area = 0.0;
    //for (Corner c : v.adjacentCorners()) {
    //    Halfedge he = c.halfedge();
    //    Halfedge he_next = he.next().next();
    //    double l_ij = norm(inputVertexPositions[he.tipVertex()] - inputVertexPositions[he.tailVertex()]);
    //    double l_ik = norm(inputVertexPositions[he_next.tipVertex()] - inputVertexPositions[he_next.tailVertex()]);
    //    area += ((cotan(he_next) * l_ik * l_ik) + (cotan(he) * l_ij * l_ij));
    //}

    //return 1.0 / (8.0 * area);

    double area = 0.0;
    // The circumcentric dual area at vertex i is given by the sum over all adjacent faces ijk of cot(alpha_j_ki) * |e_ik|^2 + cot(beta_k_ij) * |e_ij|^2, over 8.
    for (Halfedge he_ij : v.outgoingHalfedges()) {
        Halfedge he_ki = he_ij.next().next();
        double cotAlpha_j_ki = cotan(he_ki);
        double cotBeta_k_ij = cotan(he_ij);
        Vector3 e_ki = inputVertexPositions[v] - inputVertexPositions[he_ki.vertex()];
        Vector3 e_ij = inputVertexPositions[he_ij.tipVertex()] - inputVertexPositions[v];
        double l_kiSq = dot(e_ki, e_ki); // squared length of e_ki
        double l_ijSq = dot(e_ij, e_ij);
        area += (cotAlpha_j_ki * l_kiSq + cotBeta_k_ij * l_ijSq);
    }

    return (area / 8.0);
 /*   for (Face f : v.adjacentFaces()) {
        Halfedge he = f.halfedge();
        Halfedge he_next = he.next();
        Vector3 e_ij = inputVertexPositions[he.tipVertex()] - inputVertexPositions[he.tailVertex()];
        Vector3 e_ik = inputVertexPositions[he_next.tipVertex()] - inputVertexPositions[he_next.tailVertex()];
        area += (cotan(he) * norm(e_ik) * norm(e_ik)) + (cotan(he_next) * norm(e_ij) * norm(e_ij));
    }*/

}

/*
 * Computes the (pointwise) minimum and maximum principal curvature values at a vertex.
 *
 * Input: The vertex on which the principal curvatures need to be computed.
 * Returns: A std::pair containing the minimum and maximum principal curvature values at a vertex.
 */
std::pair<double, double> VertexPositionGeometry::principalCurvatures(Vertex v) const {

    // TODO
    double A = circumcentricDualArea(v);
    double H = scalarMeanCurvature(v) / A;
    double K = angleDefect(v) / A;
    double delta = sqrt(H * H - K);
    return std::make_pair(H - delta, H + delta);    
}


/*
 * Builds the sparse POSITIVE DEFINITE Laplace matrix. Do this by building the negative semidefinite Laplace matrix,
 * multiplying by -1, and shifting the diagonal elements by a small constant (e.g. 1e-8).
 *
 * Input:
 * Returns: Sparse positive definite Laplace matrix for the mesh.
 */
SparseMatrix<double> VertexPositionGeometry::laplaceMatrix() const {

    SparseMatrix<double> d0(buildExteriorDerivative0Form());
    SparseMatrix<double> star1(buildHodgeStar1Form());
    return d0.transpose() * star1 * d0 + identityMatrix<double>(mesh.nVertices()) * 1e-8;
}

/*
 * Builds the sparse diagonal mass matrix containing the barycentric dual area of each vertex.
 *
 * Input:
 * Returns: Sparse mass matrix for the mesh.
 */
SparseMatrix<double> VertexPositionGeometry::massMatrix() const {

    // TODO
    return buildHodgeStar0Form();
}

/*
 * Builds the sparse complex POSITIVE DEFINITE Laplace matrix. Do this by building the negative semidefinite Laplace
 * matrix, multiplying by -1, and shifting the diagonal elements by a small constant (e.g. 1e-8).
 *
 * Input:
 * Returns: Sparse complex positive definite Laplace matrix for the mesh.
 */
SparseMatrix<std::complex<double>> VertexPositionGeometry::complexLaplaceMatrix() const {
    SparseMatrix<double> laplace = laplaceMatrix();
    SparseMatrix<std::complex<double>> A(laplace.rows(), laplace.cols());

    std::vector<Eigen::Triplet<std::complex<double>>> triplets;

    for (int i = 0; i < laplace.outerSize(); ++i) {
        for (typename SparseMatrix<double>::InnerIterator it(laplace, i); it; ++it) {
        
        triplets.emplace_back(it.row(), it.col(), std::complex<double>(it.value(), 0));
        }
    }

    A.setFromTriplets(triplets.begin(), triplets.end());
    A.makeCompressed();
    return A;
}

/*
 * Compute the center of mass of a mesh.
 */
Vector3 VertexPositionGeometry::centerOfMass() const {

    // Compute center of mass.
    Vector3 center = {0.0, 0.0, 0.0};
    for (Vertex v : mesh.vertices()) {
        center += inputVertexPositions[v];
    }
    center /= mesh.nVertices();

    return center;
}

/*
 * Centers a mesh about the origin.
 * Also rescales the mesh to unit radius if <rescale> == true.
 */
void VertexPositionGeometry::normalize(const Vector3& origin, bool rescale) {

    // Compute center of mass.
    Vector3 center = centerOfMass();

    // Translate to origin [of original mesh].
    double radius = 0.0;
    for (Vertex v : mesh.vertices()) {
        inputVertexPositions[v] -= center;
        radius = std::max(radius, inputVertexPositions[v].norm());
    }

    // Rescale.
    if (rescale) {
        for (Vertex v : mesh.vertices()) {
            inputVertexPositions[v] /= radius;
        }
    }

    // Translate to origin [of original mesh].
    for (Vertex v : mesh.vertices()) {
        inputVertexPositions[v] += origin;
    }
}

} // namespace surface
} // namespace geometrycentral