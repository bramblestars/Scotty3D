#include <float.h>
#include <assert.h>
#include <math.h>
#include "meshEdit.h"
#include "mutablePriorityQueue.h"
#include "error_dialog.h"

#include<windows.h>

namespace CMU462 {

bool isTriangle(FaceIter f) {
    return f->degree() == 3;
}

VertexIter HalfedgeMesh::splitEdge(EdgeIter e0) {
    // This method should split the given edge and return an iterator to the
    // newly inserted vertex. The halfedge of this vertex should point along
    // the edge that was split, rather than the new edges.

    HalfedgeIter h1 = e0->halfedge();
    HalfedgeIter h1_twin = h1->twin();

    // Ensures only h1 can be a boundary halfedge
    assert(!(h1->isBoundary() && h1_twin->isBoundary()));
    if (h1_twin->isBoundary()) {
        h1 = h1_twin;
        h1_twin = h1->twin();
    }
    
    bool on_boundary = h1->isBoundary();
    bool triangle1 = isTriangle(h1->face());
    bool triangle2 = isTriangle(h1_twin->face());

    if ((!triangle1 && !on_boundary) || !triangle2) {
        return h1->vertex();
    }

    /* GOAL:
                   . <- v2
                  /|\
                 / | \ 
                /  |  \
               / e0|   \
              /    |    \
        v3 ->'-----v-----' <- v4
              \    |    /
               \ e1|   /
                \  |  /
                 \ | /
                  \|/
                   ' <-v1
    */

    //set faces' halfedges to h1 and h2 just in case
    h1->face()->halfedge() = h1;
    h1_twin->face()->halfedge() = h1_twin;

    VertexIter v1 = h1->vertex();
    VertexIter v2 = h1_twin->vertex();
    VertexIter v3 = h1->next()->twin()->vertex();
    VertexIter v4 = h1_twin->next()->twin()->vertex();

    HalfedgeIter v1_in = h1_twin;
    do {
        v1_in = v1_in->next()->twin();
    } while (v1_in->next() != h1);

    //new vertex added to midpoint
    VertexIter v = newVertex();
    v->position = e0->centroid();
    v->halfedge() = h1;

    FaceIter f1, f2;

    //create 2 new faces
    if (!on_boundary) {
        f1 = newFace();
    }
    else {
        f1 = h1->face();
    }
    f2 = newFace();

    //this new edge goes along half of the original e0
    HalfedgeIter h2 = newHalfedge();
    HalfedgeIter h2_twin = newHalfedge();
    EdgeIter e1 = newEdge();
    e1->halfedge() = h2;

    h2->setNeighbors(h1, h2_twin, v1, e1, f1);
    h2_twin->setNeighbors(h1_twin->next(), h2, v, e1, f2);

    h1->next()->next()->face() = f1;
    f1->halfedge() = h2;

    h2_twin->next()->face() = f2;
    f2->halfedge() = h2_twin;

    h1->vertex() = v;
    if (!on_boundary) {
        h1->next()->next()->next() = h2;
        h1_twin->next() = h2_twin;
    }

    else {
        v1_in->next() = h2;
    }
    
    v->halfedge() = h1;
    v1->halfedge() = h2;
    v2->halfedge() = h1_twin;

    //this new edge splits the old face of h1
    if (!on_boundary) {
        HalfedgeIter split1 = newHalfedge();
        HalfedgeIter split1_twin = newHalfedge();
        EdgeIter f1_split = newEdge();
        f1_split->halfedge() = split1;

        split1->setNeighbors(h1, split1_twin, v3, f1_split, h1->face());
        split1_twin->setNeighbors(h1->next()->next(), split1, v, f1_split, f1);

        h2->next() = split1_twin;
        h1->next()->next() = split1;
    }
    
    //this new edge splits the old face of h1's twin
    HalfedgeIter split2 = newHalfedge();
    HalfedgeIter split2_twin = newHalfedge();
    EdgeIter f2_split = newEdge();
    f2_split->halfedge() = split2;
    f2->halfedge() = split2_twin;

    split2->setNeighbors(h2_twin->next()->next(), split2_twin, v,
        f2_split, h1_twin->face());
    split2_twin->setNeighbors(h2->twin(), split2, v4, f2_split, f2);

    h1_twin->next() = split2;
    h2_twin->next()->next() = split2_twin;

    checkConsistency();

    return v;
}

VertexIter HalfedgeMesh::collapseEdge(EdgeIter e) {
    // This method should collapse the given edge and return an iterator to
    // the new vertex created by the collapse.
    HalfedgeIter h = e->halfedge();
    HalfedgeIter h_twin = h->twin();

    VertexIter v1 = h->vertex();
    VertexIter v2 = h_twin->vertex(); // to be deleted

    if (v1->degree() == 1 || v2->degree() == 1) {
        return v1;
    }

    v1->position = e->centroid();

    HalfedgeIter v1_out = h_twin->next();
    HalfedgeIter v1_in = v1_out->twin();
    do {
        v1_in = v1_in->next()->twin();
    } while (v1_in->next() != h);

    HalfedgeIter v2_out = h->next();
    HalfedgeIter v2_in = v2_out->twin();
    do {
        v2_in = v2_in->next()->twin();
    } while (v2_in->next() != h_twin);

    bool triangle_h = isTriangle(h->face());
    bool triangle_htwin = isTriangle(h_twin->face());

    HalfedgeIter iter = v2->halfedge();

    do {
        iter->vertex() = v1;
        iter = iter->twin()->next();
    } while (iter != v2->halfedge());

    v1_in->next() = v2_out;
    v2_in->next() = v1_out;

    v1->halfedge() = v1_out;

    h->face()->halfedge() = v1_in;
    h_twin->face()->halfedge() = v2_in;

    deleteEdge(e);
    deleteHalfedge(h);
    deleteHalfedge(h_twin);

    // delete some stuff if we're collapsing edges incident to triangles
    if (triangle_h) {
        v1_in->next() = v2_out->twin()->next();
        v1_in->face() = v1_in->next()->face();
        v1_in->face()->halfedge() = v1_in;

        iter = v1_in->next();
        do {
            iter = iter->next();
        } while (iter->next() != v2_out->twin());

        iter->next() = v1_in;
        v1_in->vertex()->halfedge() = v1_in;
        v1->halfedge() = v1_in->twin();

        deleteFace(v2_out->face()); 
        deleteEdge(v2_out->edge());
        deleteHalfedge(v2_out->twin());
        deleteHalfedge(v2_out);
    }

    
    if (triangle_htwin) {
        v2_in->next() = v1_out->twin()->next();
        v2_in->face() = v2_in->next()->face();
        v2_in->face()->halfedge() = v2_in;

        iter = v2_in->next();
        do {
            iter = iter->next();
        } while (iter->next() != v1_out->twin());

        iter->next() = v2_in;
        v2_in->vertex()->halfedge() = v2_in;
        v1->halfedge() = v2_in->twin();

        deleteFace(v1_out->face());
        deleteEdge(v1_out->edge());
        deleteHalfedge(v1_out->twin());
        deleteHalfedge(v1_out);
    } 

    deleteVertex(v2);

    checkConsistency();

    return v1;
}

VertexIter HalfedgeMesh::collapseFace(FaceIter f) {
  // TODO: (meshEdit)
  // This method should collapse the given face and return an iterator to
  // the new vertex created by the collapse.
  showError("collapseFace() not implemented.");
  return VertexIter();
}

FaceIter HalfedgeMesh::eraseVertex(VertexIter v) {
  // TODO: (meshEdit)
  // This method should replace the given vertex and all its neighboring
  // edges and faces with a single face, returning the new face.

  return FaceIter();
}

FaceIter HalfedgeMesh::eraseEdge(EdgeIter e) {
    // This method should erase the given edge and return an iterator to the
    // merged face.

    if (e->isBoundary()) {
        return e->halfedge()->face();
    }

    HalfedgeIter h = e->halfedge();
    HalfedgeIter h_twin = h->twin();

    VertexIter v1 = h->vertex();
    VertexIter v2 = h_twin->vertex();

    FaceIter f1 = h->face();
    FaceIter f2 = h_twin->face(); // to be deleted if different from f1

    //if edge is incident to a single face (i.e. "sticks out" into the face)
    if (f1 == f2) {

        // in order to not allow edges with 2 vertices of degree 1,
        // do nothing if the selected edge is not the last edge in path
        if (h->next() != h_twin && h_twin->next() != h) {
            checkConsistency();
            return f1;
        }

        //swap so h->next() is always h_twin
        if (h->next() != h_twin) {
            h = h_twin;
            h_twin = h->twin();
        }

        HalfedgeIter e_out = h_twin->next();
        HalfedgeIter e_in = e_out;
        do {
            e_in = e_in->next();
        } while (e_in->next() != h);

        e_in->next() = e_out;
        f1->halfedge() = e_in;
        h->vertex()->halfedge() = e_out;
        
        deleteVertex(h_twin->vertex());
        deleteEdge(e);
        deleteHalfedge(h);
        deleteHalfedge(h_twin);

        checkConsistency();

        return f1;
    }

    HalfedgeIter v1_halfedge1 = h_twin->next(); //out of v1
    HalfedgeIter v1_halfedge2 = h; //into v1

    do {
        v1_halfedge2 = v1_halfedge2->next();
    } while (v1_halfedge2->next() != h);

    HalfedgeIter v2_halfedge1 = h->next(); //out of v2
    HalfedgeIter v2_halfedge2 = h_twin; //into v2

    do {
        v2_halfedge2 = v2_halfedge2->next();
        v2_halfedge2->face() = f1;
    } while (v2_halfedge2->next() != h_twin);

    v1_halfedge2->next() = v1_halfedge1;
    v2_halfedge2->next() = v2_halfedge1;

    //in case v1's halfedge is h or v2's halfedge is h_twin
    v1->halfedge() = v1_halfedge1;
    v2->halfedge() = v2_halfedge1;

    //in case f1's halfedge is h
    f1->halfedge() = v1_halfedge1;

    deleteFace(f2);
    deleteEdge(e);
    deleteHalfedge(h);
    deleteHalfedge(h_twin);
    
    checkConsistency();

    return f1;
}

EdgeIter HalfedgeMesh::flipEdge(EdgeIter e0) {
    // TODO: (meshEdit)
    // This method should flip the given edge and return an iterator to the
    // flipped edge.

    if (e0->isBoundary()) {
        return e0;
    }

    HalfedgeIter h = e0->halfedge();
    HalfedgeIter h_twin = h->twin();

    if (h->face() == h_twin->face()) {
        return e0;
    }

    // original endpoints of the edge
    VertexIter v1 = h->vertex();
    VertexIter v2 = h_twin->vertex();

    // edges that will be affected
    HalfedgeIter v2v3 = h->next();
    HalfedgeIter v1v4 = h_twin->next();
    HalfedgeIter v3_out = h->next()->next();
    HalfedgeIter v4_out = h_twin->next()->next();

    // new endpoints of the edge
    VertexIter v3 = h->next()->twin()->vertex();
    VertexIter v4 = h_twin->next()->twin()->vertex();

    eraseEdge(e0);

    EdgeIter e = newEdge();
    FaceIter f = newFace(); //touching h1, since eraseEdge erased one face
    HalfedgeIter h1 = newHalfedge(); //v4v3
    HalfedgeIter h1_twin = newHalfedge(); //v3v4

    e->halfedge() = h1;
    f->halfedge() = h1;
    h1->setNeighbors(v3_out, h1_twin, v4, e, f);
    h1_twin->setNeighbors(v4_out, h1, v3, e, v4_out->face());

    v4_out->face()->halfedge() = h1_twin;

    v1v4->next() = h1;
    v2v3->next() = h1_twin;

    //set halfedges to new face
    HalfedgeIter iter = h1;

    do {
        iter = iter->next();
        iter->face() = f;
    } while (iter != h1);

    checkConsistency();

    return e;
}

void HalfedgeMesh::subdivideQuad(bool useCatmullClark) {
  // Unlike the local mesh operations (like bevel or edge flip), we will perform
  // subdivision by splitting *all* faces into quads "simultaneously."  Rather
  // than operating directly on the halfedge data structure (which as you've
  // seen
  // is quite difficult to maintain!) we are going to do something a bit nicer:
  //
  //    1. Create a raw list of vertex positions and faces (rather than a full-
  //       blown halfedge mesh).
  //
  //    2. Build a new halfedge mesh from these lists, replacing the old one.
  //
  // Sometimes rebuilding a data structure from scratch is simpler (and even
  // more
  // efficient) than incrementally modifying the existing one.  These steps are
  // detailed below.

  // TODO Step I: Compute the vertex positions for the subdivided mesh.  Here
  // we're
  // going to do something a little bit strange: since we will have one vertex
  // in
  // the subdivided mesh for each vertex, edge, and face in the original mesh,
  // we
  // can nicely store the new vertex *positions* as attributes on vertices,
  // edges,
  // and faces of the original mesh.  These positions can then be conveniently
  // copied into the new, subdivided mesh.
  // [See subroutines for actual "TODO"s]
  if (useCatmullClark) {
    computeCatmullClarkPositions();
  } else {
    computeLinearSubdivisionPositions();
  }

  // TODO Step II: Assign a unique index (starting at 0) to each vertex, edge,
  // and
  // face in the original mesh.  These indices will be the indices of the
  // vertices
  // in the new (subdivided mesh).  They do not have to be assigned in any
  // particular
  // order, so long as no index is shared by more than one mesh element, and the
  // total number of indices is equal to V+E+F, i.e., the total number of
  // vertices
  // plus edges plus faces in the original mesh.  Basically we just need a
  // one-to-one
  // mapping between original mesh elements and subdivided mesh vertices.
  // [See subroutine for actual "TODO"s]
  assignSubdivisionIndices();

  // TODO Step III: Build a list of quads in the new (subdivided) mesh, as
  // tuples of
  // the element indices defined above.  In other words, each new quad should be
  // of
  // the form (i,j,k,l), where i,j,k and l are four of the indices stored on our
  // original mesh elements.  Note that it is essential to get the orientation
  // right
  // here: (i,j,k,l) is not the same as (l,k,j,i).  Indices of new faces should
  // circulate in the same direction as old faces (think about the right-hand
  // rule).
  // [See subroutines for actual "TODO"s]
  vector<vector<Index> > subDFaces;
  vector<Vector3D> subDVertices;
  buildSubdivisionFaceList(subDFaces);
  buildSubdivisionVertexList(subDVertices);

  // TODO Step IV: Pass the list of vertices and quads to a routine that clears
  // the
  // internal data for this halfedge mesh, and builds new halfedge data from
  // scratch,
  // using the two lists.
  rebuild(subDFaces, subDVertices);
}

/**
 * Compute new vertex positions for a mesh that splits each polygon
 * into quads (by inserting a vertex at the face midpoint and each
 * of the edge midpoints).  The new vertex positions will be stored
 * in the members Vertex::newPosition, Edge::newPosition, and
 * Face::newPosition.  The values of the positions are based on
 * simple linear interpolation, e.g., the edge midpoints and face
 * centroids.
 */
void HalfedgeMesh::computeLinearSubdivisionPositions() {
    // TODO For each vertex, assign Vertex::newPosition to
    // its original position, Vertex::position.

    // TODO For each edge, assign the midpoint of the two original
    // positions to Edge::newPosition.

    // TODO For each face, assign the centroid (i.e., arithmetic mean)
    // of the original vertex positions to Face::newPosition.  Note
    // that in general, NOT all faces will be triangles!
  
    for (VertexIter v = verticesBegin(); v != verticesEnd(); v++) {
        v->newPosition = v->position;
    }

    for (EdgeIter e = edgesBegin(); e != edgesEnd(); e++) {
        e->newPosition = e->centroid();
    }

    for (FaceIter f = facesBegin(); f != facesEnd(); f++) {
        f->newPosition = f->centroid();
    }

}

/**
 * Compute new vertex positions for a mesh that splits each polygon
 * into quads (by inserting a vertex at the face midpoint and each
 * of the edge midpoints).  The new vertex positions will be stored
 * in the members Vertex::newPosition, Edge::newPosition, and
 * Face::newPosition.  The values of the positions are based on
 * the Catmull-Clark rules for subdivision.
 */

Vector3D averageFacePositions(VertexIter v) {
    HalfedgeIter start = v->halfedge();
    HalfedgeIter curr = start;
    Vector3D sumPositions = 0;

    do {
        sumPositions += curr->face()->newPosition;
        curr = curr->twin()->next();
    } while (curr != start);

    return sumPositions / v->degree();
}

Vector3D averageEdgePositions(VertexIter v) {
    HalfedgeIter start = v->halfedge();
    HalfedgeIter curr = start;
    Vector3D sumPositions = 0;

    do {
        sumPositions += curr->edge()->centroid();
        curr = curr->twin()->next();
    } while (curr != start);

    return sumPositions / v->degree();
}

void HalfedgeMesh::computeCatmullClarkPositions() {
    // The implementation for this routine should be
    // a lot like HalfedgeMesh::computeLinearSubdivisionPositions(),
    // except that the calculation of the positions themsevles is
    // slightly more involved, using the Catmull-Clark subdivision
    // rules. (These rules are outlined in the Developer Manual.)

    for (FaceIter f = facesBegin(); f != facesEnd(); f++) {
        f->newPosition = f->centroid();
    }

    Vector3D f_pos1, f_pos2;
    Vector3D v_pos1, v_pos2;
    for (EdgeIter e = edgesBegin(); e != edgesEnd(); e++) {
        f_pos1 = e->halfedge()->face()->newPosition;
        f_pos2 = e->halfedge()->twin()->face()->newPosition;
        v_pos1 = e->halfedge()->vertex()->position;
        v_pos2 = e->halfedge()->twin()->vertex()->position;

        e->newPosition = (f_pos1 + f_pos2 + v_pos1 + v_pos2) / 4;
    }

    Vector3D faceAverage, edgeAverage, originalPosition;
    for (VertexIter v = verticesBegin(); v != verticesEnd(); v++) {
        faceAverage = averageFacePositions(v);
        edgeAverage = averageEdgePositions(v);
        originalPosition = v->position;
     
        v->newPosition = (faceAverage + 2 * edgeAverage +
            (v->degree() - 3) * originalPosition) / v->degree();
    }
    
}

/**
 * Assign a unique integer index to each vertex, edge, and face in
 * the mesh, starting at 0 and incrementing by 1 for each element.
 * These indices will be used as the vertex indices for a mesh
 * subdivided using Catmull-Clark (or linear) subdivision.
 */
void HalfedgeMesh::assignSubdivisionIndices() {
    // TODO Start a counter at zero; if you like, you can use the
    // "Index" type (defined in halfedgeMesh.h)

    Index counter = 0;

    // TODO Iterate over vertices, assigning values to Vertex::index
    for (VertexIter v = verticesBegin(); v != verticesEnd(); v++) {
        v->index = counter;
        counter++;
    }

    // TODO Iterate over edges, assigning values to Edge::index
    for (EdgeIter e = edgesBegin(); e != edgesEnd(); e++) {
        e->index = counter;
        counter++;
    }

    // TODO Iterate over faces, assigning values to Face::index
    for (FaceIter f = facesBegin(); f != facesEnd(); f++) {
        f->index = counter;
        counter++;
    }
}

/**
 * Build a flat list containing all the vertex positions for a
 * Catmull-Clark (or linear) subdivison of this mesh.  The order of
 * vertex positions in this list must be identical to the order
 * of indices assigned to Vertex::newPosition, Edge::newPosition,
 * and Face::newPosition.
 */
void HalfedgeMesh::buildSubdivisionVertexList(vector<Vector3D>& subDVertices) {
    // Resize the vertex list so that it can hold all the vertices.
    subDVertices.resize(this->vertices.size() + this->edges.size() 
        + this->faces.size());

    // Iterate over vertices, assigning Vertex::newPosition to the
    // appropriate location in the new vertex list.
    Index index = 0;
    for (VertexIter v = verticesBegin(); v != verticesEnd(); v++) {
        subDVertices[index] = v->newPosition;
        index++;
    }

    // Iterate over edges, assigning Edge::newPosition to the appropriate
    // location in the new vertex list.
    for (EdgeIter e = edgesBegin(); e != edgesEnd(); e++) {
        subDVertices[index] = e->newPosition;
        index++;
    }

    // Iterate over faces, assigning Face::newPosition to the appropriate
    // location in the new vertex list.
    for (FaceIter f = facesBegin(); f != facesEnd(); f++) {
        subDVertices[index] = f->newPosition;
        index++;
    }

}

/**
 * Build a flat list containing all the quads in a Catmull-Clark
 * (or linear) subdivision of this mesh.  Each quad is specified
 * by a vector of four indices (i,j,k,l), which come from the
 * members Vertex::index, Edge::index, and Face::index.  Note that
 * the ordering of these indices is important because it determines
 * the orientation of the new quads; it is also important to avoid
 * "bowties."  For instance, (l,k,j,i) has the opposite orientation
 * of (i,j,k,l), and if (i,j,k,l) is a proper quad, then (i,k,j,l)
 * will look like a bowtie.
 */
void HalfedgeMesh::buildSubdivisionFaceList(vector<vector<Index> >& subDFaces) {
    // TODO This routine is perhaps the most tricky step in the construction of
    // a subdivision mesh (second, perhaps, to computing the actual Catmull-Clark
    // vertex positions).  Basically what you want to do is iterate over faces,
    // then for each for each face, append N quads to the list (where N is the
    // degree of the face).  For this routine, it may be more convenient to simply
    // append quads to the end of the list (rather than allocating it ahead of
    // time), though YMMV.  You can of course iterate around a face by starting
    // with its first halfedge and following the "next" pointer until you get
    // back to the beginning.  The tricky part is making sure you grab the right
    // indices in the right order---remember that there are indices on vertices,
    // edges, AND faces of the original mesh.  All of these should get used.  Also
    // remember that you must have FOUR indices per face, since you are making a
    // QUAD mesh!

    HalfedgeIter face_start;
    HalfedgeIter curr, next;
    

    // TODO iterate over faces
    // TODO loop around face
    // TODO build lists of four indices for each sub-quad
    // TODO append each list of four indices to face list
    for (FaceIter f = facesBegin(); f != facesEnd(); f++) {
        face_start = f->halfedge();
        curr = face_start;
        do {
            vector<Index> quad;
            next = curr->next();

            quad.push_back(f->index);
            quad.push_back(curr->edge()->index);
            quad.push_back(next->vertex()->index);
            quad.push_back(next->edge()->index);

            subDFaces.push_back(quad);
            
            curr = curr->next();
        } while (curr != face_start);
    }


}

FaceIter HalfedgeMesh::bevelVertex(VertexIter v) {
  // TODO This method should replace the vertex v with a face, corresponding to
  // a bevel operation. It should return the new face.  NOTE: This method is
  // responsible for updating the *connectivity* of the mesh only---it does not
  // need to update the vertex positions.  These positions will be updated in
  // HalfedgeMesh::bevelVertexComputeNewPositions (which you also have to
  // implement!)

  showError("bevelVertex() not implemented.");
  return facesBegin();
}

FaceIter HalfedgeMesh::bevelEdge(EdgeIter e) {
  // TODO This method should replace the edge e with a face, corresponding to a
  // bevel operation. It should return the new face.  NOTE: This method is
  // responsible for updating the *connectivity* of the mesh only---it does not
  // need to update the vertex positions.  These positions will be updated in
  // HalfedgeMesh::bevelEdgeComputeNewPositions (which you also have to
  // implement!)

  showError("bevelEdge() not implemented.");
  return facesBegin();
}

FaceIter HalfedgeMesh::bevelFace(FaceIter f) {
  // TODO This method should replace the face f with an additional, inset face
  // (and ring of faces around it), corresponding to a bevel operation. It
  // should return the new face.  NOTE: This method is responsible for updating
  // the *connectivity* of the mesh only---it does not need to update the vertex
  // positions.  These positions will be updated in
  // HalfedgeMesh::bevelFaceComputeNewPositions (which you also have to
  // implement!)

    FaceIter middle_face = newFace();
    int degree = f->degree();

    HalfedgeIter start = f->halfedge();
    
    vector<VertexIter> old_vertices;
    HalfedgeIter curr = start;
    do {
        old_vertices.push_back(curr->vertex());
        curr = curr->next();
    } while (curr != start);

    assert(old_vertices.size() == f->degree());

    vector<VertexIter> v;
    vector<FaceIter> f_new;
    vector<EdgeIter> e_radial;
    vector<EdgeIter> e_parallel;
    vector<HalfedgeIter> f_halfedges;
    vector<HalfedgeIter> h_radial;
    vector<HalfedgeIter> h_radial_twins;
    vector<HalfedgeIter> h_parallel;
    vector<HalfedgeIter> h_parallel_twins;

    f_new.push_back(f);

    for (int i = 0; i < old_vertices.size(); i++) {
        curr = curr->next();
        f_halfedges.push_back(curr);

        v.push_back(newVertex());
        if (i > 0)
            f_new.push_back(newFace());
        e_radial.push_back(newEdge());
        e_parallel.push_back(newEdge());

        h_radial.push_back(newHalfedge());
        h_radial_twins.push_back(newHalfedge());

        h_parallel.push_back(newHalfedge());
        h_parallel_twins.push_back(newHalfedge());

    }

    for (int i = 0; i < old_vertices.size(); i++) {
        f_halfedges[i]->next() = h_radial[i];
        f_halfedges[i]->face() = f_new[i];

        v[i]->halfedge() = h_radial_twins[i];
        f_new[i]->halfedge() = h_radial[i];

        e_radial[i]->halfedge() = h_radial[i];
        e_parallel[i]->halfedge() = h_parallel[i];

        h_radial[i]->setNeighbors(h_parallel[i], h_radial_twins[i],
            f_halfedges[i]->twin()->vertex(), e_radial[i], f_new[i]);
        h_radial_twins[i]->setNeighbors(f_halfedges[(i + 1) % degree], 
            h_radial[i], v[i], e_radial[i], f_new[(i + 1) % degree]);

        h_parallel[i]->setNeighbors(h_radial_twins[(i + degree - 1) % degree],
            h_parallel_twins[i], v[i], e_parallel[i], f_new[i]);
        h_parallel_twins[i]->setNeighbors(h_parallel_twins[(i + 1) % degree],
            h_parallel[i], v[(i + degree - 1) % degree], e_parallel[i], middle_face);

    }

    middle_face->halfedge() = h_parallel_twins[0];

    checkConsistency();
    
    return middle_face;
}


void HalfedgeMesh::bevelFaceComputeNewPositions(
    vector<Vector3D>& originalVertexPositions,
    vector<HalfedgeIter>& newHalfedges, double normalShift,
    double tangentialInset) {
    // Compute new vertex positions for the vertices of the beveled face.
    //
    // These vertices can be accessed via newHalfedges[i]->vertex()->position for
    // i = 1, ..., newHalfedges.size()-1.
    //
    // The basic strategy here is to loop over the list of outgoing halfedges,
    // and use the preceding and next vertex position from the original mesh
    // (in the originalVertexPositions array) to compute an offset vertex
    // position.
    //
    // Note that there is a 1-to-1 correspondence between halfedges in
    // newHalfedges and vertex positions
    // in orig.  So, you can write loops of the form
    //
    // for( int i = 0; i < newHalfedges.size(); hs++ )
    // {
    //    Vector3D pi = originalVertexPositions[i]; // get the original vertex
    //    position correponding to vertex i
    // }
    //

    size_t deg = originalVertexPositions.size();
    
    Vector3D vec_next, vec_prev, vec_center;
    Vector3D f_center = 0;
    Vector3D f_normal;

    for (size_t i = 0; i < deg; i++) {
        f_center += originalVertexPositions[i];
    }
    f_center *= 1.0 / (double)deg;

    double angle;
    VertexIter curr_v;
    double shift;

    for (size_t i = 0; i < deg; i++) {

        // somehow my vertices' indexing ended up one off for whatever reason
        assert(originalVertexPositions[i] ==
            newHalfedges[(i + deg - 1) % deg]->twin()->vertex()->position);

        curr_v = newHalfedges[(i + deg - 1) % deg]->vertex();

        vec_next = originalVertexPositions[(i + 1) % deg] 
            - originalVertexPositions[i];
        vec_prev = originalVertexPositions[(i + deg - 1) % deg] 
            - originalVertexPositions[i];
        vec_center = f_center - curr_v->position;

        angle = asin(cross(vec_prev, vec_next).norm() / 
            (vec_next.norm() * vec_prev.norm()));
        f_normal = cross(vec_prev, vec_next)/cross(vec_prev, vec_next).norm();

        shift = (tangentialInset / sin(angle / 2.0)) / vec_center.norm();

        curr_v->position = originalVertexPositions[i];
        curr_v->position += shift * vec_center;
        curr_v->position += normalShift * f_normal;

    }

}

void HalfedgeMesh::bevelVertexComputeNewPositions(
    Vector3D originalVertexPosition, vector<HalfedgeIter>& newHalfedges,
    double tangentialInset) {
  // TODO Compute new vertex positions for the vertices of the beveled vertex.
  //
  // These vertices can be accessed via newHalfedges[i]->vertex()->position for
  // i = 1, ..., hs.size()-1.
  //
  // The basic strategy here is to loop over the list of outgoing halfedges,
  // and use the preceding and next vertex position from the original mesh
  // (in the orig array) to compute an offset vertex position.

}

void HalfedgeMesh::bevelEdgeComputeNewPositions(
    vector<Vector3D>& originalVertexPositions,
    vector<HalfedgeIter>& newHalfedges, double tangentialInset) {
  // TODO Compute new vertex positions for the vertices of the beveled edge.
  //
  // These vertices can be accessed via newHalfedges[i]->vertex()->position for
  // i = 1, ..., newHalfedges.size()-1.
  //
  // The basic strategy here is to loop over the list of outgoing halfedges,
  // and use the preceding and next vertex position from the original mesh
  // (in the orig array) to compute an offset vertex position.
  //
  // Note that there is a 1-to-1 correspondence between halfedges in
  // newHalfedges and vertex positions
  // in orig.  So, you can write loops of the form
  //
  // for( int i = 0; i < newHalfedges.size(); i++ )
  // {
  //    Vector3D pi = originalVertexPositions[i]; // get the original vertex
  //    position correponding to vertex i
  // }
  //

}

void HalfedgeMesh::splitPolygons(vector<FaceIter>& fcs) {
  for (auto f : fcs) splitPolygon(f);
}

void HalfedgeMesh::splitPolygon(FaceIter f) {
    // TODO: (meshedit) 
    // Triangulates a polygonal face recursively

    if (f->degree() == 3) return;

    HalfedgeIter start = f->halfedge();
    HalfedgeIter temp = start->next()->next();

    FaceIter new_face = newFace();
    EdgeIter new_edge = newEdge();
    HalfedgeIter new_h = newHalfedge();
    HalfedgeIter new_twin = newHalfedge();

    new_edge->halfedge() = new_h;
    start->next()->next() = new_h;
    new_h->setNeighbors(start, new_twin, temp->vertex(), new_edge, f);
    new_twin->setNeighbors(temp, new_h, start->vertex(), new_edge, new_face);

    while (temp->next() != start) {
        temp->face() = new_face;
        temp = temp->next();
    }

    assert(temp->next() == start);

    temp->face() = new_face;
    temp->next() = new_twin;

    new_face->halfedge() = temp->next()->next();

    checkConsistency();

    splitPolygon(new_face);
}

EdgeRecord::EdgeRecord(EdgeIter& _edge) : edge(_edge) {
  // TODO: (meshEdit)
  // Compute the combined quadric from the edge endpoints.
  // -> Build the 3x3 linear system whose solution minimizes the quadric error
  //    associated with these two endpoints.
  // -> Use this system to solve for the optimal position, and store it in
  //    EdgeRecord::optimalPoint.
  // -> Also store the cost associated with collapsing this edg in
  //    EdgeRecord::Cost.
}

void MeshResampler::upsample(HalfedgeMesh& mesh)
// This routine should increase the number of triangles in the mesh using Loop
// subdivision.
{
  // TODO: (meshEdit)
  // Compute new positions for all the vertices in the input mesh, using
  // the Loop subdivision rule, and store them in Vertex::newPosition.
  // -> At this point, we also want to mark each vertex as being a vertex of the
  //    original mesh.
  // -> Next, compute the updated vertex positions associated with edges, and
  //    store it in Edge::newPosition.
  // -> Next, we're going to split every edge in the mesh, in any order.  For
  //    future reference, we're also going to store some information about which
  //    subdivided edges come from splitting an edge in the original mesh, and
  //    which edges are new, by setting the flat Edge::isNew. Note that in this
  //    loop, we only want to iterate over edges of the original mesh.
  //    Otherwise, we'll end up splitting edges that we just split (and the
  //    loop will never end!)
  // -> Now flip any new edge that connects an old and new vertex.
  // -> Finally, copy the new vertex positions into final Vertex::position.

  // Each vertex and edge of the original surface can be associated with a
  // vertex in the new (subdivided) surface.
  // Therefore, our strategy for computing the subdivided vertex locations is to
  // *first* compute the new positions
  // using the connectity of the original (coarse) mesh; navigating this mesh
  // will be much easier than navigating
  // the new subdivided (fine) mesh, which has more elements to traverse.  We
  // will then assign vertex positions in
  // the new mesh based on the values we computed for the original mesh.

  // Compute updated positions for all the vertices in the original mesh, using
  // the Loop subdivision rule.

  // Next, compute the updated vertex positions associated with edges.

  // Next, we're going to split every edge in the mesh, in any order.  For
  // future
  // reference, we're also going to store some information about which
  // subdivided
  // edges come from splitting an edge in the original mesh, and which edges are
  // new.
  // In this loop, we only want to iterate over edges of the original
  // mesh---otherwise,
  // we'll end up splitting edges that we just split (and the loop will never
  // end!)

  // Finally, flip any new edge that connects an old and new vertex.

  // Copy the updated vertex positions to the subdivided mesh.
  showError("upsample() not implemented.");
}

void MeshResampler::downsample(HalfedgeMesh& mesh) {
  // TODO: (meshEdit)
  // Compute initial quadrics for each face by simply writing the plane equation
  // for the face in homogeneous coordinates. These quadrics should be stored
  // in Face::quadric
  // -> Compute an initial quadric for each vertex as the sum of the quadrics
  //    associated with the incident faces, storing it in Vertex::quadric
  // -> Build a priority queue of edges according to their quadric error cost,
  //    i.e., by building an EdgeRecord for each edge and sticking it in the
  //    queue.
  // -> Until we reach the target edge budget, collapse the best edge. Remember
  //    to remove from the queue any edge that touches the collapsing edge
  //    BEFORE it gets collapsed, and add back into the queue any edge touching
  //    the collapsed vertex AFTER it's been collapsed. Also remember to assign
  //    a quadric to the collapsed vertex, and to pop the collapsed edge off the
  //    top of the queue.
  showError("downsample() not implemented.");
}

void MeshResampler::resample(HalfedgeMesh& mesh) {
  // TODO: (meshEdit)
  // Compute the mean edge length.
  // Repeat the four main steps for 5 or 6 iterations
  // -> Split edges much longer than the target length (being careful about
  //    how the loop is written!)
  // -> Collapse edges much shorter than the target length.  Here we need to
  //    be EXTRA careful about advancing the loop, because many edges may have
  //    been destroyed by a collapse (which ones?)
  // -> Now flip each edge if it improves vertex degree
  // -> Finally, apply some tangential smoothing to the vertex positions
  showError("resample() not implemented.");
}

}  // namespace CMU462
