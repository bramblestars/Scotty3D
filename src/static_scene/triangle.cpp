#include "triangle.h"

#include "CMU462/CMU462.h"
#include "GL/glew.h"

namespace CMU462 {
namespace StaticScene {

Triangle::Triangle(const Mesh* mesh, vector<size_t>& v) : mesh(mesh), v(v) {}
Triangle::Triangle(const Mesh* mesh, size_t v1, size_t v2, size_t v3)
    : mesh(mesh), v1(v1), v2(v2), v3(v3) {}

BBox Triangle::get_bbox() const {
  // TODO (PathTracer):
  // compute the bounding box of the triangle

  return BBox();
}

bool Triangle::intersect(const Ray& r) const {
    // TODO (PathTracer): implement ray-triangle intersection`

    // Procedure and variable names from wiki

    Vector3D p1 = mesh->positions[v1];
    Vector3D p2 = mesh->positions[v2];
    Vector3D p3 = mesh->positions[v3];

    Vector3D e1 = p2 - p1;
    Vector3D e2 = p3 - p1;
    Vector3D s = r.o - p1;

    if (dot(cross(e1, r.d), e2) == 0)
        return false;

    double coeff = 1 / (dot(cross(e1, r.d), e2));

    double u = coeff * (-1.0 * (dot(cross(s, e2), r.d)));
    double v = coeff * (dot(cross(e1, r.d), s));
    double t = coeff * (-1.0 * (dot(cross(s, e2), e1)));

    if (u > 0 && v > 0 && u + v < 1.0) {
        return true;
    }

    return false;
}

bool Triangle::intersect(const Ray& r, Intersection* isect) const {
    // TODO (PathTracer):
    // implement ray-triangle intersection. When an intersection takes
    // place, the Intersection data should be updated accordingly

    if (!intersect(r)) {
        return false;
    }

    else {
        Vector3D p1 = mesh->positions[v1];
        Vector3D p2 = mesh->positions[v2];
        Vector3D p3 = mesh->positions[v3];

        Vector3D e1 = p2 - p1;
        Vector3D e2 = p3 - p1;
        Vector3D s = r.o - p1;

        Vector3D normal = cross(e1, e2);
        normal.normalize();

        double coeff = 1 / (dot(cross(e1, r.d), e2));

        isect->n = normal;
        isect->primitive = this;
        isect->t = coeff * (-1.0 * (dot(cross(s, e2), e1)));
        isect->bsdf = mesh->get_bsdf();

        return true;
    }
}

void Triangle::draw(const Color& c) const {
  glColor4f(c.r, c.g, c.b, c.a);
  glBegin(GL_TRIANGLES);
  glVertex3d(mesh->positions[v1].x, mesh->positions[v1].y,
             mesh->positions[v1].z);
  glVertex3d(mesh->positions[v2].x, mesh->positions[v2].y,
             mesh->positions[v2].z);
  glVertex3d(mesh->positions[v3].x, mesh->positions[v3].y,
             mesh->positions[v3].z);
  glEnd();
}

void Triangle::drawOutline(const Color& c) const {
  glColor4f(c.r, c.g, c.b, c.a);
  glBegin(GL_LINE_LOOP);
  glVertex3d(mesh->positions[v1].x, mesh->positions[v1].y,
             mesh->positions[v1].z);
  glVertex3d(mesh->positions[v2].x, mesh->positions[v2].y,
             mesh->positions[v2].z);
  glVertex3d(mesh->positions[v3].x, mesh->positions[v3].y,
             mesh->positions[v3].z);
  glEnd();
}

}  // namespace StaticScene
}  // namespace CMU462
