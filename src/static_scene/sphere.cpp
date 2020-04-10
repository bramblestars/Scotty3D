#include "sphere.h"

#include <cmath>

#include "../bsdf.h"
#include "../misc/sphere_drawing.h"

namespace CMU462 {
namespace StaticScene {

bool Sphere::test(const Ray& r, double& t1, double& t2) const {
  // TODO (PathTracer):
  // Implement ray - sphere intersection test.
  // Return true if there are intersections and writing the
  // smaller of the two intersection times in t1 and the larger in t2.
  
  return false;
}

bool Sphere::intersect(const Ray& R) const {
    // TODO (PathTracer):
    // Implement ray - sphere intersection.
    // Note that you might want to use the the Sphere::test helper here.

    // From Wikipedia - "Line-Sphere Intersection"
    Vector3D c = o; //origin of sphere
    double under_root = (dot(R.d, R.o - c)) * (dot(R.d, R.o - c))
        - (((R.o - c).norm()) * ((R.o - c).norm()) - r * r);

    return under_root > 0;
}

bool Sphere::intersect(const Ray& R, Intersection* isect) const {
    // TODO (PathTracer):
    // Implement ray - sphere intersection.
    // Note again that you might want to use the the Sphere::test helper here.
    // When an intersection takes place, the Intersection data should be updated
    // correspondingly.

    Vector3D c = o; //origin of sphere
    double under_root = (dot(R.d, R.o - c)) * (dot(R.d, R.o - c))
        - (((R.o - c).norm()) * ((R.o - c).norm()) - r * r);

    double t = (-2 * (dot(R.d, R.o - c))) - sqrt(under_root);
    t /= (2 * R.d.norm2());

    Vector3D normal = (R.o + t * R.d) - c;
    normal.normalize();

    if (under_root > 0) {
        isect->n = normal;
        isect->primitive = this;
        isect->t = t;
        isect->bsdf = get_bsdf();
    }

    return under_root > 0;
}

void Sphere::draw(const Color& c) const { Misc::draw_sphere_opengl(o, r, c); }

void Sphere::drawOutline(const Color& c) const {
  // Misc::draw_sphere_opengl(o, r, c);
}

}  // namespace StaticScene
}  // namespace CMU462
