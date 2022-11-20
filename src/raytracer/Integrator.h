#pragma once

#include "UtilConcepts.h"
#include "Vector.h"
#include "Primitives.h"
#include "IntersectionTests.h"
#include "Material.h"

//"Non-Raytracing" variant of TraceRay function for testing
//simply returns albedo texture value of closest primitive of input ray, at the
//point of intersection
template<Float_t space_type, Float_t color_type>
Pixel<color_type> TraceRay(
    const Ray<space_type>& ray,
    const BVH<space_type,color_type>& bvh,
    int32_t remaining_bounces,
    const Sphere<space_type,color_type>& background) {
  //This will will be passed into the IntersectWithScene() and if it returns
  //true will contain the required information
  IntersectionResults<space_type,color_type> intersection_results;

  //not unused see
  remaining_bounces = remaining_bounces; // :)

  bool hit_a_primitive =
      IntersectBVH(ray,bvh,intersection_results);

  if (!hit_a_primitive) {
    //If no primitive is hit in BVH find intersection point of ray with background
    //sphere and retrieve the coresponding albedo pixel from the background texture
    return background.material.get()->GetAlbedoSphere(
        ray.direction,
        background.origin);
  }

  return intersection_results.hit_triangle ?
      intersection_results.hit_material.get()->GetAlbedoTri(
          intersection_results.intersection_point,
          intersection_results.hit_triangle.value()) :
      intersection_results.hit_material.get()->GetAlbedoSphere(
          intersection_results.intersection_point,
          intersection_results.hit_sphere_origin.value());
}