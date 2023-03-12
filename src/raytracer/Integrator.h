#pragma once

#include <cmath>
#include <optional>

#include "UtilConcepts.h"
#include "Vector.h"
#include "Primitives.h"
#include "IntersectionTests.h"
#include "Material.h"

//"Non-Raytracing" variant of TraceRay function for testing
//simply returns albedo texture value of closest primitive of input ray, at the
//point of intersection
template<Float_t space_type, Float_t color_type>
Pixel<color_type> TraceRasterRay(
    const Ray<space_type>& ray,
    const BVH<space_type,color_type>& bvh,
    const std::optional<Sphere<space_type,color_type>>& background) {
  //This will will be passed into the IntersectWithScene() and if it returns
  //true will contain the required information
  IntersectionResults<space_type,color_type> intersection_results;

  bool hit_a_primitive =
      IntersectBVH(ray,bvh,intersection_results);

  if (!hit_a_primitive) {
    //If no primitive is hit in BVH find intersection point of ray with background
    //sphere and retrieve the coresponding albedo pixel from the background texture
    return background ? background.value().material.get()->GetAlbedoSphere(
        ray.direction,
        background.value().origin) :
        Pixel<color_type>(0.0,0.0,0.0);
  }

  return GetAlbedo<space_type,color_type>(intersection_results);
}