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
Pixel<color_type> TraceRayRaster(
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

template<Float_t space_type, Float_t color_type>
struct LocalSamplingParameters {
  Vec3<space_type> shading_normal;
  Vec3<space_type> geometry_normal;
  Vec3<space_type> inverted_ray_direction;
  space_type alpha;
  Pixel<color_type> F0_metal;
};

template<Float_t space_type, Float_t color_type>
struct LocalShadingParametersNoTransparency {
  LocalSamplingParameters<space_type,color_type> sampling_para;
  Pixel<color_type> local_emission;
  space_type alpha;
  space_type cos_angle;
};

template<Float_t space_type, Float_t color_type>
LocalShadingParameters<space_type,color_type> GetLocalShadingParameters(
    const IntersectionResults<space_type,color_type>& intersection_results,
    const Vec3<space_type>& incoming_ray_direction) {
  LocalShadingParameters<space_type,color_type> para;

  para.shading_normal = 
      GetNormal<space_type,color_type>(intersection_results);
  para.geometry_normal = intersection_results.hit_triangle ? 
      intersection_results.hit_triangle.value().geometry_normal :
      intersection_results.intersection_point -
      intersection_results.hit_sphere_origin.value();
  
  para.ray_originates_inside_medium =
      incoming_ray_direction * para.geometry_normal < 0;
  //Flip normals if inside medium
  if (para.ray_originates_inside_medium) {
    para.shading_normal = -1 * para.shading_normal;
    para.geometry_normal = -1 * para.geometry_normal;
  }

  para.cos_angle = incoming_ray_direction * para.shading_normal;

  space_type shading_angle = std::acos(para.cos_angle);
  para.local_emission = 
      GetEmitance<space_type,color_type>(intersection_results, shading_angle);
  para.transmissity =
      GetTransmissity<space_type,color_type>(intersection_results);
  
  space_type refractive_index_medium = 
      GetRefractiveIndex<space_type,color_type>(intersection_results);
  para.refractive_index_1 =
      para.ray_originates_inside_medium ? refractive_index_medium : 1.0;
  para.refractive_index_2 =
      para.ray_originates_inside_medium ? 1.0 : refractive_index_medium;
  
  //Compute F0 depending on whether the ray comes from inside of the medium or
  //hits the medium. The other medium is assumed to be air with n roughly = 1.0.
  
  space_type tmp = (para.refractive_index_1 - para.refractive_index_2) /
                   (para.refractive_index_1 + para.refractive_index_2);
  para.F0_refractive = tmp*tmp;

  Pixel<color_type> albedo =
      GetAlbedo<space_type,color_type>(intersection_results);
  space_type metalness =
      GetMetalness<space_type,color_type>(intersection_results);
  constexpr space_type min_dielectrics_F0 = 0.04;
  para.F0_metal =
    Pixel<color_type>(min_dielectrics_F0,min_dielectrics_F0,min_dielectrics_F0) *
    (1.0 - metalness) + (albedo * metalness);

  //alpha = roughness**2
  para.alpha = GetRoughness<space_type,color_type>(intersection_results);
  para.alpha *= para.alpha;

  return para;
}

template<Float_t space_type, Float_t color_type>
Pixel<color_type> TraceRayNoTransparency(
    const Ray<space_type>& ray,
    const BVH<space_type,color_type>& bvh,
    int32_t remaining_bounces,
    RNGLutEngine<space_type>& rng_engine,
    const std::optional<Sphere<space_type,color_type>>& background = std::nullopt) {
  //Max recursion depth reached
  if (remaining_bounces == 0) {
    return Pixel<color_type>(0.0,0.0,0.0);
  }

  //This will be passed into the IntersectWithScene() and if it returns
  //true will contain the required information
  IntersectionResults<space_type,color_type> intersection_results;

  bool hit_a_primitive =
      IntersectBVH<space_type,color_type>(ray,bvh,intersection_results);
  
  //Rays that did not intersect geometry return background emittence if
  //a background was supplied
  if (!hit_a_primitive) {
    
    return background ? 
        background.value().material.get()->GetEmitanceSphere(
            ray.direction,
            background.value().origin,
            0.0) :
        Pixel<color_type>(0.0,0.0,0.0);
  }
  
  //Convenience var; incoming ray dir but points away from geometry
  Vec3<space_type> ray_dir_inv = ray.direction;
  ray_dir_inv = ray_dir_inv * (-1.0f);

  LocalShadingParameters<space_type,color_type> local_para = 
    GetLocalShadingParameters<space_type,color_type>(
      intersection_results,
      ray_dir_inv);
  
  //Flip in case that we approach primitive from below (shading and geometry
  //normal) are already adjusted by GetLocalShadingParameters()
  ray_dir_inv = local_para.ray_originates_inside_medium ? 
      -1.0f * ray_dir_inv : ray_dir_inv;

  Pixel<color_type> attenuence;
  Ray<space_type> bounce_ray;
  bounce_ray.origin = ray.origin;

  //Sample reflection ray.
  bool spawned_reflection_ray = SampleRay(
      local_para.shading_normal,
      local_para.geometry_normal,
      ray_dir_inv,
      local_para.alpha,
      local_para.F0_metal,
      attenuence,
      bounce_ray.direction,
      rng_engine);
  
  //Can happen when the resulting ray is emitted below shading plane but not
  //geometry plane or incoming ray is below shading normal
  if (!spawned_reflection_ray) {
    return Pixel<color_type>(0.0,0.0,0.0);
  }

  return local_para.local_emission +
         attenuence * TraceRayNoTransparency(
            bounce_ray,bvh,remaining_bounces-1,rng_engine,background);
}