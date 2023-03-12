#pragma once

#include <cmath>
#include <algorithm>

#include "UtilConcepts.h"
#include "Vector.h"
#include "Primitives.h"
#include "RNG.h"

//Is called after an intersection to compute the new outgoing ray and the
//acording attenuence
//We use:
//-Importance sampling
//-Microfacet model with GGX distribution
//-GGX height corection
//-Metalness workflow
//
//This function assumes as is common that all normals and rays point away from
//the plane of the geometry
//Also assumes that shading_normal * incoming_ray_direction <= 0
template<Float_t space_type, Float_t color_type>
bool SampleRay(const LocalSamplingParameters<space_type,color_type>& para, 
               Pixel<color_type>& attenuence,
               Vec3<space_type>& out_going_ray_direction,
               RNGLutEngine<space_type>& rng_engine){
  //The first task is to rotate the ray into the shading space defined by the
  //shading normal
  Vec3<space_type> ray_in_shading_space;

  Vec3<space_type> e_y = Vec3<space_type>(0.0,1.0,0.0);
  space_type e_y_dot_shading_normal = e_y * para.shading_normal;
  
  //Quaterions used for transform into and out of shading space
  Quaternion<space_type> look_at_quat{};
  Quaternion<space_type> inverse_look_at_quat{};

  if (fabs(e_y_dot_shading_normal) <= 0.001) {
    //Dont bother to do rotation if already close enough
    ray_in_dir_shading_space = para.inverted_ray_direction;
  } else {
    space_type angle_shading_normal = std::acos(e_y_dot_shading_normal);

    look_at_quat = Quaternion<space_type>(
        ComputeCrossProduct<space_type,space_type>(para.shading_normal, e_y),
        angle_shading_normal);
    
    inverse_look_at_quat = look_at_quat;
    inverse_look_at_quat.Invert();

    ray_in_dir_shading_space =
        Vec3<space_type>(((inverse_look_at_quat) *
            Quaternion<space_type>(para.inverted_ray_direction, 0.0)) * look_at_quat);
  }
  
  space_type alpha_squared = para.alpha * para.alpha;

  //Sample Half vector via VNDF method
  Vec3<space_type> tmp1(para.alpha*ray_in_dir_shading_space.x,
                        (-1)*para.alpha*ray_in_dir_shading_space.z,
                        ray_in_dir_shading_space.y);
  tmp1.Normalize();
  space_type tmp2 = tmp1.x*tmp1.x + tmp1.y*tmp1.y;
  Vec3<space_type> tmp3 = tmp2 > 0.0 ? Vec3<space_type>(-tmp1.y,tmp1.x,0.0) :
                                       Vec3<space_type>(1.0,0.0,0.0);
  Vec3<space_type> tmp4 = ComputeCrossProduct(tmp1, tmp3);

  space_type r = std::sqrt(rng_engine.GetUniformRandom(0.0,false,1.0,false));
  space_type phi = rng_engine.GetUniformRandom(0.0,false,2.0*M_PI,false);
  space_type tmp5 = r * std::cos(phi);
  space_type tmp6 = r * std::sin(phi);
  space_type tmp7 = 0.5 * (1.0 + tmp1.z);
  tmp6 = std::sqrt(1.0 - tmp5 * tmp5) * (1.0 - tmp7) + tmp6 * tmp7;

  Vec3<space_type> tmp8 =
      tmp5 * tmp3 +
      tmp6 * tmp4 +
      std::sqrt(std::max<space_type>(0.0, 1.0 - tmp5 * tmp5 - tmp6 * tmp6)) * tmp1;

  Vec3<space_type> tmp9 =
      {para.alpha * tmp8.x, para.alpha * tmp8.y, std::max<space_type>(0.0, tmp8.z)}; 
  tmp9.Normalize();
  Vec3<space_type> half_vec = {tmp9.x, tmp9.z, (-1)*tmp9.y};
  space_type normal_dot_half_vec = e_y * half_vec;


  //Reflect incomming ray on plane for which the half vector is the normal
  //i.e. reflecht on microfacet plane
  Vec3<space_type> ray_out_dir =
    ray_in_dir_shading_space - 2.0 * normal_dot_half_vec * e_y;
  space_type half_vec_dot_ray_out_dir = half_vec * ray_out_dir;

  //Rotate ray back from shading space to intial
  if (fabs(e_y_dot_shading_normal) <= 0.001) {
    //Again, dont bother if already close enough
    out_going_ray_direction = ray_out_dir;
  } else {
    look_at_quat.x1 = -look_at_quat.x1;
    inverse_look_at_quat.x1 = -inverse_look_at_quat.x1;
    out_going_ray_direction =
        Vec3<space_type>((inverse_look_at_quat *
            Quaternion<space_type>(ray_out_dir, 0.0)) * look_at_quat);
  }

  //Reject rays that would be reflected below geometry plane
  //this can happen in the case of a shallow reflection and the geometry normal
  //and the shading normal being significantly different
  if (out_going_ray_direction * para.geometry_normal <= 0.0) {
    return false;
  }
  
  space_type normal_dot_ray_out_dir = e_y * ray_out_dir;
  space_type normal_dot_ray_in_dir =
      e_y * ray_in_dir_shading_space;
  space_type normal_dot_ray_in_dir_squared =
      normal_dot_ray_in_dir * normal_dot_ray_in_dir;
  space_type normal_dot_ray_out_dir_squared =
      normal_dot_ray_out_dir * normal_dot_ray_out_dir;
  
  //Compute sample weight using optimized formulation of height correlated GGX
  //calculation
  space_type G1_ray_in = 
      2.0 / (1.0 + std::sqrt((alpha_squared *
      (1.0 - normal_dot_ray_in_dir_squared) + normal_dot_ray_in_dir_squared) /
      normal_dot_ray_in_dir_squared));
  space_type G1_ray_out =
      2.0 / (1.0 + std::sqrt((alpha_squared *
      (1.0 - normal_dot_ray_out_dir_squared) + normal_dot_ray_out_dir_squared) /
      normal_dot_ray_out_dir_squared));
  space_type sample_weight =
      G1_ray_in / (G1_ray_in + G1_ray_out - G1_ray_in * G1_ray_out);

  //Compute Attenuence from sample weight and Fresenel contribution
  //Fresenel contribution computed via Fresnel-Schlick approximation
  attenuence = sample_weight * 
               (para.F0_metal + (1.0 - para.F0_metal) *
                std::pow(1.0 - half_vec_dot_ray_out_dir, 5));
  return true;
}