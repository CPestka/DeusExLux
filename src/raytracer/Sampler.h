
//Is called after an intersection to compute the new outgoing ray and the
//acording attenuence
//We use:
//-Importance sampling
//-Microfacet model with GGX distribution
//-GGX height corection
//-Metalness workflow
template<Float_t space_type, Float_t color_type>
bool SampleRay(const Vec3<space_type>& shading_normal,
               const Vec3<space_type>& geometry_normal,
               const Vec3<space_type>& incomming_ray_direction,
               space_type alpha, //alpha = roughness**2
               //Fresnel coef computed from basecolor and metallness
               const Pixel<color_type>& F0,
               Pixel<color_type>& attenuence,
               Vec3<space_type>& out_going_ray_direction){
  //rotate incoming ray to plane that has the shading normal as its normal
  //return if ray comes from below the "shading plane"
  if (shading_normal * incomming_ray_direction <= 0.0) {
    return false;
  }
  Vec3<space_type> ray_in_dir_shading_space;
  space_type angle_shading_normal = std::acos(shading_normal*Vec3<space_type>(0.0,1.0,0.0));
  Quaternion<space_type> quat;
  Quaternion<space_type> quat_inv;
  if (fabs(angle_shading_normal) <= 0.01) {
    ray_in_dir_shading_space = incomming_ray_direction;
  } else {
    quat = {ComputeCrossProduct(shading_normal, {0.0,1.0,0.0}),
            angle_shading_normal};
    quat_inv = quat.GetInverted();
    ray_in_dir_shading_space =
        Vec3<space_type>(((quat_inv *
            Quaternion<space_type>(incomming_ray_direction, 0.0)) * quat));
  }
  
  Vec3<space_type> normal_shading_space = {0.0,1.0,0.0};
  space_type alpha_squared = alpha * alpha;


  //Sample Half vector via VNDF method
  Vec3<space_type> tmp1(alpha*ray_in_dir_shading_space.data[0],
                        (-1)*alpha*ray_in_dir_shading_space.data[2],
                        ray_in_dir_shading_space.data[1]);
  tmp1.Normalize();
  space_type tmp2 = tmp1.data[0]*tmp1.data[0] + tmp1.data[1]*tmp1.data[1];
  Vec3<space_type> tmp3 = tmp2 > 0.0 ? Vec3<space_type>(-tmp1.data[1],tmp1.data[0],0.0) :
                                       Vec3<space_type>(1.0,0.0,0.0);
  Vec3<space_type> tmp4 = ComputeCrossProduct(tmp1, tmp3);

  space_type r = std::sqrt(GetUniformRandom<space_type>(0.0,false,1.0,false));
  space_type phi = GetUniformRandom<space_type>(0.0,false,2.0*M_PI,false);
  space_type tmp5 = r * std::cos(phi);
  space_type tmp6 = r * std::sin(phi);
  space_type tmp7 = 0.5 * (1.0 + tmp1.data[3]);
  tmp6 = std::sqrt(1.0 - tmp5 * tmp5) * (1.0 - tmp7) + tmp6 * tmp7;

  Vec3<space_type> tmp8 =
      tmp5 * tmp3 +
      tmp6 * tmp4 +
      std::sqrt(max(0.0, 1.0 - tmp5 * tmp5 - tmp6 * tmp6)) * tmp1;

  Vec3<space_type> tmp9 =
      {alpha * tmp8.data[0], alpha * tmp8.data[1], max(0.0, tmp8.data[2])}; 
  tmp9.Normalize();
  Vec3<space_type> half_vec = {tmp9.data[0], tmp9.data[2], (-1)*tmp9.data[1]};
  space_type normal_dot_half_vec = normal_shading_space * half_vec;


  //Reflect incomming ray on plane for which the half vector is the normal
  //i.e. reflecht on microfacet plane
  Vec3<space_type> ray_out_dir =
    ray_in_dir_shading_space - 2.0 * normal_dot_half_vec * normal_shading_space;
  space_type half_vec_dot_ray_out_dir = half_vec * ray_out_dir;

  //Rotate ray back with previously computed quaternions
  if (angle_shading_normal <= 0.01) {
    out_going_ray_direction = ray_out_dir;
  } else {
    out_going_ray_direction =
        Vec3<space_type>(((quat *
            Quaternion<space_type>(ray_out_dir, 0.0)) * quat_inv));
  }

  //Reject rays that would be reflected below geometry plane
  //this can happen in the case of a shallow reflection and the geometry normal
  //and the shading normal being significantly different
  if (out_going_ray_direction * geometry_normal <= 0.0) {
    return false;
  }
  
  //Precompute frequently used dot products
  space_type normal_dot_ray_out_dir = normal_shading_space * ray_out_dir;
  space_type normal_dot_ray_in_dir =
      normal_shading_space * ray_in_dir_shading_space;
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
               (F0 + (1.0 - F0) * std::pow(1.0 - half_vec_dot_ray_out_dir, 5));
  return true;
}


template<Float_t space_type, Float_t color_type>
Pixel<color_type> GetEmissiveContribution(
    IntersectionResults<space_type,color_type>& result,
    space_type angle) {

  Pixel<color_type> emissiveness;                                      
  if (!result.hit_material.emmittence_map) {
    emissiveness = result.hit_material.emmittence;
  } else {
    if (result.hit_triangle) {
      emissiveness = 
          result.hit_material.emmittence_map.value().GetPixel(
              result.hit_matrial.emittence_uv_tri,
              result.hit_triangle.value(),
              result.intersection_point);
    } else {
      emissiveness = 
          result.hit_material.emmittence_map.value().GetPixel(
              result.intersection_point,
              result.hit_material.emittence_sphere_ey.value());
    }
  }

  if (result.hit_material.emit_type == COSINE) {
    return std::cos(angle) * emissiveness;
  }
  if (result.hit_material.emit_type == DISCRETE_BEAM_5) {
    return angle < 2 * M_PI * (5 / 360) ? emissiveness : Vec3<color_type>(0.0,0.0,0.0);
  }
  if (result.hit_material.emit_type == DISCRETE_BEAM_10) {
    return angle < 2 * M_PI * (10 / 360) ? emissiveness : Vec3<color_type>(0.0,0.0,0.0);
  }
  if (result.hit_material.emit_type == DISCRETE_BEAM_30) {
    return angle < 2 * M_PI * (30 / 360) ? emissiveness : Vec3<color_type>(0.0,0.0,0.0);
  }
  
  return emissiveness;
}

template<Float_t space_type, Float_t color_type>
Vec3<space_type> GetShadingNormal(
    IntersectionResults<space_type,color_type>& result,
    Vec3<space_type>& geometry_normal) {

  if (result.hit_material.normal_map) {
    if (result.hit_triangle) {
      return result.hit_material.normal_map.value().GetPixel(
          result.hit_material.normal_uv_tri.value(),
          result.hit_triangle.value(),
          result.intersection_point);
    } else {
      return result.hit_material.normal_map.value().GetPixel(
          result.intersection_point,
          result.hit_material.normal_sphere_ey.value());
    }
  } else {
    return geometry_normal;
  }
}

template<Float_t space_type, Float_t color_type>
space_type GetAlpha(
    IntersectionResults<space_type,color_type>& result) {

  if (result.hit_material.roughness_map) {
    if (result.hit_triangle) {
      space_type tmp = result.hit_material.roughness_map.value().GetPixel(
          result.hit_material.roughness_uv_tri.value(),
          result.hit_triangle.value(),
          result.intersection_point);
      return tmp*tmp;
    } else {
      space_type tmp = result.hit_material.roughness_map.value().GetPixel(
          result.intersection_point,
          result.hit_material.roughness_sphere_ey.value());
      return tmp*tmp;
    }
  } else {
    return result.hit_material.roughness * result.hit_material.roughness;
  }
}

template<Float_t space_type, Float_t color_type>
Pixel<color_type> GetFresnel0(
    IntersectionResults<space_type,color_type>& result) {
  constexpr space_type min_dielectrics_f0 = 0.04;
  Pixel<color_type> albedo;
  
  if (result.hit_material.albedo_map) {
    if (result.hit_triangle) {
      albedo = result.hit_material.albedo_map.value().GetPixel(
          result.hit_material.albedo_uv_tri.value(),
          result.hit_triangle.value(),
          result.intersection_point);
    } else {
      albedo = result.hit_material.albedo_map.value().GetPixel(
          result.intersection_point,
          result.hit_material.albedo_sphere_ey.value());
    }
  } else {
    return result.hit_material.albedo * result.hit_material.albedo;
  }

  space_type metalness;
  
  if (result.hit_material.metalness_map) {
    if (result.hit_triangle) {
      albedo = result.hit_material.metalness_map.value().GetPixel(
          result.hit_material.metalness_uv_tri.value(),
          result.hit_triangle.value(),
          result.intersection_point);
    } else {
      albedo = result.hit_material.metalness_map.value().GetPixel(
          result.intersection_point,
          result.hit_material.metalness_sphere_ey.value());
    }
  } else {
    return result.hit_material.metalness * result.hit_material.metalness;
  }

  return ((Pixel<color_type>(min_dielectrics_f0,min_dielectrics_f0,min_dielectrics_f0) *
           (1 - metalness)) + (albedo * metalness));
}



template<Float_t space_type, Float_t color_type>
Pixel<color_type> TraceRay(
    const Ray<space_type>& ray,
    const BVH<space_type,color_type>& bvh,
    int32_t remaining_bounces,
    const std::optional<Texture<Pixel<color_type>,
                        space_type,color_type>>& world_hdr = std::nullopt) {
  //Max recursion depth reached
  if (remaining_bounces == 0) {
    return {0.0,0.0,0.0};
  }

  //This will will be passed into the IntersectWithScene() and if it returns
  //true will contain the required information
  IntersectionResults<space_type,color_type> intersection_results;

  bool hit_a_primitive =
      IntersectWithScene(ray,bvh,intersection_results);

  if (!hit_a_primitive) {
    if (world_hdr) { 
      //In case there is a background HDR as a global emissive object
      return world_hdr.value().GetPixel(Vec3<space_type>(ray.direction), Vec3<space_type>(0.0,1.0,0.0));
    } else { 
      //In case there isnt -> nothing was hit -> no emission -> black
      return {0.0,0.0,0.0};
    }
  }

  space_type ray_primitive_angle =
      std::acos(intersection_results.geometry_normal*ray.direction);

  Pixel<color_type> attenuence;
  Ray<space_type> bounce_ray;
  //Computes and sets attenuence and bounce_ray
  SampleRay<space_type,color_type>(
            GetShadingNormal(intersection_results,
                             intersection_results.geometry_normal),
            intersection_results.geometry_normal,
            ray.direction,
            GetAlpha(intersection_results),
            GetFresnel0(intersection_results),
            attenuence,
            bounce_ray.direction);
  
  //Recursive call for next bounce
  return GetEmissiveContribution(intersection_results, ray_primitive_angle) +
         attenuence * TraceRay(bounce_ray,
                               bvh,
                               remaining_bounces - 1);
}