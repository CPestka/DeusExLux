#pragma once
//Code that is responsible for the intersection testing with the primitives,
//generation of bounce rays and calculation of the color for the given chain of
//rays

#include <optional>
#include <vector>
#include <limits>

#include "UtilConcepts.h"
#include "BVH.h"
#include "Primitives.h"
#include "Mesh.h"
#include "Material.h"
#include "RNG.h"
#include "Vector.h"
#include "Image.h"

//Helper struct used to return Intersection test results
template<Float_t space_type, Float_t color_type>
struct IntersectionResults {
  Vec3<space_type> intersection_point;
  Vec3<space_type> geometry_normal; // NEEDED????
  std::optional<Tri<space_type,color_type>> hit_triangle;
  std::optional<Vec3<space_type>> hit_sphere_origin;
  std::shared_ptr<PBRMaterial<space_type, color_type>> hit_material;
};

//Möller-Trumbore intersection algorithm
template<Float_t space_type, Float_t color_type>
void IntersectTriangle(const Tri<space_type,color_type>& tri,
                       const Ray<space_type>& ray,
                       space_type& distance_of_closest_hit_untill_now,
                       IntersectionResults<space_type,color_type>& results){
  constexpr space_type epsilon = 10e-7;
  Vec3<space_type> edge_0 = tri.vertex[1] - tri.vertex[0];
  Vec3<space_type> edge_1 = tri.vertex[2] - tri.vertex[0];
  Vec3<space_type> ray_dir_cross_edge_1 =
      ComputeCrossProduct(ray.direction, edge_1);
  space_type a = edge_0 * ray_dir_cross_edge_1;
  
  if (a > epsilon && a < epsilon) {
    return;
  }
  
  space_type b = 1.0 / a;
  Vec3<space_type> c = ray.origin - tri.vertex[0];
  space_type d = b * (c * ray_dir_cross_edge_1);

  if (d < 0.0 || d > 1.0) {
    return;
  }

  Vec3<space_type> e = ComputeCrossProduct(c, edge_0);
  space_type f = b * (ray.direction * e);

  if (f < 0.0 || f + d > 1.0) {
    return;
  }

  space_type dist = b * (edge_1 * e);
  
  if (dist > epsilon && dist < distance_of_closest_hit_untill_now) {
    //If hit is actually the closest yet modify the hit parameter
    distance_of_closest_hit_untill_now = dist;
    results.intersection_point = (dist * ray.direction) + ray.origin;
    results.hit_material = tri.material;
    results.geometry_normal = tri.geometry_normal;
    results.hit_triangle = tri;
  }
}

template<Float_t space_type, Float_t color_type>
void IntersectSphere(const Sphere<space_type,color_type>& sph,
                     const Ray<space_type>& ray,
                     space_type& distance_of_closest_hit_untill_now,
                     IntersectionResults<space_type,color_type>& results){
  constexpr space_type epsilon = 10e-7;
  Vec3<space_type> a = ray.origin - sph.origin;
  space_type b = ray.direction * ray.direction;
  space_type c = a * ray.direction;
  space_type d = (a * a) - (sph.r * sph.r);
  space_type discriminant = (c * c) - (b * d);
  if (discriminant > 0) {
    space_type dist = (-c - std::sqrt(discriminant)) / b;
    if (dist > epsilon && dist < distance_of_closest_hit_untill_now) {
      //If hit is actually the closest yet modify the hit parameter
      distance_of_closest_hit_untill_now = dist;
      results.intersection_point = (dist * ray.direction) + ray.origin;
      results.hit_material = sph.material;
      results.geometry_normal = sph.ComputeSurfaceNormal(results.intersection_point);
      results.hit_sphere_origin = sph.origin;
    }
  }
}

//Utility function to wrap intersect tests of triangle of a bounding box and
//return only the distance of the closest intersection point, as only that is
//needed fopr the decision whether to enter the bvh of the bounding box for
//further intersection tests
template<Float_t space_type, Float_t color_type>
space_type IntersectBox(const Box<space_type,color_type>& bounding_box,
                        const Ray<space_type>& ray) {
  space_type tri_dist = std::numeric_limits<space_type>::max();
  IntersectionResults<space_type,color_type> intersection_result;
  
  for(int i=0; i<12; i++){
    IntersectTriangle(bounding_box.tris[i],
                      ray,
                      tri_dist,
                      intersection_result);
  }
  return tri_dist;
}

template<Float_t space_type,
         Float_t color_type>
class BVHBoxHelper {
  public:
  BVHNode<space_type,color_type>* bvh_node = nullptr;
  space_type dist_to_bounding_box = 0;

  // BVHBoxHelper() {};
  // BVHBoxHelper(
  //     const BVHBoxHelper<space_type,color_type>& other) :
  //         bvh_node(other.bvh_node),
  //         dist_to_bounding_box(other.dist_to_bounding_box) {};
  // BVHBoxHelper(
  //     const BVHNode<space_type,color_type>* bvh_node,
  //     space_type dist) :
  //     bvh_node( bvh_node),
  //     dist_to_bounding_box(dist) {};
  // BVHBoxHelper<space_type,color_type>& operator=(
  //     const BVHBoxHelper<space_type,color_type>& other) :
  //         bvh_node(other.bvh_node),
  //         dist_to_bounding_box(other.dist_to_bounding_box) {};
};

template<Float_t space_type,
         Float_t color_type>
bool TestBVHLvl(
      const Ray<space_type>& ray,
      const BVHNode<space_type,color_type>& bvh_node,
      space_type& distance_of_closest_hit_untill_now,
      IntersectionResults<space_type,color_type>& results){
  //Only leaf bvh have triangles in it -> non leaf bvh wont set distance_of_..
  //in these two loops here

  //Test against triangles in this volume
  for(size_t i=0; i<bvh_node.tri_count; i++){
    IntersectTriangle(bvh_node.tris[i],
                      ray,
                      distance_of_closest_hit_untill_now,
                      results);
  }

  //Test against spheres in this volume
  for(size_t i=0; i<bvh_node.sphere_count; i++){
    IntersectSphere(bvh_node.spheres[i],
                    ray,
                    distance_of_closest_hit_untill_now,
                    results);
  }

  //Test against other bvhs in this lvl -> traverse recursively untill leave bvh
  //is hit and triangles and sphere intersects start
  size_t amount_of_child_bvh = bvh_node.child_count;
  std::vector<BVHBoxHelper<space_type,color_type>> test_priority_list;
  
  //Test all boxes of child bvhs and put them in vec in sorted form
  //Non-Hits are indicated by dist= int_max and are thus at the end of the list
  for(size_t i=0; i<amount_of_child_bvh; i++){
    space_type dist_to_box_intersection =
        IntersectBox(*((*(bvh_node.child_ptrs[i])).bounding_box.get()), ray);
    
    //Insert child box in priority list 
    for(size_t j=0; j<test_priority_list.size(); j++){
      if (test_priority_list[j].dist_to_bounding_box > dist_to_box_intersection) {
        test_priority_list.insert(
            std::next(test_priority_list.begin(), j),
            BVHBoxHelper<space_type,color_type>(
                bvh_node.child_ptrs[i],
                dist_to_box_intersection));
        break;
      }
    }
    //In case dist_to_boundingbox of this child is larger than any before it
    //has not been inserted above and is thus now appended to the end
    if (test_priority_list.size() == i) {
      test_priority_list.push_back(
          BVHBoxHelper<space_type,color_type>(
              bvh_node.child_ptrs[i],
              dist_to_box_intersection));
    }
  }
  
  bool hit_in_lvl_below = false;
  for(size_t i=0; i<amount_of_child_bvh; i++) {
    if (test_priority_list[i].dist_to_bounding_box !=
        std::numeric_limits<space_type>::max()) {
      hit_in_lvl_below = TestBVHLvl(
        ray,
        *(test_priority_list[i].bvh_node),
        distance_of_closest_hit_untill_now,
        results);
    } else { //At the point in the list where we dont have hits with the box we can stop
      break;
    }
    
    //If we hit a primitive in the currently check box we dont have to check
    //boxes that are further away
    if (hit_in_lvl_below) {
      break;
    }
  }

  return hit_in_lvl_below;
}

template<Float_t space_type,
         Float_t color_type>
bool IntersectBVH(
      const Ray<space_type>& ray,
      const BVH<space_type, color_type>& bvh,
      IntersectionResults<space_type,color_type>& results){
  //These will be set by the Intersect..() functions in case they hit
  space_type distance_of_closest_hit = std::numeric_limits<space_type>::max();
  
  //Test against everything in top-lvl bvh
  TestBVHLvl(ray,
             bvh.nodes[0],
             distance_of_closest_hit,
             results);
  
  if (distance_of_closest_hit < std::numeric_limits<space_type>::max()) {
    return true;
  } else {
    return false;
  }
}
