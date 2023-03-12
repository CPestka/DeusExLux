#pragma once
//Contains the Boundingbox-Volume-Hirearchy i.e. BVH class, an acceleration data
//structure that harnesses the power of deviding and conquering to speed up the
//primitive intersection task that dominates the runtime of the renderer
//Also has some naive functions to construct BVHs from sceneces

#include <vector>
#include <memory>
#include <cmath>

#include "UtilConcepts.h"
#include "Primitives.h"
#include "Mesh.h"
#include "Timer.h"

template<Float_t space_type, Float_t color_type>
class BVHNode {
public:
  std::unique_ptr<Box<space_type,color_type>> bounding_box = nullptr;
  //Triangles are required to have a ascociated material but it will be ignored
  //here
  std::shared_ptr<PBRMaterial<space_type,color_type>> box_material = nullptr;

  std::unique_ptr<Tri<space_type,color_type>[]> tris = nullptr;
  uint32_t tri_count = 0;
  std::unique_ptr<Sphere<space_type,color_type>[]> spheres = nullptr;
  uint32_t sphere_count = 0;

  BVHNode<space_type,color_type>* parent = nullptr;
  std::unique_ptr<BVHNode<space_type,color_type>*[]> child_ptrs = nullptr;
  uint16_t child_count = 0;
  uint16_t child_reserved = 0;

  bool IsLeafNode() const {
    return child_ptrs == nullptr ? true : false;
  }
  bool HasNoPrimitives() const {
    return (tri_count + sphere_count) == 0;
  }
  bool HasParent() const {
    return parent != nullptr;
  }

  void AddChildPtr(BVHNode<space_type,color_type>* child_ptr) {
    if (child_count < child_reserved) {
      child_ptrs[child_count] = child_ptr;
      child_count++;
    } else {
      auto tmp = std::make_unique<BVHNode<space_type,color_type>*[]>(child_reserved + 10);
      child_reserved += 10;
      tmp[child_count] = child_ptr;
      child_count++;
      child_ptrs.swap(tmp);
    }
  }

  BVHNode() {};
  BVHNode(BVHNode<space_type,color_type>&& other) {
    bounding_box = other.bounding_box == nullptr ? 
        nullptr : std::move(other.bounding_box);
    box_material = other.box_material;

    tris = other.tris == nullptr ? 
        nullptr : std::move(other.tris);
    tri_count = other.tri_count;
    spheres = other.spheres == nullptr ? 
        nullptr : std::move(other.spheres);
    sphere_count = other.sphere_count;
    parent = other.parent;
    child_ptrs = other.child_ptrs == nullptr ?
        nullptr : std::move(other.child_ptrs);
    child_count = other.child_count;
  }
  BVHNode(
      Box<space_type,color_type> box,
      std::shared_ptr<PBRMaterial<space_type,color_type>> box_material,
      const std::vector<Tri<space_type,color_type>>& tris_,
      const std::vector<Sphere<space_type,color_type>>& spheres_,
      BVHNode<space_type,color_type>* parent_) :
        box_material(box_material),
        tri_count(tris_.size()),
        sphere_count(spheres_.size()),
        parent(parent_) {
    bounding_box = std::make_unique<Box<space_type,color_type>>();
    *(bounding_box.get()) = box;

    tris = std::make_unique<Tri<space_type,color_type>[]>(tri_count);
    spheres = std::make_unique<Sphere<space_type,color_type>[]>(sphere_count);

    for(uint32_t i=0; i<tri_count; i++){
      tris[i] = tris_[i];
    }
    for(uint32_t i=0; i<sphere_count; i++){
      spheres[i] = spheres_[i];
    }
  }
  BVHNode(
      Box<space_type,color_type> box,
      std::shared_ptr<PBRMaterial<space_type,color_type>> box_material,
      std::optional<BVHNode<space_type,color_type>*> parent_) :
        box_material(box_material){
    parent = parent_ ? parent_.value() : nullptr;
    bounding_box = std::make_unique<Box<space_type,color_type>>();
    *(bounding_box.get()) = box;
  }
  BVHNode<space_type,color_type>& operator=(BVHNode&& other) {
    bounding_box = other.bounding_box == nullptr ? 
        nullptr : std::move(other.bounding_box);
    box_material = other.box_material;

    tris = other.tris == nullptr ? 
        nullptr : std::move(other.tris);
    tri_count = other.tri_count;
    spheres = other.spheres == nullptr ? 
        nullptr : std::move(other.spheres);
    sphere_count = other.sphere_count;
    parent = other.parent;
    child_ptrs = other.child_ptrs == nullptr ?
        nullptr : std::move(other.child_ptrs);
    child_count = other.child_count;

    return *this;
  }

  void PrintNode() const {
    std::cout << "Tris: " << tri_count
              <<"\tSpheres: " << sphere_count
              << "\tChilds: " << child_count << "\n\n";
    return;
  }
  void PrintNodeAndChilds(int lvl) const {
    std::cout << "Current lvl: " << lvl << "\n" << std::endl;
    lvl++;

    PrintNode();
 
    for(uint16_t i=0; i<child_count; i++){
      child_ptrs[i]->PrintNodeAndChilds(lvl);
    }
  }
};

template<Float_t space_type, Float_t color_type>
class BVH {
  public:
  std::unique_ptr<BVHNode<space_type,color_type>[]> nodes = nullptr;
  uint32_t size = 0;
  uint32_t reserved = 0;

  uint32_t tri_count = 0;
  uint32_t sphere_count = 0;

  //returns a ptr to the added node
  BVHNode<space_type,color_type>* AddNode(BVHNode<space_type,color_type>&& node) {
    if (size < reserved) {
      nodes[size] = std::move(node);
      size++;
    } else {
      if (nodes == nullptr) {
        nodes = std::make_unique<BVHNode<space_type,color_type>[]>(10);
        reserved = 10;

        nodes[0] = std::move(node);
        size = 1;
      } else {
        auto tmp = std::make_unique<BVHNode<space_type,color_type>[]>(reserved + 10);
        reserved += 10;

        for(uint32_t i=0; i<size; i++){
          tmp[i] = std::move(nodes[i]);
        }
        tmp[size] = std::move(node);
        size++;
        nodes.swap(tmp);
      }
    }
    return &(nodes[size-1]);
  }

  void PrintBVH() const {
    if (size == 0) {
      return;
    }
    
    std::cout << "\nBVH: \n" << std::endl;
    std::cout << "Node count: " << size 
              << "\tTri count: " << tri_count
              << "\tSphere count: " << sphere_count << std::endl;
    nodes[0].PrintNodeAndChilds(0);
  }
};

template<Float_t space_type, Float_t color_type>
Box<space_type,color_type> GetFittingBox(
    const std::vector<Tri<space_type,color_type>>& tris,
    const std::vector<Sphere<space_type,color_type>>& spheres,
    std::shared_ptr<PBRMaterial<space_type,color_type>> box_material) {
  space_type max_x = std::numeric_limits<space_type>::min();
  space_type max_y = std::numeric_limits<space_type>::min();
  space_type max_z = std::numeric_limits<space_type>::min();
  space_type min_x = std::numeric_limits<space_type>::max();
  space_type min_y = std::numeric_limits<space_type>::max();
  space_type min_z = std::numeric_limits<space_type>::max();

  for(int64_t i=0; i<static_cast<int64_t>(tris.size()); i++){
    for(int j=0; j<3; j++){
    if (tris[i].vertex[j].x > max_x)
        max_x = tris[i].vertex[j].x;
    if (tris[i].vertex[j].x < min_x)
        min_x = tris[i].vertex[j].x;
    if (tris[i].vertex[j].y > max_y)
        max_y = tris[i].vertex[j].y;
    if (tris[i].vertex[j].y < min_y)
        min_y = tris[i].vertex[j].y;
    if (tris[i].vertex[j].z > max_z)
        max_z = tris[i].vertex[j].z;
    if (tris[i].vertex[j].z < min_z)
        min_z = tris[i].vertex[j].z;
    }
  }
  for(int64_t i=0; i<static_cast<int64_t>(spheres.size()); i++){
    if (spheres[i].origin.x + spheres[i].r > max_x)
        max_x = spheres[i].origin.x + spheres[i].r;
    if (spheres[i].origin.x - spheres[i].r < min_x)
        min_x = spheres[i].origin.x - spheres[i].r;
    if (spheres[i].origin.y + spheres[i].r > max_y)
        max_y = spheres[i].origin.y + spheres[i].r;
    if (spheres[i].origin.y - spheres[i].r < min_y)
        min_y = spheres[i].origin.y - spheres[i].r;
    if (spheres[i].origin.z + spheres[i].r > max_z)
        max_z = spheres[i].origin.z + spheres[i].r;
    if (spheres[i].origin.z - spheres[i].r < min_z)
        min_z = spheres[i].origin.z - spheres[i].r;
  }

  return Box<space_type,color_type>(Vec3<space_type>(min_x,min_y,min_z),
                                    std::fabs(max_x-min_x),
                                    std::fabs(max_y-min_y),
                                    std::fabs(max_z-min_z),
                                    box_material);
}

//quick and dirty check (treats sphere like box :))
template<Float_t space_type, Float_t color_type>
bool IsInBox(const Vec3<space_type>& my_vertex,
             const Box<space_type,color_type>& my_box){
  space_type max_x = my_box.tris[4].vertex[0].x;
  space_type max_y = my_box.tris[0].vertex[0].y;
  space_type max_z = my_box.tris[6].vertex[0].z;
  space_type min_x = my_box.tris[10].vertex[0].x;
  space_type min_y = my_box.tris[10].vertex[0].y;
  space_type min_z = my_box.tris[10].vertex[0].z;
  
  return (my_vertex.x <= max_x &&
          my_vertex.x >= min_x &&
          my_vertex.y <= max_y &&
          my_vertex.y >= min_y &&
          my_vertex.z <= max_z &&
          my_vertex.z >= min_z);
}

//quick and dirty check (treats sphere like box :))
template<Float_t space_type, Float_t color_type>
bool IsInBox(const Sphere<space_type,color_type>& my_sphere,
             const Box<space_type,color_type>& my_box){
  space_type max_x = my_box.tris[4].vertex[0].x;
  space_type max_y = my_box.tris[0].vertex[0].y;
  space_type max_z = my_box.tris[6].vertex[0].z;
  space_type min_x = my_box.tris[10].vertex[0].x;
  space_type min_y = my_box.tris[10].vertex[0].y;
  space_type min_z = my_box.tris[10].vertex[0].z;
  
  return (my_sphere.origin.x - my_sphere.r <= max_x &&
          my_sphere.origin.x + my_sphere.r >= min_x &&
          my_sphere.origin.y - my_sphere.r <= max_y &&
          my_sphere.origin.y + my_sphere.r >= min_y &&
          my_sphere.origin.z - my_sphere.r <= max_z &&
          my_sphere.origin.z + my_sphere.r >= min_z);
}

template<Float_t space_type, Float_t color_type>
bool IsInBox(const Tri<space_type,color_type>& tri,
             const Box<space_type,color_type>& box){
  return IsInBox(tri.vertex[0], box) || 
         IsInBox(tri.vertex[1], box) ||
         IsInBox(tri.vertex[2], box);
}



template<Float_t space_type, Float_t color_type>
std::vector<Box<space_type,color_type>> SplitBoxIn8Evenly(
    const Box<space_type,color_type>& parent,
    std::shared_ptr<PBRMaterial<space_type,color_type>> box_mat) {
  std::vector<Box<space_type,color_type>> child_boxes;
  Vec3<space_type> origin = parent.tris[10].vertex[0];
  //sizes vector is 0.5 * Vec3<space_type>(width,height,depth) of box
  Vec3<space_type> sizes =
      0.5 * (Vec3<space_type>(parent.tris[9].vertex[2]) - origin);

  child_boxes.push_back(Box<space_type,color_type>(
      origin,
      sizes.x,
      sizes.y,
      sizes.z,
      box_mat));
  child_boxes.push_back(Box<space_type,color_type>(
      origin + Vec3<space_type>(sizes.x,0.0,0.0),
      sizes.x,
      sizes.y,
      sizes.z,
      box_mat));
  child_boxes.push_back(Box<space_type,color_type>(
      origin + Vec3<space_type>(0.0,sizes.y,0.0),
      sizes.x,
      sizes.y,
      sizes.z,
      box_mat));
  child_boxes.push_back(Box<space_type,color_type>(
      origin + Vec3<space_type>(sizes.x,sizes.y,0.0),
      sizes.x,
      sizes.y,
      sizes.z,
      box_mat));
  child_boxes.push_back(Box<space_type,color_type>(
      origin + Vec3<space_type>(0.0,0.0,sizes.z),
      sizes.x,
      sizes.y,
      sizes.z,
      box_mat));
  child_boxes.push_back(Box<space_type,color_type>(
      origin + Vec3<space_type>(sizes.x,0.0,sizes.z),
      sizes.x,
      sizes.y,
      sizes.z,
      box_mat));
  child_boxes.push_back(Box<space_type,color_type>(
      origin + Vec3<space_type>(0.0,sizes.y,sizes.z),
      sizes.x,
      sizes.y,
      sizes.z,
      box_mat));
  child_boxes.push_back(Box<space_type,color_type>(
      origin + Vec3<space_type>(sizes.x,sizes.y,sizes.z),
      sizes.x,
      sizes.y,
      sizes.z,
      box_mat));
  
  return child_boxes;
}

template<Float_t space_type, Float_t color_type>
void SubdivideBVH(
    BVH<space_type,color_type>& bvh,
    Box<space_type,color_type> box,
    std::vector<Tri<space_type,color_type>>&& tris,
    std::vector<Sphere<space_type,color_type>>&& spheres,
    std::shared_ptr<PBRMaterial<space_type,color_type>> box_material,
    std::optional<BVHNode<space_type,color_type>*> parent_of_parent_ptr,
    uint16_t current_bvh_depth,
    uint16_t max_bvh_depth) {
  // std::cout << "current lvl " << current_bvh_depth << std::endl;
  // std::cout << "tris to distribute: " << tris.size()
  //           << "\nspheres to distribute: " << spheres.size() << std::endl;

  //Add parent node based on passed box
  std::optional<BVHNode<space_type,color_type>*> parent_ptr = 
      bvh.AddNode(BVHNode(
          box,
          box_material,
          parent_of_parent_ptr));

  auto child_boxes = SplitBoxIn8Evenly(box, box_material);

  for(int i=0; i<8; i++){
    std::vector<Tri<space_type,color_type>> tris_in_child;
    std::vector<Sphere<space_type,color_type>> spheres_in_child;

    for(size_t j=0; j<tris.size(); j++){
      //If tri is in box put it in tris_in_child and remove from tris
      if (IsInBox(tris[j], child_boxes[i])) { 
        tris_in_child.push_back(tris[j]);
        tris.erase(std::next(tris.begin(),j));
        j--;
      }
    }
    for(size_t j=0; j<spheres.size(); j++){
      if (IsInBox(spheres[j], child_boxes[i])) { 
        spheres_in_child.push_back(spheres[j]);
        spheres.erase(std::next(spheres.begin(),j));
        j--;
      }
    }
    // std::cout << "child box id: " << i
    //           << "\ttris: " << tris_in_child.size()
    //           << "\tsph: " << spheres_in_child.size() << std::endl;
    
    uint32_t prim_in_child = (tris_in_child.size() + spheres_in_child.size());

    if (prim_in_child != 0) {
      auto current_box =
          GetFittingBox(tris_in_child,spheres_in_child,box_material);

      if (prim_in_child < 12 || current_bvh_depth == max_bvh_depth) {
        //Add leave node to bvh if only a small number of primitives are in the
        //child or if max recursion depth has been reached
        //std::cout << "adding child" << std::endl;
        auto ptr_to_this_node = bvh.AddNode(BVHNode(
            current_box,
            box_material,
            tris_in_child,
            spheres_in_child,
            parent_ptr.value()));

        (parent_ptr.value())->AddChildPtr(ptr_to_this_node);
        
      } else {
        SubdivideBVH(
            bvh,
            current_box,
            std::move(tris_in_child),
            std::move(spheres_in_child),
            box_material,
            parent_ptr,
            current_bvh_depth + 1,
            max_bvh_depth);
      }
    }
    tris_in_child.clear();
    spheres_in_child.clear();
  }
}

//Used to construct the BVH for a scene
//Returns unique ptr array that holds all "BVH-nodes" of the entire BVH
//Primitives that are passed in will be owned by the leaf BVH nodes that they
//are contained in
template<Float_t space_type, Float_t color_type>
BVH<space_type,color_type> ConstructOptimizedBVH(
    std::vector<Tri<space_type,color_type>>&& tris,
    std::vector<Sphere<space_type,color_type>>&& spheres,
    std::shared_ptr<PBRMaterial<space_type,color_type>> box_material,
    uint16_t max_bvh_depth){
  std::cout << "Constructing Bounding Box Volume Hirearchy (BVH) acceleration structure\n\n"
            << "From Scene with:\n"
            << "Triangles:\t" << tris.size() << "\n"
            << "Spheres:\t" << spheres.size() << "\n\n"
            << "Maximum BVH depth:\t" << max_bvh_depth << "\n" << std::endl;

  IntervallTimer timer;

  Box<space_type,color_type> scene_box = GetFittingBox(
      tris,
      spheres,
      box_material);
  
  BVH<space_type,color_type> bvh;
  bvh.tri_count = tris.size();
  bvh.sphere_count = spheres.size();

  SubdivideBVH<space_type,color_type>(
      bvh,
      scene_box,
      std::move(tris),
      std::move(spheres),
      box_material,
      std::nullopt,
      1,
      max_bvh_depth);
  
  double time = timer.getTimeInMicroseconds() * 0.000001;
  
  std::cout << "BVH successfully constructed.\n"
            << "Time elapsed:\t" << time << "s\n"
            << "Total BVH nodes:\t" << bvh.size << "\n\n" << std::endl;

  return bvh;
}