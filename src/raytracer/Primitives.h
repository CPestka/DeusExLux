#pragma once
//Contains all primitives considered for the raytracing and the ray class
//These are atm triangles i.e. tris and spheres
//For containers/helper classes to deal with more complex geometry see
//Mesh.h

#include <memory>

#include "UtilConcepts.h"
#include "Vector.h"
#include "Material.h"

template<Float_t space_type, Float_t color_type>
class PBRMaterial;

template<Float_t space_type>
class Ray {
public:
  Vec3<space_type> origin;
  Vec3<space_type> direction;
  
  Ray() {};
  //UB for origin == direction
  Ray(Vec3<space_type> origin, Vec3<space_type> direction) :
      origin(origin), direction(direction) {
    direction.Normalize();
  }
};

template<Float_t space_type, Float_t color_type>
class Primitive {
public:
  std::shared_ptr<PBRMaterial<space_type, color_type>> material;
};

//Does not support r == 0
template<Float_t space_type, Float_t color_type>
class Sphere: public Primitive<space_type, color_type> {
public:
  Vec3<space_type> origin;
  space_type r;
  
  Sphere() {};
  Sphere(
      space_type r,
      const Vec3<space_type>& origin,
      std::shared_ptr<PBRMaterial<space_type,color_type>> mat) : 
          origin(origin), r(r) {
    this->material = mat;
  };
  Sphere(space_type r,
         space_type origin_x,
         space_type origin_y,
         space_type origin_z,
         std::shared_ptr<PBRMaterial<space_type,color_type>> mat) :
      origin(Vec3<space_type>(origin_x,origin_y,origin_z)), r(r) {
    this->material = mat;
  };
  
  //Returns surface normal if point is on the surface and otherwise the normal
  //of the surface point that is above/below the specified point
  Vec3<space_type> ComputeSurfaceNormal(const Vec3<space_type>& point) const {
    return ComputeNormal(Vec3<space_type>(point) - origin);
  }

  void Translate(const Vec3<space_type>& offset_vector) {
    origin += offset_vector;
  }
  void Scale(space_type scale) {
    r *= scale;
  }
};

template<Float_t space_type, Float_t color_type>
class Tri: public Primitive<space_type, color_type> {
public:
  Vec3<space_type> vertex[3] =
      {Vec3<space_type>(0.0,0.0,0.0),
       Vec3<space_type>(1.0,0.0,0.0),
       Vec3<space_type>(0.0,1.0,0.0)};
  Vec3<space_type> geometry_normal = Vec3<space_type>(0.0,0.0,1.0);
  
  Tri() {};
  Tri(std::shared_ptr<PBRMaterial<space_type,color_type>> mat) {
    this->material = mat;
  };
  Tri(const Vec3<space_type>& p0,
      const Vec3<space_type>& p1,
      const Vec3<space_type>& p2,
      std::shared_ptr<PBRMaterial<space_type,color_type>> mat) {
    this->material = mat;
    
    vertex[0] = p0;
    vertex[1] = p1;
    vertex[2] = p2;

    ComputeNormal();
  };
  
  void ComputeNormal(){
    geometry_normal =
        ComputeCrossProduct(vertex[1]-vertex[0],
                            vertex[2]-vertex[0]);
    geometry_normal.Normalize();
  }
  //Util func to test if triangle is actually valid 
  bool AreAnyVerteciesIdentical() {
    return (vertex[0] == vertex[1] ||
            vertex[0] == vertex[2] ||
            vertex[1] == vertex[2]);
  }
  
  //rotates the triangel so that the normal aligns with the y axis
  void RotateTriToNormalUp() {
    vertex[0].PerformAlignRotation(geometry_normal);
    vertex[1].PerformAlignRotation(geometry_normal);
    vertex[2].PerformAlignRotation(geometry_normal);
    geometry_normal = {0.0,1.0,0.0};
  }

  void RotatePolar(space_type theta) {
    vertex[0].RotateCartesianVecPolar(theta);
    vertex[1].RotateCartesianVecPolar(theta);
    vertex[2].RotateCartesianVecPolar(theta);
    ComputeNormal();
  }
  void RotateAzimuth(space_type phi) {
    vertex[0].RotateCartesianVecAzimuthal(phi);
    vertex[1].RotateCartesianVecAzimuthal(phi);
    vertex[2].RotateCartesianVecAzimuthal(phi);
    ComputeNormal();
  }
  void Rotate(space_type theta, space_type phi) {
    vertex[0].RotateCartisianVec(theta, phi);
    vertex[1].RotateCartisianVec(theta, phi);
    vertex[2].RotateCartisianVec(theta, phi);
    ComputeNormal();
  }

  void Translate(const Vec3<space_type>& offset_vector) {
    vertex[0] += offset_vector;
    vertex[1] += offset_vector;
    vertex[2] += offset_vector;
  }
  void ScaleWidth(space_type scale) {
    vertex[0].data[0] *= scale;
    vertex[1].data[0] *= scale;
    vertex[2].data[0] *= scale;
  }
  void ScaleHeight(space_type scale) {
    vertex[0].data[1] *= scale;
    vertex[1].data[1] *= scale;
    vertex[2].data[1] *= scale;
  }
  void ScaleDepth(space_type scale) {
    vertex[0].data[2] *= scale;
    vertex[1].data[2] *= scale;
    vertex[2].data[2] *= scale;
  }
  void ScaleUniformly(space_type scale) {
    vertex[0].data *= scale;
    vertex[1].data *= scale;
    vertex[2].data *= scale;
  }
};