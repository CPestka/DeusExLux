#pragma once
//Includes some composites of primitives for ease of use when constructing a 
//scene and for the bvh

#include "UtilConcepts.h"
#include "Primitives.h"

template<Float_t space_type, Float_t color_type>
class Box {
public:
  Tri<space_type, color_type> tris[12];
  
  Box() {};
  //with no offset provided by origin the box is in the positive quadrant
  Box(const Vec3<space_type>& origin, 
      space_type width,
      space_type height,
      space_type depth,
      std::shared_ptr<PBRMaterial<space_type,color_type>> mat) {
          
    tris[0] = 
        Tri<space_type,color_type>(origin + Vec3<space_type>(0.0,height,0.0),     //top face
        origin + Vec3<space_type>(width,height,0.0),
        origin + Vec3<space_type>(0.0,height,depth),
        mat);
    
    tris[1] = 
        Tri<space_type,color_type>(origin + Vec3<space_type>(0.0,height,depth),
        origin + Vec3<space_type>(width,height,0.0),
        origin + Vec3<space_type>(width,height,depth),
        mat);

    tris[2] = 
        Tri<space_type,color_type>(origin + Vec3<space_type>(0.0,height,0.0),    //left face
        origin + Vec3<space_type>(0.0,0.0,depth),
        origin + Vec3<space_type>(0.0,height,depth),
        mat);

    tris[3] = 
        Tri<space_type,color_type>(origin + Vec3<space_type>(0.0,height,0.0),
        origin + Vec3<space_type>(0.0,0.0,depth),
        origin + Vec3<space_type>(0.0,0.0,0.0),
        mat);

    tris[4] = 
        Tri<space_type,color_type>(origin + Vec3<space_type>(width,height,0.0),   //right face
        origin + Vec3<space_type>(width,height,depth),
        origin + Vec3<space_type>(width,0.0,depth),
        mat);

    tris[5] = 
        Tri<space_type,color_type>(origin + Vec3<space_type>(width,height,0.0),
        origin + Vec3<space_type>(width,0.0,0.0),
        origin + Vec3<space_type>(width,0.0,depth),
        mat);

    tris[6] = 
        Tri<space_type,color_type>(origin + Vec3<space_type>(0.0,0.0,depth),      //front face
        origin + Vec3<space_type>(0.0,height,depth),
        origin + Vec3<space_type>(width,0.0,depth),
        mat);

    tris[7] =
        Tri<space_type,color_type>(origin + Vec3<space_type>(0.0,height,depth),
        origin + Vec3<space_type>(width,height,depth),
        origin + Vec3<space_type>(width,0.0,depth),
        mat);

    tris[8] =
        Tri<space_type,color_type>(origin + Vec3<space_type>(0.0,0.0,depth),     //back face
        origin + Vec3<space_type>(width,0.0,depth),
        origin + Vec3<space_type>(0.0,height,depth),
        mat);

    tris[9] =
        Tri<space_type,color_type>(origin + Vec3<space_type>(0.0,height,depth),
        origin + Vec3<space_type>(width,0.0,depth),
        origin + Vec3<space_type>(width,height,depth),
        mat);

    tris[10] = 
        Tri<space_type,color_type>(origin + Vec3<space_type>(0.0,0.0,0.0),       //bottom face
        origin + Vec3<space_type>(0.0,0.0,depth),
        origin + Vec3<space_type>(width,0.0,0.0),
        mat);

    tris[11] =
        Tri<space_type,color_type>(origin + Vec3<space_type>(0.0,0.0,depth),
        origin + Vec3<space_type>(width,0.0,depth),
        origin + Vec3<space_type>(width,0.0,0.0),
        mat);

  }
  
  void RotateAzimuth(space_type phi) {
    for(int i=0; i<12; i++){
      tris[i].RotateAzimuth(phi);
    } 
  }
  void RotatePolar(space_type theta) {
    for(int i=0; i<12; i++){
      tris[i].RotatePolar(theta);
    } 
  }
  void Rotate(space_type theta, space_type phi) {
    for(int i=0; i<12; i++){
      tris[i].Rotate(theta, phi);
    }
  }

  void Translate(const Vec3<space_type>& offset_vector) {
    for(int i=0; i<12; i++){
      tris[i].Translate(offset_vector);
    }
  }

  void ScaleWidth(space_type scale) {
    for(int i=0; i<12; i++){
      tris[i].ScaleWidth(scale);
    }
  }
  void ScaleHeight(space_type scale) {
    for(int i=0; i<12; i++){
      tris[i].ScaleHeight(scale);
    }
  }
  void ScaleDepht(space_type scale) {
    for(int i=0; i<12; i++){
      tris[i].ScaleDepht(scale);
    }
  }
  void ScaleUniformly(space_type scale) {
    for(int i=0; i<12; i++){
      tris[i].ScaleUniformly(scale);
    }
  }
};

