#pragma once

#include <cmath>

#include "UtilConcepts.h"
#include "Vector.h"

//For this entire project the default coordinate systems are cartesian if not
//state otherwise!
//
//All coord sys are of course "Rechtssysteme" i.e. right hand rule aplies
//
//Convention is that if x or x0 = "view-dir" -> y or x1 = "up"

//A Cartesian coordinate system relative to the "base cartesian coord-sys" 
//ONB + origin + util to transform it
template<Float_t space_type>
class Coordinates {
public:
  Vec3<space_type> origin = Vec3<space_type>(0.0,0.0,0.0);
  Vec3<space_type> e1_view_dir = Vec3<space_type>(1.0,0.0,0.0);
  Vec3<space_type> e2_up = Vec3<space_type>(0.0,1.0,0.0);
  Vec3<space_type> e3_right = Vec3<space_type>(0.0,0.0,1.0);
  space_type pitch = 0.0;
  space_type yaw = 0.0;
  space_type roll = 0.0;
  
  void Translate(const Vec3<space_type>& offset) {
    origin += offset;
  }
  void SetOrigin(const Vec3<space_type>& new_origin) {
    origin = new_origin;
  }

  void Pitch(space_type phi) {
    e1_view_dir.RotateAroundVec(e3_right, phi);
    e2_up.RotateAroundVec(e3_right, phi);
    pitch = GetAngleBetween(Vec3<space_type>(1.0,0.0,0.0), e1_view_dir);
  }
  void Yaw(space_type phi) {
    e1_view_dir.RotateAroundVec(e2_up, phi);
    e3_right.RotateAroundVec(e2_up, phi);
    yaw = GetAngleBetween(Vec3<space_type>(0.0,0.0,1.0), e3_right);
  }
  void Roll(space_type phi) {
    e2_up.RotateAroundVec(e1_view_dir, phi);
    e3_right.RotateAroundVec(e1_view_dir, phi);
    roll = GetAngleBetween(Vec3<space_type>(0.0,1.0,0.0), e2_up);
  }
};

