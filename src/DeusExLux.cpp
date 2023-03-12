#include <iostream>

#include "UtilConcepts.h"
#include "Camera.h"
#include "Material.h"
#include "Primitives.h"

//TODO:
// check rotations; diffuse brdf no transmission; cooler scene; + transmission; specular + t
//blended + t
//      - basic scene
//      - scene import
//      - gpu
//      - volumetrics

int main(){
  typedef float space_t;
  typedef float color_t;

  constexpr int32_t img_width = 192;
  constexpr int32_t img_height = 108;

  constexpr space_t f_stop = 11;
  constexpr space_t focal_length = 20;
  constexpr space_t focus_distance = 1.0;
  constexpr space_t sensor_size = 35; //35 == Full frame sensor

  Camera<space_t,color_t> cam(
      img_width,
      img_height, 
      f_stop,
      focal_length,
      focus_distance,
      sensor_size);
 
  //Make example Scene
  //
  //Make primitves
  std::shared_ptr<PBRMaterial<space_t,color_t>> basic_material =
      std::make_shared<PBRMaterial<space_t,color_t>>();
  std::shared_ptr<PBRMaterial<space_t,color_t>> background_material =
      std::make_shared<PBRMaterial<space_t,color_t>>();
  background_material.get()->albedo_const = Pixel<color_t>(0.1,0.1,0.1);
  background_material.get()->emitance_const = Pixel<color_t>(0.05,0.0,0.0);

  Sphere<space_t,color_t> background(
      std::numeric_limits<space_t>::max(),
      0,0,0,
      background_material);

  std::vector<Sphere<space_t,color_t>> spheres;
  spheres.push_back(
      Sphere<space_t,color_t>(1.0, Vec3<float>(7.0,1.0,0.0), basic_material));
  spheres.push_back(
      Sphere<space_t,color_t>(0.5, Vec3<float>(7.0,2.0,0.0), basic_material));

  std::vector<Tri<space_t,color_t>> tris;
  tris.push_back(Tri<space_t,color_t>(Vec3<space_t>(10.0,0,-10.0),
                                      Vec3<space_t>(10.0,2.0,-2.0),
                                      Vec3<space_t>(10.0,2.0,2.0),
                                      basic_material));
  tris.push_back(Tri<space_t,color_t>(Vec3<space_t>(10.0,4.0,-10.0),
                                      Vec3<space_t>(10.0,12.0,2.0),
                                      Vec3<space_t>(10.0,12.0,-8.0),
                                      basic_material));


  return EXIT_SUCCESS;
}