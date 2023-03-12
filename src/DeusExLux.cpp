#include <iostream>

#include "UtilConcepts.h"
#include "Camera.h"

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

  return EXIT_SUCCESS;
}