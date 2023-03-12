#include <iostream>

#include "UtilConcepts.h"
#include "Camera.h"
#include "Material.h"
#include "Primitives.h"
#include "BVH.h"
#include "RenderImage.h"
#include "Tonemapper.h"
#include "PPMFileWriter.h"

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
  typedef uint8_t descrete_color_t;

  constexpr int32_t img_width = 192;
  constexpr int32_t img_height = 108;

  constexpr space_t f_stop = 11;
  constexpr space_t focal_length = 20;
  constexpr space_t focus_distance = 1.0;
  constexpr space_t sensor_size = 35; //35 == Full frame sensor

  constexpr uint16_t max_bvh_depth = 10; 

  uint16_t render_thread_count = std::thread::hardware_concurrency();
  constexpr uint32_t tile_size_x = 64;
  constexpr uint32_t tile_size_y = 64;

  constexpr uint32_t samples_per_pixel = 20;
  constexpr uint32_t max_ray_bounces = 1;

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

  //Construct BVH that will hold primitives and is used for the raytracing
  auto bvh = ConstructOptimizedBVH<space_t,color_t>(
      std::move(tris),
      std::move(spheres),
      basic_material,
      max_bvh_depth);

  //RenderImage or ContinueRenderImage are used to Render the scene via
  //raytracing
  auto render = RenderImage<space_t,color_t,RasterRay>(
      cam,
      samples_per_pixel,
      max_ray_bounces,
      tile_size_x,
      tile_size_y,
      bvh,
      background,
      render_thread_count);

  auto descrete_img =
      PerformGlobalReinhardTonemapping<color_t,descrete_color_t>(render);

  if (!WriteImageToPPMFile<descrete_color_t>(descrete_img, "./out.ppm")) {
    std::cout << "Failed to write image to file" << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}