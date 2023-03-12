#pragma once

#include <cmath>
#include <optional>
#include <mutex>
#include <thread>

#include "UtilConcepts.h"
#include "Vector.h"
#include "Camera.h"
#include "Image.h"
#include "Material.h"
#include "BVH.h"
#include "RNG.h"
#include "Primitives.h"
#include "Integrator.h"
#include "Timer.h"

enum RenderMode {
  RasterRay, //Most akin to regular rasterization = only samples albedo on hit, no bounces
  AlbedoRay, //Proper bounces used, but only considers albedo over full PBR properties
  FullRayNoTransparency, //Considers all PBR properties execpt transperency
  FullRay
};

template<Float_t space_type, Float_t color_type, RenderMode mode>
Pixel<color_type> RenderPixel(
    const Camera<space_type,color_type>& camera,
    uint32_t x,
    uint32_t y,
    uint32_t samples_per_pixel,
    [[maybe_unused]] uint32_t max_bounces, //unused in RasterRay path
    const BVH<space_type,color_type>& bvh,
    RNGLutEngine<space_type>& rng_engine,
    const std::optional<Sphere<space_type,color_type>>& background) {
  Pixel<color_type> pxl(0.0,0.0,0.0);
  for(uint32_t i=0; i<samples_per_pixel; i++){
    //Free AA yeahhhh
    space_type inter_pixel_delta_x = rng_engine.GetUniformRandom();
    space_type inter_pixel_delta_y = rng_engine.GetUniformRandom();

    Vec3<space_type> viewport_intersection_point = 
        ((x + inter_pixel_delta_x) * camera.pixel_size *
         camera.viewport_coordinates.e3_right) +
        ((y + inter_pixel_delta_y) * camera.pixel_size * (-1) *
         camera.viewport_coordinates.e2_up) +
        camera.viewport_top_left;
      
    //Depth of field :)
    space_type r = rng_engine.GetUniformRandom(0.0, camera.lense_radius);
    space_type phi = rng_engine.GetUniformRandom(0.0, 2 * M_PI);
    space_type lense_delta_x = r * std::cos(phi);
    space_type lense_delta_y = r * std::sin(phi);

    Vec3<space_type> camera_lense_exit_point = 
        (lense_delta_x * camera.camera_coordinates.e3_right) +
        (lense_delta_y * camera.camera_coordinates.e2_up) +
        camera.camera_coordinates.origin;
    
    //Select required raytracing function based on selected mode at COMPILETIME
    if constexpr (mode == RasterRay) {
      pxl += TraceRasterRay<space_type,color_type>(
        Ray<space_type>(camera_lense_exit_point, viewport_intersection_point),
        bvh,
        background);
    } else {
      std::terminate();
    }
    
  }
  return pxl * (1.0 / static_cast<space_type>(samples_per_pixel));
}

template<Float_t space_type, Float_t color_type, RenderMode mode>
void RenderTile(
      Image<Pixel<color_type>>& render,
      const Camera<space_type,color_type>& camera,
      uint32_t x,
      uint32_t x_max,
      uint32_t y,
      uint32_t y_max,
      uint32_t samples_per_pixel,
      uint32_t max_bounces,
      const BVH<space_type,color_type>& bvh,
      RNGLutEngine<space_type>& rng_engine,
      const std::optional<Sphere<space_type,color_type>>& background) {
    for(uint32_t i=y; i<=y_max; i++){
      for(uint32_t j=x; j<=x_max; j++){
        render.pixel_ptr[j + (camera.image_width * i)] +=
            RenderPixel<space_type,color_type,mode>(
                camera,j,i,samples_per_pixel,max_bounces,bvh,rng_engine,background);
      }
    }
    //std::cout << "Finished a Tile :)\n";
    return;
  }

struct TileInfo {
  uint32_t x;
  uint32_t x_max;
  uint32_t y;
  uint32_t y_max;
  TileInfo(uint32_t x, uint32_t x_max, uint32_t y, uint32_t y_max) :
      x(x), x_max(x_max), y(y), y_max(y_max) {};
};

class TileQueue {
  private:
  std::vector<TileInfo> remaining_tiles;
  std::mutex mtx;

  public:
  void InsertTiles(const std::vector<TileInfo>& todo_tiles) {
    std::lock_guard<std::mutex> lck_grd(mtx);

    for(size_t i=0; i<todo_tiles.size(); i++){
      remaining_tiles.push_back(todo_tiles[i]);
    }
  }

  std::optional<TileInfo> AcquireTile() {
    std::lock_guard<std::mutex> lck_grd(mtx);

    if (remaining_tiles.size() == 0) {
      return std::nullopt;
    }
    TileInfo tmp = remaining_tiles.back();
    remaining_tiles.pop_back();
    return tmp;
  }

  size_t GetTileCount() {
    std::lock_guard<std::mutex> lck_grd(mtx);
    
    return remaining_tiles.size();
  }
};

template<Float_t space_type, Float_t color_type, RenderMode mode>
void RenderTiles(
    Image<Pixel<color_type>>& render,
    const Camera<space_type,color_type>& camera,
    uint32_t samples_per_pixel,
    uint32_t max_bounces,
    const BVH<space_type,color_type>& bvh,
    const std::optional<Sphere<space_type,color_type>>& background,
    TileQueue& tile_queue) {
  RNGLutEngine<space_type> rng_engine{};

  while(true) {
    auto current_tile = tile_queue.AcquireTile();

    if (!current_tile) {
      return;
    }

    RenderTile<space_type,color_type,mode>(
               render,
               camera,
               current_tile.value().x,
               current_tile.value().x_max,
               current_tile.value().y,
               current_tile.value().y_max,
               samples_per_pixel,
               max_bounces,
               bvh,
               rng_engine,
               background);
  }
}

template<Float_t space_type, Float_t color_type, RenderMode mode>
void ContinueRenderingImage(
    Image<Pixel<color_type>>& render,
    const Camera<space_type,color_type>& camera,
    uint32_t samples_per_pixel,
    uint32_t max_bounces,
    uint32_t tile_size_x,
    uint32_t tile_size_y,
    const BVH<space_type,color_type>& bvh,
    const std::optional<Sphere<space_type,color_type>>& background,
    uint16_t thread_count){
  
  std::vector<TileInfo> tiles;
  for(uint32_t i=0; i<camera.image_height; i+=tile_size_y){
    for(uint32_t j=0; j<camera.image_width; j+=tile_size_x){
      uint32_t tile_width =
          j + tile_size_x >= camera.image_width ?
              camera.image_width - j : tile_size_x;
      uint32_t tile_height =
          i + tile_size_y >= camera.image_height ?
              camera.image_height - i : tile_size_y;
      tiles.push_back(TileInfo(j, j + tile_width - 1, i, i + tile_height - 1));
    }
  }

  TileQueue queue;
  queue.InsertTiles(tiles);

  std::vector<std::thread> render_threads;

  std::cout << "Starting to render: \n"
            << "Resolution:\t" << camera.image_width << " x " << camera.image_height << "\n"
            << "Total Pixels:\t" << camera.image_width * camera.image_height << "\n"
            << "Tile Size:\t" << tile_size_x << " x " << tile_size_y << "\n"
            << "Total Tiles:\t" << queue.GetTileCount() << "\n"
            << "Samples Per Pixel:\t" << samples_per_pixel << "\n"
            << "Maximum Bounces:\t" << max_bounces << "\n"
            << "Number of Tris in Scene:\t" << bvh.tri_count << "\n"
            << "Number of Spheres in Scene:\t" << bvh.sphere_count << "\n"
            << "Amount of Render threads:\t" << thread_count << std::endl;
  
  IntervallTimer timer;

  //Start the other threads
  for(int i=0; i<thread_count-1; i++){
    render_threads.push_back(std::thread(
      RenderTiles<space_type,color_type,mode>,
      std::ref(render),
      std::ref(camera),
      samples_per_pixel,
      max_bounces,
      std::ref(bvh),
      std::ref(background),
      std::ref(queue)));
  }
  
  //Main thread works itself
  RenderTiles<space_type,color_type,mode>(
    render,
    camera,
    samples_per_pixel,
    max_bounces,
    bvh,
    background,
    queue);

  for(int i=0; i<thread_count-1; i++){
    render_threads[i].join();
  }

  double time = timer.getTimeInMicroseconds() * 0.000001;

  std::cout << "Render completed!\n\n Time elapsed:\t" << time << "s" << std::endl;
}

template<Float_t space_type, Float_t color_type, RenderMode mode>
Image<Pixel<color_type>> RenderImage(
    const Camera<space_type,color_type>& camera,
    uint32_t samples_per_pixel,
    uint32_t max_bounces,
    uint32_t tile_size_x,
    uint32_t tile_size_y,
    const BVH<space_type,color_type>& bvh,
    const std::optional<Sphere<space_type,color_type>>& background,
    uint16_t thread_count) {
  auto render =
      Image<Pixel<color_type>>(camera.image_width, camera.image_height, Pixel<color_type>(0,0,0));
  
  ContinueRenderingImage<space_type,color_type,mode>(
      render,
      camera,
      samples_per_pixel,
      max_bounces,
      tile_size_x,
      tile_size_y,
      bvh,
      background,
      thread_count);
    
    return render;
  }