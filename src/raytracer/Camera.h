#pragma once

#include <cmath>
#include <optional>

#include "UtilConcepts.h"
#include "Vector.h"
#include "AlgebraUtil.h"

template<Float_t space_type,
         Float_t color_type>
class Camera {
public:
  uint32_t image_width;
  uint32_t image_height;
  space_type aspect_ratio;

  space_type fov_horizontal;
  space_type fov_vertical;
  space_type viewport_half_width;
  space_type viewport_half_height;
  space_type pixel_size;
  RelativeCoordinates<space_type> viewport_coordinates;
  Vec3<space_type> viewport_top_left;

  RelativeCoordinates<space_type> camera_coordinates;
  
  space_type lense_radius;

  Camera(int32_t image_width = 1920,
         int32_t image_height = 1080,
         space_type f_stop = 2.8,
         space_type focal_length = 100, //in mm 
         space_type focus_distance = 1.0,
         space_type sensor_size = 35 //in mm
         ) :
      image_width(image_width),
      image_height(image_height),
      aspect_ratio(image_width/image_height),
      lense_radius((focal_length * 0.001) / f_stop) {
    space_type sensor_width = 0.001 * sensor_size;
    space_type sensor_height = sensor_width / aspect_ratio;

    fov_horizontal = 2 * std::atan((sensor_width) / (2*0.001*focal_length));
    fov_vertical = 2 * std::atan((sensor_height) / (2*0.001*focal_length));
    
    viewport_half_width = focus_distance * std::tan(fov_horizontal*0.5);
    viewport_half_height = viewport_half_width / aspect_ratio;

    pixel_size = 2 * (viewport_half_width / image_width);

    viewport_coordinates.origin = Vec3<space_type>(focus_distance, 0.0, 0.0);
    viewport_top_left =
        viewport_coordinates.origin +
        Vec3<space_type>(0.0,viewport_half_height,-viewport_half_width);
  }

  space_type GetFocusDistance() {
    return (viewport_coordinates.origin - camera_coordinates.origin).GetLength();
  }

  void MoveFocus(space_type offset) {
    if (offset + GetFocusDistance() <= 0.0) {
      return;
    }
    
    camera_coordinates.origin += (camera_coordinates.e1_view_dir * offset);
    viewport_half_width = GetFocusDistance() * std::tan(fov_horizontal*0.5);
    viewport_half_height = viewport_half_width / aspect_ratio;
    pixel_size = viewport_half_width / image_width;
  }

  void MoveCameraForward(space_type offset) {
    camera_coordinates.origin += (camera_coordinates.e1_view_dir * offset);
    viewport_coordinates.origin += (viewport_coordinates.e1_view_dir * offset);
  }
  
  //Util funcs to move camera that also deal with the required trafos for the
  //viewport
  void Jaw(space_type phi) {
    space_type current_focus_distance = GetFocusDistance();
    camera_coordinates.Yaw(phi);
    viewport_coordinates.Yaw(phi);
    viewport_coordinates.origin =
        camera_coordinates.origin +
        (viewport_coordinates.e1_view_dir * current_focus_distance);
    viewport_top_left =
        viewport_coordinates.origin + 
        (viewport_coordinates.e2_up * viewport_half_height) -
        (viewport_coordinates.e3_right * viewport_half_width);
  }

  void Pitch(space_type theta) {
    space_type current_focus_distance = GetFocusDistance();
    camera_coordinates.Pitch(theta);
    viewport_coordinates.Pitch(theta);
    viewport_coordinates.origin =
        camera_coordinates.origin +
        (viewport_coordinates.e1_view_dir * current_focus_distance);
    viewport_top_left =
        viewport_coordinates.origin + 
        (viewport_coordinates.e2_up * viewport_half_height) -
        (viewport_coordinates.e3_right * viewport_half_width);
  }

  void Roll(space_type theta) {
    space_type current_focus_distance = GetFocusDistance();
    camera_coordinates.Roll(theta);
    viewport_coordinates.Roll(theta);
    viewport_coordinates.origin =
        camera_coordinates.origin +
        (viewport_coordinates.e1_view_dir * current_focus_distance);
    viewport_top_left =
        viewport_coordinates.origin + 
        (viewport_coordinates.e2_up * viewport_half_height) -
        (viewport_coordinates.e3_right * viewport_half_width);
  }
};