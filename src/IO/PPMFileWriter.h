#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <limits>

#include "UtilConcepts.h"
#include "Image.h"

template<Color_discrete_t color_channel_type>
bool WriteImageToPPMFile(Image<Pixel<color_channel_type>> image_in,
                         std::string file_path){
  std::ofstream file(file_path, std::ios_base::out | std::ios_base::binary);

  if (file.is_open()) {
    //Write PPM header
    file << "P6" << std::endl
         << image_in.width << " "
         << image_in.height << std::endl
         << std::to_string(std::numeric_limits<color_channel_type>::max())
         << std::endl;
    //Write pixels
    for(int i=0; i<image_in.height; i++){
      for(int j=0; j<image_in.width; j++){
          file << image_in.pixel_ptr[j + image_in.width * i].channel.x
               << image_in.pixel_ptr[j + image_in.width * i].channel.y
               << image_in.pixel_ptr[j + image_in.width * i].channel.z;
      }
    }
    return true;
  } else {
    std::cout << "Failed to write image to file: " << file_path << std::endl;
    return false;
  }
}
