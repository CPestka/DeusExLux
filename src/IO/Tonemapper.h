#pragma once

//Maps floating point represenatation of an image into the integer representation

#include <cmath>
#include <cstdint>

#include "Image.h"
#include "UtilConcepts.h"

//BT.709 luma
template<Scalar_t color_type>
color_type ComputeLuminance(const Pixel<color_type>& pixel_in) {
  return (0.2126 * pixel_in.channel.x +
          0.7152 * pixel_in.channel.y +
          0.0722 * pixel_in.channel.z);
}

template<Float_t continuous_color_type>
continuous_color_type FindWhitePoint(
    const Image<Pixel<continuous_color_type>>& input_image) {
  continuous_color_type white_point = 0.0;

  for(int32_t i=0; i<input_image.width*input_image.height; i++){
    continuous_color_type tmp = ComputeLuminance(input_image.pixel_ptr[i]);
    white_point = tmp > white_point ? tmp : white_point;
  }

  return white_point;
}

template<Float_t continuous_color_type, Color_discrete_t discrete_color_type>
Image<Pixel<discrete_color_type>> PerformGlobalReinhardTonemapping(
    const Image<Pixel<continuous_color_type>>& input_image,
    std::optional<continuous_color_type> clamp_from = 1.0) {
  Image<Pixel<discrete_color_type>> output_image(input_image.width,
                                                 input_image.height,
                                                 Pixel<discrete_color_type>(0,0,0));

  continuous_color_type white_point =
      clamp_from ? clamp_from.value() : FindWhitePoint(input_image);

  continuous_color_type scale_factor =
      std::numeric_limits<discrete_color_type>::max() / white_point;
  
  for(int32_t i=0; i<input_image.width*input_image.height; i++){
    continuous_color_type old_luminance = ComputeLuminance(input_image.pixel_ptr[i]);

    continuous_color_type new_luminance =
        (old_luminance * (1.0 + (old_luminance / (white_point * white_point)))) /
        (1.0 + old_luminance);
    
    output_image.pixel_ptr[i] =
        (input_image.pixel_ptr[i] * (new_luminance / old_luminance) * scale_factor);
  }

  return output_image;
}
