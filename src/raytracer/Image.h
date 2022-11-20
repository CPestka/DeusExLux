#pragma once
//Contains the different classes of images and pixels used during rendering for
//the output

#include <cstdint>
#include <memory>
#include <cstring>

#include "UtilConcepts.h"
#include "Vector.h"

//E.g. T=Pixel<colortype> for image of intensities or rgb, or T=Vec3<space_type>
//for normal map
template<typename T>
class Image {
public:
  std::unique_ptr<T[]> pixel_ptr = nullptr;
  int32_t width = 0;
  int32_t height = 0;
  
  Image() {};
  Image(int32_t width, int32_t height, T init_value) : width(width), height(height) {
    pixel_ptr = std::make_unique<T[]>(width*height);
    for(int32_t i=0; i<width*height; i++){
      pixel_ptr[i] = init_value;
    }
  };
  Image(const Image<T>& other) {
    width = other.width;
    height = other.height;
    pixel_ptr = std::make_unique<T[]>(width * height);
    std::memcpy(
        pixel_ptr.get(),
        other.pixel_ptr.get(),
        width * height * sizeof(T));
  }
  Image(Image<T>&& other) {
    width = other.width;
    height = other.height;
    pixel_ptr = std::move(other.pixel_ptr);
  }

  Image<T>& operator=(const Image<T>& other) {
    width = other.width;
    height = other.height;
    pixel_ptr = std::make_unique<T[]>(width * height);
    std::memcpy(pixel_ptr.get(), other.pixel_ptr.get(), width*height*sizeof(T));
  }

  void SetUniform(T value) {
    for(int32_t i=0; i<width*height; i++){
      pixel_ptr[i] = value;
    }
  }
};

//Basically restricted version of Vec3 class for use as pixel
template<Scalar_t color_type>
class Pixel {
public:
  Vec3<color_type> channel;

  Pixel() {};
  Pixel(color_type r, color_type g, color_type b) : channel(Vec3<color_type>(r,g,b)) {};
  Pixel(Vec3<color_type> other) : channel(other) {};
  
  template<Scalar_t color_type_other>
  bool operator==(const Pixel<color_type_other>& rhs) const noexcept {
    return channel == rhs.channel;
  }
  template<Scalar_t color_type_other>
  bool operator!=(const Pixel<color_type_other>& rhs) const noexcept {
    return channel != rhs.channel;
  }
  template<Scalar_t color_type_other>
  Pixel<color_type> operator+(const Pixel<color_type_other>& rhs) const noexcept {
    return Pixel<color_type>(channel) += rhs.channel;
  }
  template<Scalar_t color_type_other>
  Pixel<color_type> operator+(Pixel<color_type_other>&& rhs) const noexcept {
    return Pixel<color_type>(channel) += rhs.channel;
  }
  template<Scalar_t color_type_other>
  Pixel<color_type> operator*(color_type_other rhs) const noexcept {
    return Pixel<color_type>(channel) *= rhs;
  }
  template<Scalar_t color_type_other>
  Pixel<color_type>& operator+=(const Pixel<color_type_other>& rhs) noexcept {
    channel += rhs.channel;
    return *this;
  }
  template<Scalar_t color_type_other>
  Pixel<color_type>& operator+=(Pixel<color_type_other>&& rhs) noexcept {
    channel += rhs.channel;
    return *this;
  }
  template<Scalar_t color_type_other>
  Pixel<color_type>& operator*=(color_type_other rhs) noexcept {
    channel *= rhs;
    return *this;
  }
  template<Scalar_t color_type_other>
  Pixel<color_type> operator-(const Pixel<color_type_other>& rhs) const noexcept {
    return Pixel<color_type>(channel) -= rhs.channel;
  }
  template<Scalar_t color_type_other>
  Pixel<color_type> operator-(Pixel<color_type_other>&& rhs) const noexcept {
    return Pixel<color_type>(channel) -= rhs.channel;
  }
  template<Scalar_t color_type_other>
  Pixel<color_type>& operator-=(const Pixel<color_type_other>& rhs) noexcept {
    channel -= rhs.channel;
    return *this;
  }
  template<Scalar_t color_type_other>
  Pixel<color_type>& operator-=(Pixel<color_type_other>&& rhs) noexcept {
    channel -= rhs.channel;
    return *this;
  }

  template<Scalar_t color_type_other>
  Pixel<color_type>& operator=(const Pixel<color_type_other>& other) {
    channel = other.channel;
    return *this;
  }
  
  template<Scalar_t color_type_other>
  color_type_other GetLinearGreyscale() const {
    return static_cast<color_type_other>(channel.x + channel.y + channel.z) /
           static_cast<color_type_other>(3);
  }
};