#pragma once

#include <random>
#include <cmath>

#include "UtilConcepts.h"

template<Float_t float_type>
float_type GetUniformRandom(){
  //I think intializing the distro and generator is expensive -> static
  static std::uniform_real_distribution<float_type> distribution(0.0, 1.0);
  static std::mt19937 generator;
  return distribution(generator);
}

template<Float_t float_type>
float_type GetUniformRandom(float_type min, float_type max){
  return min + (max-min) * GetUniformRandom<float_type>();
}

//Overload to exclude boundaries of the intervall 
template<Float_t float_type>
float_type GetUniformRandom(float_type min,
                            bool including_min,
                            float_type max,
                            bool including_max){
  float_type tmp = GetUniformRandom(min, max);
  bool checked_boundaries = false;

  while (!checked_boundaries) {
    if (!including_min || tmp == min ||
        !including_max || tmp == max) {
      tmp = GetUniformRandom(min,max); //redraw if boundaries hit
    } else {
      checked_boundaries = true;
    }
  }
  return tmp;
}