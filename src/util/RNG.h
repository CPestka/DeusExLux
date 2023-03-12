#pragma once

#include <random>
#include <cmath>
#include <memory>

#include "UtilConcepts.h"

std::mt19937 InitGenerator() {
  std::random_device random_dev;
  return std::mt19937(random_dev());
  //return std::mt19937();
}

template<Float_t float_type>
float_type GetInitialUniformRandom(){
  //I think intializing the distro and generator is expensive -> static
  static std::uniform_real_distribution<float_type> distribution(0.0, 1.0);
  static std::mt19937 generator = InitGenerator();
  return distribution(generator);
}

template<Float_t float_type>
std::unique_ptr<float_type[]> MakeRngTable(int size) {
  std::unique_ptr<float_type[]> rng_table = 
      std::make_unique<float_type[]>(size);
  for(int i=0; i<size; i++){
    rng_table[i] = GetInitialUniformRandom<float_type>();
  }
  return rng_table;
}

template<Float_t float_type>
class RNGLutEngine {
  public:
  int size;
  std::unique_ptr<float_type[]> rng_table;
  int i;
  RNGLutEngine(int size = 1024*16) : size(size), 
                              rng_table(MakeRngTable<float_type>(size)),
                              i(0) {};

  float_type GetUniformRandom() {
    float_type result = rng_table[i];
    i++;
    i = i == size ? 0 : i;
    return result;
  }

  float_type GetUniformRandom(float_type min, float_type max){
    return min + (max-min) * GetUniformRandom();
  }

  float_type GetUniformRandom(
      float_type min,
      bool including_min,
      float_type max,
      bool including_max){
  float_type tmp = GetUniformRandom(min, max);
  bool checked_boundaries = false;

  while (!checked_boundaries) {
    if ((!including_min && tmp == min) ||
        (!including_max && tmp == max)) {
      tmp = GetUniformRandom(min,max); //redraw if boundaries hit
    } else {
      checked_boundaries = true;
    }
  }
  return tmp;
  }
};