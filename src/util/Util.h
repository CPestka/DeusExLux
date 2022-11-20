#pragma once

#include "UtilConcepts.h"

template<Float_t T>
T RadToDegrees(T angle) {
  return angle * 180.0 / M_PI;
}

template<Float_t T>
T DegreesToRad(T angle) {
  return angle * M_PI / 180.0;
}