#pragma once

#include <concepts>
#include <type_traits>

//Used to constrain the templates to sensible data types all over the app

template<typename T>
concept Integer_t = std::is_integral<T>::value;

template<typename T>
concept UInteger_t = std::is_integral<T>::value && std::is_unsigned<T>::value;

template<typename T>
concept Float_t = std::is_floating_point<T>::value;

template<typename T>
concept Scalar_t = std::is_arithmetic<T>::value;

template<typename T>
concept Color_discrete_t = UInteger_t<T> && (sizeof(T) <= 2);
