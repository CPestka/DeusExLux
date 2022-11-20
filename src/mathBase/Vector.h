#pragma once

#include <cmath>
#include <optional>

#include "UtilConcepts.h"

template<Scalar_t T>
class Quaternion;

template<Scalar_t T>
class Vec2 {
  public:
  T x;
  T y;

  Vec2(){};
  Vec2(T x, T y) : x(x), y(y) {};
  template<Scalar_t V>
  Vec2(const Vec2<V>& other) : x(other.x), y(other.y) {};
  
  template<Scalar_t V>
  bool operator==(const Vec2<V>& other) const noexcept {
    return x == other.x && y == other.y;
  }
  template<Scalar_t V>
  bool operator!=(const Vec2<V>& other) const noexcept {
    return x != other.x && y != other.y;
  }

  template<Scalar_t V>
  Vec2<T> operator+(const Vec2<V>& other) const noexcept {
    return Vec2<T>(x+other.x,y+other.y);
  }
  template<Scalar_t V>
  Vec2<T> operator+(Vec2<V>&& other) const noexcept {
    return Vec2<T>(x+other.x,y+other.y);
  }
  template<Scalar_t V>
  Vec2<T>& operator+=(const Vec2<V>& other) noexcept {
    x += other.x;
    y += other.y;
    return *this;
  }
  template<Scalar_t V>
  Vec2<T>& operator+=(Vec2<V>&& other) noexcept {
    x += other.x;
    y += other.y;
    return *this;
  }
  template<Scalar_t V>
  Vec2<T> operator-(const Vec2<V>& other) const noexcept {
    return Vec2<T>(x-other.x,y-other.y);
  }
  template<Scalar_t V>
  Vec2<T> operator-(Vec2<V>&& other) const noexcept {
    return Vec2<T>(x-other.x,y-other.y);
  }
  template<Scalar_t V>
  Vec2<T>& operator-=(const Vec2<V>& other) noexcept {
    x -= other.x;
    y -= other.y;
    return *this;
  }
  template<Scalar_t V>
  Vec2<T>& operator-=(Vec2<V>&& other) noexcept {
    x -= other.x;
    y -= other.y;
    return *this;
  }
  //Scalar Product
  template<Scalar_t V>
  T operator*(const Vec2<V>& other) const noexcept {
    return x * other.x + y * other.y;
  }
  template<Scalar_t V>
  T operator*(Vec2<V>&& other) const noexcept {
    return x * other.x + y * other.y;
  }
  
  template<Scalar_t V>
  Vec2<T> operator+(V other) const noexcept {
    return Vec2<T>(x+other,y+other);
  }
  template<Scalar_t V>
  Vec2<T>& operator+=(V other) noexcept {
    x += other;
    y += other;
    return *this;
  }
  template<Scalar_t V>
  Vec2<T> operator-(V other) const noexcept {
    return Vec2<T>(x-other,y-other);
  }
  template<Scalar_t V>
  Vec2<T>& operator-=(V other) noexcept {
    x -= other;
    y -= other;
    return *this;
  }
  //Scalar multiplication
  template<Scalar_t V>
  Vec2<T> operator*(V other) const noexcept {
    return Vec2<T>(x*other,y*other);
  }
  template<Scalar_t V>
  Vec2<T>& operator*=(V other) noexcept {
    x *= other;
    y *= other;
    return *this;
  }
  //Scalar devision
  template<Scalar_t V>
  Vec2<T> operator/(V other) const noexcept {
    return Vec2<T>(x/other,y/other);
  }

  //Util functions

  T GetLength() const {
    return std::sqrt((*this) * (*this));
  }

  bool TryNormalize() {
    T r = GetLength();
    if (r != 0) {
      *this = (*this) / r;
      return true; 
    } else {
      return false;
    }
  }
  void Normalize() {
    *this = (*this) / GetLength();
  }

  void ConvertToPolarCoords() {
    x = GetLength();
    y = asin(y);
  }
  void ConvertToKartesianCoords() {
    T r = x;
    x = r * cos(y);
    y = r * sin(y);
  }

  void RotateCartesian(T phi) {
    ConvertToPolarCoords();
    y += phi;
    ConvertToKartesianCoords();
  }
};

//Lhs versions of operator overloads 
template<Scalar_t T, Scalar_t V>
Vec2<T> operator+(V lhs, const Vec2<T>& rhs) noexcept {
  return Vec2<T>(lhs + rhs.x, lhs + rhs.y);
}
template<Scalar_t T, Scalar_t V>
Vec2<T> operator+(V lhs, Vec2<T>&& rhs) noexcept {
  return Vec2<T>(lhs + rhs.x, lhs + rhs.y);
}
//+= and -= lhs versions ommited as they dont make sense
template<Scalar_t T, Scalar_t V>
Vec2<T> operator-(V lhs, const Vec2<T>& rhs) noexcept {
  return Vec2<T>(lhs - rhs.x, lhs - rhs.y);
}
template<Scalar_t T, Scalar_t V>
Vec2<T> operator-(V lhs, Vec2<T>&& rhs) noexcept {
  return Vec2<T>(lhs - rhs.x, lhs - rhs.y);
}
//No scalar devided vec overload as it doesnt make sense
template<Scalar_t T, Scalar_t V>
Vec2<T> operator*(V lhs, Vec2<T>&& rhs) noexcept {
  return Vec2<T>(lhs * rhs.x, lhs * rhs.y);
}
template<Scalar_t T, Scalar_t V>
Vec2<T> operator*(V lhs, const Vec2<T>& rhs) noexcept {
  return Vec2<T>(lhs * rhs.x, lhs * rhs.y);
}

template<Scalar_t T>
Vec2<T> ComputeNormal(const Vec2<T>& v) {
  return Vec2<T>(v) / v.GetLength();
}
template<Scalar_t T>
Vec2<T> ComputeNormal(Vec2<T>&& v) {
  return Vec2<T>(v) / v.GetLength();
}
template<Scalar_t T>
std::optional<Vec2<T>> TryComputeNormal(const Vec2<T>& v) {
  T r = v.GetLength();
  return r != 0 ? Vec2<T>(v) / r : std::nullopt;
}
template<Scalar_t T>
std::optional<Vec2<T>> TryComputeNormal(Vec2<T>&& v) {
  T r = v.GetLength();
  return r != 0 ? Vec2<T>(v) / r : std::nullopt;
}

template<Scalar_t T>
class Vec3 {
  public:
  T x;
  T y;
  T z;

  Vec3(){};
  Vec3(T x, T y, T z) : x(x), y(y), z(z) {};
  template<Scalar_t V>
  Vec3(const Vec3<V>& other) : x(other.x), y(other.y), z(other.z) {};
  template<Scalar_t V>
  explicit Vec3(const Quaternion<V>& q) : x(q.x2), y(q.x3), z(q.x4) {};

  template<Scalar_t V>
  bool operator==(const Vec3<V>& other) const noexcept {
    return x == other.x && y == other.y && z == other.z;
  }
  template<Scalar_t V>
  bool operator!=(const Vec3<V>& other) const noexcept {
    return x != other.x && y != other.y && z != other.z;
  }

  template<Scalar_t V>
  Vec3<T> operator+(const Vec3<V>& other) const noexcept {
    return Vec3<T>(x+other.x,y+other.y,z+other.z);
  }
  template<Scalar_t V>
  Vec3<T> operator+(Vec3<V>&& other) const noexcept {
    return Vec3<T>(x+other.x,y+other.y,z+other.z);
  }
  template<Scalar_t V>
  Vec3<T>& operator+=(const Vec3<V>& other) noexcept {
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
  }
  template<Scalar_t V>
  Vec3<T>& operator+=(Vec3<V>&& other) noexcept {
    x += other.x;
    y += other.y;
    z += other.z;
    return *this;
  }
  template<Scalar_t V>
  Vec3<T> operator-(const Vec3<V>& other) const noexcept {
    return Vec3<T>(x-other.x,y-other.y,z-other.z);
  }
  template<Scalar_t V>
  Vec3<T> operator-(Vec3<V>&& other) const noexcept {
    return Vec3<T>(x-other.x,y-other.y,z-other.z);
  }
  template<Scalar_t V>
  Vec3<T>& operator-=(const Vec3<V>& other) noexcept {
    x -= other.x;
    y -= other.y;
    z -= other.z;
    return *this;
  }
  template<Scalar_t V>
  Vec3<T>& operator-=(Vec3<V>&& other) noexcept {
    x -= other.x;
    y -= other.y;
    z -= other.z;
    return *this;
  }
  //Scalar Product
  template<Scalar_t V>
  T operator*(const Vec3<V>& other) const noexcept {
    return x * other.x + y * other.y + z * other.z;
  }
  template<Scalar_t V>
  T operator*(Vec3<V>&& other) const noexcept {
    return x * other.x + y * other.y + z * other.z;
  }
  
  template<Scalar_t V>
  Vec3<T> operator+(V other) const noexcept {
    return Vec3<T>(x+other,y+other,z+other);
  }
  template<Scalar_t V>
  Vec3<T>& operator+=(V other) noexcept {
    x += other;
    y += other;
    z += other;
    return *this;
  }
  template<Scalar_t V>
  Vec3<T> operator-(V other) const noexcept {
    return Vec3<T>(x-other,y-other,z-other);
  }
  template<Scalar_t V>
  Vec3<T>& operator-=(V other) noexcept {
    x -= other;
    y -= other;
    z -= other;
    return *this;
  }
  //Scalar multiplication
  template<Scalar_t V>
  Vec3<T> operator*(V other) const noexcept {
    return Vec3<T>(x*other,y*other,z*other);
  }
  template<Scalar_t V>
  Vec3<T>& operator*=(V other) noexcept {
    x *= other;
    y *= other;
    z *= other;
    return *this;
  }
  //Scalar devision
  template<Scalar_t V>
  Vec3<T> operator/(V other) const noexcept {
    return Vec3<T>(x/other,y/other,z/other);
  }

  T GetLength() const {
    return std::sqrt((*this) * (*this));
  }

  bool TryNormalize() {
    T r = GetLength();
    if (r != 0) {
      *this = (*this) / r;
      return true; 
    } else {
      return false;
    }
  }
  void Normalize() {
    *this = (*this) / GetLength();
  }

  void ConvertToSphericalCoords() {
    T r = GetLength();
    T theta = std::atan(std::sqrt(x*x + y*y) / z);
    T phi = std::atan2(y,x);
    x = r;
    y = theta;
    z = phi;
  }
  void ConvertToCartisianCoords() {
    T tmp_x = x * std::sin(y) * std::cos(z);
    T tmp_y = x * std::sin(y) * std::sin(z);
    T tmp_z = x * std::cos(y);
    x = tmp_x;
    y = tmp_y;
    z = tmp_z;
  }

  void RotateCartisianVecPolar(T theta) {
    ConvertToSphericalCoords();
    y += theta;
    ConvertToCartisianCoords();
  }
  void RotateCartisianVecAzimuth(T phi) {
    ConvertToSphericalCoords();
    z += phi;
    ConvertToCartisianCoords();
  }
  void RotateCartisianVec(T theta, T phi) {
    ConvertToSphericalCoords();
    y += theta;
    z += phi;
    ConvertToCartisianCoords();
  }

  template<Scalar_t V, Scalar_t Z>
  bool TryRotateAroundVec(const Vec3<V>& look_at, Z phi) {
    if (GetLength() == 0 || look_at.GetLength() == 0) {
      return false;
    }

    Quaternion<T> vec_q(*this);
    Quaternion<T> look_at_q(ComputeNormal(look_at), phi);
    Quaternion<T> inv_look_at_q(look_at_q);
    inv_look_at_q.Invert();

    *this = Vec3<T>((inv_look_at_q * vec_q) * look_at_q);
    return true;
  }
  template<Scalar_t V, Scalar_t Z>
  void RotateAroundVec(const Vec3<V>& look_at, Z phi) {
    Quaternion<T> vec_q(*this);
    Quaternion<T> look_at_q(ComputeNormal(look_at), phi);
    Quaternion<T> inv_look_at_q(look_at_q);
    inv_look_at_q.Invert();

    *this = Vec3<T>((inv_look_at_q * vec_q) * look_at_q);
  }
  template<Scalar_t V>
  void RotateAroundVec(const Quaternion<V>& look_at_q) {
    Quaternion<T> vec_q(*this);
    Quaternion<T> inv_look_at_q(look_at_q);
    inv_look_at_q.Invert();

    *this = Vec3<T>((inv_look_at_q * vec_q) * look_at_q);
  }
  template<Scalar_t V, Scalar_t Z>
  void RotateAroundVec(
      const Quaternion<V>& look_at_q,
      const Quaternion<Z>& inv_look_at_q) {
    Quaternion<T> vec_q(*this);

    *this = Vec3<T>((inv_look_at_q * vec_q) * look_at_q);
  }
  
  //Performs the rotation on this vector that would align the vector "other"
  //with the y axis (0.0,1.0,0.0)
  template<Scalar_t V>
  void PerformAlignRotation(Vec3<V> other_in_cartesian_coords) {
    other_in_cartesian_coords.ConvertToSphericalCoords();
    ConvertToSphericalCoords();
    y -= other_in_cartesian_coords.y; //add theta
    z -= other_in_cartesian_coords.z; //add phi
    ConvertToCartisianCoords();
  }
  template<Scalar_t V>
  bool TryPerformAlignRotation(Vec3<V> other_in_cartesian_coords) {
    if (other_in_cartesian_coords.GetLength() == 0 ||
        GetLength() == 0) {
      return false;
    }
    other_in_cartesian_coords.ConvertToSphericalCoords();
    ConvertToSphericalCoords();
    y -= other_in_cartesian_coords.y; //add theta
    z -= other_in_cartesian_coords.z; //add phi
    ConvertToCartisianCoords();
    return true;
  }
};

template<Scalar_t T, Scalar_t V, Scalar_t Z, Scalar_t U>
Vec3<T> GetRotatedVec(
      const Quaternion<V>& vec_q,
      const Quaternion<Z>& look_at_q,
      const Quaternion<U>& inv_look_at_q) {
  return Vec3<T>((inv_look_at_q * vec_q) * look_at_q);
}
template<Scalar_t T, Scalar_t V, Scalar_t Z, Scalar_t U>
void RotateVec(
      Vec3<V>& vec,
      const Quaternion<Z>& look_at_q,
      const Quaternion<U>& inv_look_at_q) {
  Quaternion<T> vec_q(vec);
  vec = Vec3<T>((inv_look_at_q * vec_q) * look_at_q);
}
template<Scalar_t T, Scalar_t V, Scalar_t Z, Scalar_t U>
void RotateVec(
      Vec3<V>& vec,
      Z phi,
      const Vec3<U>& look_at) {
  Quaternion<T> vec_q(vec);
  Quaternion<T> look_at_q(ComputeNormal(look_at), phi);
  Quaternion<T> inv_look_at_q(look_at_q);
  inv_look_at_q.Invert();
  vec = Vec3<T>((inv_look_at_q * vec_q) * look_at_q);
}

//Lhs versions of operator overloads 
template<Scalar_t T, Scalar_t V>
Vec3<T> operator+(V lhs, const Vec3<T>& rhs) noexcept {
  return Vec3<T>(lhs + rhs.x, lhs + rhs.y);
}
template<Scalar_t T, Scalar_t V>
Vec3<T> operator+(V lhs, Vec3<T>&& rhs) noexcept {
  return Vec3<T>(lhs + rhs.x, lhs + rhs.y);
}
//+= and -= lhs versions ommited as they are probably likely to be mistakes
template<Scalar_t T, Scalar_t V>
Vec3<T> operator-(V lhs, const Vec3<T>& rhs) noexcept {
  return Vec3<T>(lhs - rhs.x, lhs - rhs.y);
}
template<Scalar_t T, Scalar_t V>
Vec3<T> operator-(V lhs, Vec3<T>&& rhs) noexcept {
  return Vec3<T>(lhs - rhs.x, lhs - rhs.y);
}
//No scalar devided vec overload as it doesnt make sense
template<Scalar_t T, Scalar_t V>
Vec3<T> operator*(V lhs, Vec3<T>&& rhs) noexcept {
  return Vec3<T>(lhs * rhs.x, lhs * rhs.y, lhs * rhs.z);
}
template<Scalar_t T, Scalar_t V>
Vec3<V> operator*(V lhs, const Vec3<T>& rhs) noexcept {
  return Vec3<V>(lhs * rhs.x, lhs * rhs.y, lhs * rhs.z);
}

template<Scalar_t T>
Vec3<T> ComputeNormal(const Vec3<T>& v) {
  return Vec3<T>(v) / v.GetLength();
}

template<Scalar_t T>
std::optional<Vec3<T>> TryComputeNormal(const Vec3<T>& v) {
  T r = v.GetLength();
  return r != 0 ? Vec3<T>(v) / r : std::nullopt;
}

template<Scalar_t T, Scalar_t V>
Vec3<T> ComputeCrossProduct(const Vec3<T>& v1,const Vec3<V>& v2) {
  return Vec3<T>(
      v1.y * v2.z - v1.z * v2.y,
      v1.z * v2.x - v1.x * v2.z,
      v1.x * v2.y - v1.y * v2.x);
}

template<Scalar_t T, Scalar_t V>
T GetAngleBetween(const Vec2<T>& lhs, const Vec2<V>& rhs) {
  return std::acos((lhs * rhs) / (lhs.GetLength() * rhs.GetLength()));
}
template<Scalar_t T, Scalar_t V>
std::optional<T> TryGetAngleBetween(const Vec2<T>& lhs, const Vec2<V>& rhs) {
  if (lhs.GetLength() == 0 || rhs.GetLength() == 0) {
    return std::nullopt;
  }
  return std::acos((lhs * rhs) / (lhs.GetLength() * rhs.GetLength()));
}
template<Scalar_t T, Scalar_t V>
T GetAngleBetween(const Vec3<T>& lhs, const Vec3<V>& rhs) {
  return std::acos((lhs * rhs) / (lhs.GetLength() * rhs.GetLength()));
}
template<Scalar_t T, Scalar_t V>
std::optional<T> TryGetAngleBetween(const Vec3<T>& lhs, const Vec3<V>& rhs) {
  if (lhs.GetLength() == 0 || rhs.GetLength() == 0) {
    return std::nullopt;
  }
  return std::acos((lhs * rhs) / (lhs.GetLength() * rhs.GetLength()));
}

template<Scalar_t T>
class Quaternion {
  public:
  T x1,x2,x3,x4;
  explicit Quaternion(T x1, T x2, T x3, T x4) : x1(x1), x2(x2), x3(x3), x4(x4) {}
  //For quaternion that is to be rotated
  template<Scalar_t V>
  explicit Quaternion(Vec3<V> v) {
    x1 = 0;
    x2 = v.x;
    x3 = v.y;
    x4 = v.z;
  }
  //For Rotation quaternion
  template<Scalar_t V, Scalar_t Z>
  explicit Quaternion(Vec3<V> v, Z phi) {
    x1 = std::cos(phi / 2.0);
    x2 = v.x * std::sin(phi / 2.0);
    x3 = v.y * std::sin(phi / 2.0);
    x4 = v.z * std::sin(phi / 2.0);
  }
  template<Scalar_t V>
  explicit Quaternion(const Quaternion<V>& other) {
    x1 = other.x1;
    x2 = other.x2;
    x3 = other.x3;
    x4 = other.x4;
  }

  void Invert() {
    x2 = -x2;
    x3 = -x3;
    x4 = -x4;
  }
  
  template<Scalar_t V>
  Quaternion<T> operator*(const Quaternion<V>& rhs) const noexcept {
    return Quaternion(
        x1 * rhs.x1 - x2 * rhs.x2 - x3 * rhs.x3 - x4 * rhs.x4,
        x1 * rhs.x2 + x2 * rhs.x1 - x3 * rhs.x4 + x4 * rhs.x3,
        x1 * rhs.x3 + x2 * rhs.x4 + x3 * rhs.x1 - x4 * rhs.x2,
        x1 * rhs.x4 - x2 * rhs.x3 + x3 * rhs.x2 + x4 * rhs.x1);
  }
  template<Scalar_t V>
  Quaternion<T>& operator*=(const Quaternion<V>& rhs) noexcept {
    T tmp1 = x1 * rhs.x1 - x2 * rhs.x2 - x3 * rhs.x3 - x4 * rhs.x4;
    T tmp2 = x1 * rhs.x2 + x2 * rhs.x1 - x3 * rhs.x4 + x4 * rhs.x3;
    T tmp3 = x1 * rhs.x3 + x2 * rhs.x4 + x3 * rhs.x1 - x4 * rhs.x2;
    T tmp4 = x1 * rhs.x4 - x2 * rhs.x3 + x3 * rhs.x2 + x4 * rhs.x1;
    x1 = tmp1;
    x2 = tmp2;
    x3 = tmp3;
    x4 = tmp4;
    return *this;
  }
};