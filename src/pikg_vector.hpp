#ifndef HPP_PIKG_VECTOR
#define HPP_PIKG_VECTOR

#ifdef PIKG_USE_FDPS_VECTOR
#include <particle_simulator.hpp>
namespace PIKG = ParticleSimulator;

namespace ParticleSimulator{

#else

#include <iostream>
#include <iomanip>
#include <cassert>

namespace PIKG{
  using F64 = double;
  using F32 = float;
  //using F16 = float16_t;
  using U64 = uint64_t;
  using U32 = uint32_t;
  //using U16 = uint16_t;
  using S64 = int64_t;
  using S32 = int32_t;
  //using S16 = int16_t;

  template<class T>
  class Vector3{
  public:
    T x, y, z;
    Vector3() : x(T(0)), y(T(0)), z(T(0)) {}
    Vector3(const T _x, const T _y, const T _z) : x(_x), y(_y), z(_z) {}
    Vector3(const T s) : x(s), y(s), z(s) {}
    Vector3(const Vector3 & src) : x(src.x), y(src.y), z(src.z) {}

    typedef T DataType;
    static const int DIM = 3;

    const Vector3 & operator = (const Vector3 & rhs){
      x = rhs.x;
      y = rhs.y;
      z = rhs.z;
      return (*this);
    }

    const Vector3 & operator = (const T s){
      x = y = z = s;
      return (*this);
    }

    Vector3 operator + (const Vector3 & rhs) const{
      return Vector3(x + rhs.x, y + rhs.y, z + rhs.z);
    }
    const Vector3 & operator += (const Vector3 & rhs){
      (*this) = (*this) + rhs;
      return (*this);
    }
    Vector3 operator - (const Vector3 & rhs) const{
      return Vector3(x - rhs.x, y - rhs.y, z - rhs.z);
    }

    const Vector3 & operator -= (const Vector3 & rhs){
      (*this) = (*this) - rhs;
      return (*this);
    }

    Vector3 operator * (const T s) const{
      return Vector3(x * s, y * s, z * s);
    }
    const Vector3 & operator *= (const T s){
      (*this) = (*this) * s;
      return (*this);
    }
    friend Vector3 operator * (const T s, const Vector3 & v){
      return (v * s);
    }
    Vector3 operator / (const T s) const{
      return Vector3(x / s, y / s, z / s);
    }
    const Vector3 & operator /= (const T s){
      (*this) = (*this) / s;
      return (*this);
    }

    const Vector3 & operator + () const {
      return (* this);
    }

    const Vector3 operator - () const {
      return Vector3(-x, -y, -z);
    }

    T operator * (const Vector3 & rhs) const{
      return (x * rhs.x) + (y * rhs.y) + (z * rhs.z);
    }

    Vector3 operator ^ (const Vector3 & rhs) const{
      return Vector3( (y * rhs.z - z * rhs.y),
		      (z * rhs.x - x * rhs.z),
		      (x * rhs.y - y * rhs.x) );
    }

    template <typename U>
    operator Vector3<U> () const {
      return Vector3<U> (static_cast<U>(x),
			 static_cast<U>(y),
			 static_cast<U>(z));
    }

    friend std::ostream & operator <<(std::ostream & c, const Vector3 & u){
      c<<u.x<<"   "<<u.y<<"    "<<u.z;
      return c;
    }

    friend std::istream & operator >>(std::istream & c, Vector3 & u){
      c>>u.x; c>>u.y; c>>u.z;
      return c;
    }

    const T & operator[](const int i) const {
      assert(i>=0 && i<3);
      if(0==i) return x;
      if(1==i) return y;
      return z;
    }

    T & operator[](const int i){
      assert(i>=0 && i<3);
      if(0==i) return x;
      if(1==i) return y;
      return z;
    }

    bool operator == (const Vector3 & u) const {
      return ( (x==u.x) && (y==u.y) && (z==u.z) );
    }
    bool operator != (const Vector3 & u) const {
      return ( (x!=u.x) || (y!=u.y) || (z!=u.z) );
    }
  };

  template<typename T>
  inline T max(const Vector3<T>& v){
    T max_val = (v.x > v.y) ? v.x : v.y;
    max_val = (max_val > v.z ) ? max_val : v.z;
    return max_val;
  }

  template<typename T>
  inline T min(const Vector3<T>& v){
    T min_val = (v.x < v.y) ? v.x : v.y;
    min_val = (min_val < v.z ) ? min_val : v.z;
    return min_val;
  }

  template <>
  inline Vector3<float> Vector3<float>::operator / (const float s) const {
    const float inv_s = 1.0f/s;
    return Vector3(x * inv_s, y * inv_s, z * inv_s);
  }
  template <>
  inline Vector3<double> Vector3<double>::operator / (const double s) const {
    const double inv_s = 1.0/s;
    return Vector3(x * inv_s, y * inv_s, z * inv_s);
  }

  using F64vec = Vector3<F64>;
  using F32vec = Vector3<F32>;
  using F64vec3 = Vector3<F64>;
  using F32vec3 = Vector3<F32>;
  //using F16vec = Vector3<F16>;

  template <typename T>
  class Vector2{
  public:
    T x, y;
    Vector2() : x(T(0)), y(T(0)) {}
    Vector2(const T _x, const T _y) : x(_x), y(_y) {}
    Vector2(const T s) : x(s), y(s) {}
    Vector2(const Vector2 & src) : x(src.x), y(src.y) {}

    static const int DIM = 2;
	
    const Vector2 & operator = (const Vector2 & rhs){
      x = rhs.x;
      y = rhs.y;
      return (*this);
    }

    const Vector2 & operator = (const T s){
      x = y = s;
      return (*this);
    }

    Vector2 operator + (const Vector2 & rhs) const{
      return Vector2(x + rhs.x, y + rhs.y);
    }
    const Vector2 & operator += (const Vector2 & rhs){
      (*this) = (*this) + rhs;
      return (*this);
    }
    Vector2 operator - (const Vector2 & rhs) const{
      return Vector2(x - rhs.x, y - rhs.y);
    }
    const Vector2 & operator -= (const Vector2 & rhs){
      (*this) = (*this) - rhs;
      return (*this);
    }

    // vector scholar products
    Vector2 operator * (const T s) const{
      return Vector2(x * s, y * s);
    }
    const Vector2 & operator *= (const T s){
      (*this) = (*this) * s;
      return (*this);
    }
    friend Vector2 operator * (const T s, const Vector2 & v){
      return (v * s);
    }
    Vector2 operator / (const T s) const{
      return Vector2(x / s, y / s);
    }
    const Vector2 & operator /= (const T s){
      (*this) = (*this) / s;
      return (*this);
    }

    const Vector2 & operator + () const {
      return (* this);
    }

    const Vector2 operator - () const {
      return Vector2(-x, -y);
    }

    // inner product
    T operator * (const Vector2 & rhs) const{
      return (x * rhs.x) + (y * rhs.y);
    }

    // outer product (retruned value is scholar)
    T operator ^ (const Vector2 & rhs) const{
      const T z = (x * rhs.y) - (y * rhs.x);
      return z;
    }

    //cast to Vector2<U>
    template <typename U>
    operator Vector2<U> () const {
      return Vector2<U> (static_cast<U>(x),
			 static_cast<U>(y));
    }

    T getMax() const {
      return x > y ? x : y;
    }

    T getMin() const {
      return x < y ? x : y;
    }

    template <class F>
    Vector2 applyEach(F f) const {
      return Vector2(f(x), f(y));
    }
    template <class F>
    friend Vector2 ApplyEach(F f, const Vector2 & arg1, const Vector2 & arg2){
      return Vector2( f(arg1.x, arg2.x), f(arg1.y, arg2.y) );
    }

    friend std::ostream & operator <<(std::ostream & c, const Vector2 & u){
      c<<u.x<<"   "<<u.y;
      return c;
    }

    friend std::istream & operator >>(std::istream & c, Vector2 & u){
      c>>u.x; c>>u.y;
      return c;
    }

    T& operator[](const int i){
      assert(i>=0 && i<2);
      if(0==i) return x;
      return y;
    }
    const T & operator[](const int i) const {
      assert(i>=0 && i<2);
      if(0==i) return x;
      return y;
    }


    bool operator == (const Vector2 & u) const {
      return ( (x==u.x) && (y==u.y) );
    }
    bool operator != (const Vector2 & u) const {
      return ( (x!=u.x) || (y!=u.y) );
    }

  };

  template <>
  inline Vector2<float> Vector2<float>::operator / (const float s) const {
    const float inv_s = 1.0f/s;
    return Vector2(x * inv_s, y * inv_s);
  }
  template <>
  inline Vector2<double> Vector2<double>::operator / (const double s) const {
    const double inv_s = 1.0/s;
    return Vector2(x * inv_s, y * inv_s);
  }

  using F64vec2 = Vector2<F64>;
  using F32vec2 = Vector2<F32>;
  //using F16vec2 = Vector2<float16_t>;

#endif

  template<class T>
  class Vector4{
  public:
    T x, y, z, w;
    Vector4() : x(T(0)), y(T(0)), z(T(0)), w(T(0)) {}
    Vector4(const T _x, const T _y, const T _z, const T _w) : x(_x), y(_y), z(_z), w(_w) {}
    Vector4(const T s) : x(s), y(s), z(s), w(s) {}
    Vector4(const Vector4 & src) : x(src.x), y(src.y), z(src.z), w(src.w) {}

    typedef T DataType;
    static const int DIM = 4;

    const Vector4 & operator = (const Vector4 & rhs){
      x = rhs.x;
      y = rhs.y;
      z = rhs.z;
      w = rhs.w;
      return (*this);
    }

    const Vector4 & operator = (const T s){
      x = y = z = w = s;
      return (*this);
    }

    Vector4 operator + (const Vector4 & rhs) const{
      return Vector4(x + rhs.x, y + rhs.y, z + rhs.z, w + rhs.z);
    }
    const Vector4 & operator += (const Vector4 & rhs){
      (*this) = (*this) + rhs;
      return (*this);
    }
    Vector4 operator - (const Vector4 & rhs) const{
      return Vector4(x - rhs.x, y - rhs.y, z - rhs.z, w - rhs.w);
    }

    const Vector4 & operator -= (const Vector4 & rhs){
      (*this) = (*this) - rhs;
      return (*this);
    }

    Vector4 operator * (const T s) const{
      return Vector4(x * s, y * s, z * s, w * s);
    }
    const Vector4 & operator *= (const T s){
      (*this) = (*this) * s;
      return (*this);
    }
    friend Vector4 operator * (const T s, const Vector4 & v){
      return (v * s);
    }
    Vector4 operator / (const T s) const{
      return Vector4(x / s, y / s, z / s, w / s);
    }
    const Vector4 & operator /= (const T s){
      (*this) = (*this) / s;
      return (*this);
    }

    const Vector4 & operator + () const {
      return (* this);
    }

    const Vector4 operator - () const {
      return Vector4(-x, -y, -z, -w);
    }

    T operator * (const Vector4 & rhs) const{
      return (x * rhs.x) + (y * rhs.y) + (z * rhs.z) + (w * rhs.w);
    }


    template <typename U>
    operator Vector4<U> () const {
      return Vector4<U> (static_cast<U>(x),
			 static_cast<U>(y),
			 static_cast<U>(z),
			 static_cast<U>(w));
    }

    friend std::ostream & operator <<(std::ostream & c, const Vector4 & u){
      c<<u.x<<"   "<<u.y<<"    "<<u.z<<"   "<<u.w;
      return c;
    }

    friend std::istream & operator >>(std::istream & c, Vector4 & u){
      c>>u.x; c>>u.y; c>>u.z; c>>u.w;
      return c;
    }

    const T & operator[](const int i) const {
      assert(i>=0 && i<3);
      if(0==i) return x;
      if(1==i) return y;
      if(2==i) return z;
      return w;
    }

    T & operator[](const int i){
      assert(i>=0 && i<3);
      if(0==i) return x;
      if(1==i) return y;
      if(2==i) return z;
      return w;
    }

    bool operator == (const Vector4 & u) const {
      return ( (x==u.x) && (y==u.y) && (z==u.z) && (w==u.w));
    }
    bool operator != (const Vector4 & u) const {
      return ( (x!=u.x) || (y!=u.y) || (z!=u.z) || (w!=u.w));
    }
  };

  template<typename T>
  inline T max(const Vector4<T>& v){
    T max_val = (v.x > v.y) ? v.x : v.y;
    max_val = (max_val > v.z ) ? max_val : v.z;
    max_val = (max_val > v.w ) ? max_val : v.w;
    return max_val;
  }

  template<typename T>
  inline T min(const Vector4<T>& v){
    T min_val = (v.x < v.y) ? v.x : v.y;
    min_val = (min_val < v.z ) ? min_val : v.z;
    min_val = (min_val < v.w ) ? min_val : v.w;
    return min_val;
  }

  template <>
  inline Vector4<float> Vector4<float>::operator / (const float s) const {
    const float inv_s = 1.0f/s;
    return Vector4(x * inv_s, y * inv_s, z * inv_s, w * inv_s);
  }
  template <>
  inline Vector4<double> Vector4<double>::operator / (const double s) const {
    const double inv_s = 1.0/s;
    return Vector4(x * inv_s, y * inv_s, z * inv_s, w * inv_s);
  }

  using F64vec4 = Vector4<F64>;
  using F32vec4 = Vector4<F32>;
  //using F16vec = Vector4<F16>;
};

#endif
