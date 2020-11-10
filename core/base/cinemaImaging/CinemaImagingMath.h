#ifndef CINEMAIMAGINGMATH_H
#define CINEMAIMAGINGMATH_H
namespace ttk {
class CinemaImagingMath {
 public:
  static void multiplyByScalar(float* out, const float* a, const float& scalar);
  static void addVectors(float* out, const float* a, const float* b);
  static void subVectors(float* out, const float* a, const float* b);
  static void cross(float* out, const float* a, const float* b);
  static float dot(const float* a, const float* b);
};
}
#endif
