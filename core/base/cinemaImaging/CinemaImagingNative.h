/// \ingroup base
/// \class ttk::CinemaImagingNative
/// \authors Jonas Lukasczyk <jl@jluk.de>, Rosty Hnatyshyn
/// <rostyslav.hnatyshyn@gmail.com> \date 10.11.2020
///
/// \brief Native renderer. Uses a bounding volume hierarchy for
/// its acceleration structure.

#pragma once
#include "BVH.h"
#include "Ray.h"
#include <Debug.h>
namespace ttk {
  class CinemaImagingNative : virtual public Debug {
  public:
    CinemaImagingNative() {
      this->setDebugMsgPrefix("CinemaImaging(Native)");
    }
    ~CinemaImagingNative(){};

    template <typename IT>
    int renderImage(float *depthBuffer,
                    unsigned int *primitiveIds,
                    float *barycentricCoordinates,

                    const size_t &nVertices,
                    const float *vertexCoords,
                    const size_t &nTriangles,
                    const IT *connectivityList,

                    const BVH<IT> &bvh,

                    const double resolution[2],
                    const double camPos[3],
                    const double camDirRaw[3],
                    const double camUp[3],
                    const double &camHeight,
                    const bool &orthographicProjection,
                    const double &viewAngle) const;

    template <typename DT, typename IT>
    int interpolateArray(DT *outputArray,

                         const unsigned int *primitiveIds,
                         const float *barycentricCoordinates,
                         const IT *connectivityList,

                         const DT *inputArray,
                         const size_t &nTuples,
                         const size_t &nComponents = 1,
                         const DT &missingValue
                         = std::numeric_limits<DT>::has_quiet_NaN
                             ? std::numeric_limits<DT>::quiet_NaN()
                             : std::numeric_limits<DT>::max()) const;

    template <typename DT>
    int lookupArray(DT *outputArray,

                    const unsigned int *primitiveIds,

                    const DT *inputArray,
                    const size_t &nTuples,
                    const size_t &nComponents = 1,
                    const DT &missingValue
                    = std::numeric_limits<DT>::has_quiet_NaN
                        ? std::numeric_limits<DT>::quiet_NaN()
                        : std::numeric_limits<DT>::max()) const;
  };

}; // namespace ttk

template <typename DT, typename IT>
int ttk::CinemaImagingNative::interpolateArray(
  DT *outputArray,

  const unsigned int *primitiveIds,
  const float *barycentricCoordinates,
  const IT *connectivityList,

  const DT *inputArray,
  const size_t &nTuples,
  const size_t &nComponents,
  const DT &missingValue) const {

  if(nComponents != 1) {
    this->printErr(
      "Current implementation only supports interpolation of scalars.");
    return 0;
  }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
  for(size_t i = 0; i < nTuples; i++) {
    const unsigned int &cellId = primitiveIds[i];
    if(cellId == -1) {
      outputArray[i] = missingValue;
      continue;
    }

    const size_t cellIndex = cellId * 3;
    const IT &v0 = connectivityList[cellIndex + 0];
    const IT &v1 = connectivityList[cellIndex + 1];
    const IT &v2 = connectivityList[cellIndex + 2];

    const size_t bcIndex = i * 2;
    const float &u = barycentricCoordinates[bcIndex + 0];
    const float &v = barycentricCoordinates[bcIndex + 1];
    const float w = 1 - u - v;

    outputArray[i]
      = w * inputArray[v0] + u * inputArray[v1] + v * inputArray[v2];
  }

  return 1;
};

template <typename DT>
int ttk::CinemaImagingNative::lookupArray(DT *outputArray,

                                          const unsigned int *primitiveIds,

                                          const DT *inputArray,
                                          const size_t &nTuples,
                                          const size_t &nComponents,
                                          const DT &missingValue) const {

  for(size_t i = 0; i < nTuples; i++) {
    size_t outputOffset = i * nComponents;
    const unsigned int &cellId = primitiveIds[i];

    if(cellId == -1) {
      for(size_t j = 0; j < nComponents; j++)
        outputArray[outputOffset++] = missingValue;
      continue;
    } else {
      size_t inputOffset = cellId * nComponents;
      for(size_t j = 0; j < nComponents; j++)
        outputArray[outputOffset++] = inputArray[inputOffset++];
    }
  }

  return 1;
};

template <typename IT>
int ttk::CinemaImagingNative::renderImage(float *depthBuffer,
                                          unsigned int *primitiveIds,
                                          float *barycentricCoordinates,
                                          const size_t &nVertices,
                                          const float *vertexCoords,
                                          const size_t &nTriangles,
                                          const IT *connectivityList,
                                          const BVH<IT> &bvh,
                                          const double resolution[2],
                                          const double camPos[3],
                                          const double camDirRaw[3],
                                          const double camUp[3],
                                          const double &camHeight,
                                          const bool &orthographicProjection,
                                          const double &viewAngle) const {
  ttk::Timer timer;
  int resX = resolution[0];
  int resY = resolution[1];

  this->printMsg("Rendering Image ("
                   + std::string(orthographicProjection ? "O" : "P") + "|"
                   + std::to_string(resX) + "x" + std::to_string(resY) + ")",
                 0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);

  // Compute camera size
  const double aspect = resolution[0] / resolution[1];
  const double camSize[2] = {aspect * camHeight, camHeight};

  const auto normalize = [](double out[3], const double in[3]) {
    double temp = sqrt(in[0] * in[0] + in[1] * in[1] + in[2] * in[2]);
    out[0] = in[0] / temp;
    out[1] = in[1] / temp;
    out[2] = in[2] / temp;
  };

  double camDir[3]{0, 0, 0};
  normalize(camDir, camDirRaw);

  // Compute camRight = camDir x CamUp
  double camRight[3]{camDir[1] * camUp[2] - camDir[2] * camUp[1],
                     camDir[2] * camUp[0] - camDir[0] * camUp[2],
                     camDir[0] * camUp[1] - camDir[1] * camUp[0]};
  normalize(camRight, camRight);

  // Compute true up std::vector
  double camUpTrue[3]{camDir[1] * (-camRight[2]) - camDir[2] * (-camRight[1]),
                      camDir[2] * (-camRight[0]) - camDir[0] * (-camRight[2]),
                      camDir[0] * (-camRight[1]) - camDir[1] * (-camRight[0])};
  normalize(camUpTrue, camUpTrue);

  // Compute pixel size in world coordinates
  double pixelWidthWorld = camSize[0] / resolution[0];
  double pixelHeightWorld = camSize[1] / resolution[1];

  // Optimization: precompute half of the camera size to reduce the number of
  // operations in the for loop. Include a half pixel offset (-0.5) to center
  // vertices at pixel centers
  double camWidthWorldHalf = 0.5 * camSize[0] - 0.5 * pixelWidthWorld;
  double camHeightWorldHalf = 0.5 * camSize[1] - 0.5 * pixelHeightWorld;

  // Optimization: reorient camera model to bottom left corner to reduce
  // operations in for loop
  double camPosCorner[3] = {camPos[0] - camRight[0] * camWidthWorldHalf
                              - camUpTrue[0] * camHeightWorldHalf,
                            camPos[1] - camRight[1] * camWidthWorldHalf
                              - camUpTrue[1] * camHeightWorldHalf,
                            camPos[2] - camRight[2] * camWidthWorldHalf
                              - camUpTrue[2] * camHeightWorldHalf};

  float nan = std::numeric_limits<float>::quiet_NaN();
  if(orthographicProjection) {
    #ifdef TTK_ENABLE_OPENMP
    #pragma omp parallel for num_threads(this->threadNumber_)
    #endif
    for(int y = 0; y < resY; y++) {
      double v = ((double)y) * pixelHeightWorld;

      size_t pixelIndex = y * resX;
      size_t bcIndex = 2 * pixelIndex;

      for(int x = 0; x < resX; x++) {
        double u = ((double)x) * pixelWidthWorld;

        depthBuffer[pixelIndex] = nan;
        primitiveIds[pixelIndex] = -1;
        barycentricCoordinates[bcIndex] = nan;
        barycentricCoordinates[bcIndex + 1] = nan;

        // set origin
        float org_x = camPosCorner[0] + u * camRight[0] + v * camUpTrue[0];
        float org_y = camPosCorner[1] + u * camRight[1] + v * camUpTrue[1];
        float org_z = camPosCorner[2] + u * camRight[2] + v * camUpTrue[2];

        float ray_origin[3] = {org_x, org_y, org_z};

        // set dir
        float dir_x = camDir[0];
        float dir_y = camDir[1];
        float dir_z = camDir[2];

        float ray_dir[3] = {dir_x, dir_y, dir_z};

        Ray ray(ray_dir, ray_origin);
        bool wasHit = false;
        int triIdx;
        float distance;
        wasHit = bvh.intersect(
          ray, connectivityList, vertexCoords, &triIdx, &distance);
        if(wasHit) {
          depthBuffer[pixelIndex] = distance;
          primitiveIds[pixelIndex] = triIdx;
          barycentricCoordinates[bcIndex] = ray.u;
          barycentricCoordinates[bcIndex + 1] = ray.v;
        }
        pixelIndex++;
        bcIndex += 2;
      }
    }
  } else {
    double focalLength = camSize[0] / 2 / tan(viewAngle / 360.0 * 3.141592653589793);
    #ifdef TTK_ENABLE_OPENMP
    #pragma omp parallel for num_threads(this->threadNumber_)
    #endif
    for(int y = 0; y < resY; y++) {
      double v = ((double)y) * pixelHeightWorld;
      size_t pixelIndex = y * resX;
      size_t bcIndex = 2 * pixelIndex;

      for(int x = 0; x < resX; x++) {
        double u = ((double)x) * pixelWidthWorld;

        depthBuffer[pixelIndex] = nan;
        primitiveIds[pixelIndex] = -1;
        barycentricCoordinates[bcIndex] = nan;
        barycentricCoordinates[bcIndex + 1] = nan;

        // set origin
        float org_x = camPos[0];
        float org_y = camPos[1];
        float org_z = camPos[2];

        float ray_origin[3] = {org_x, org_y, org_z};
        // set dir
        float dir_x = camPosCorner[0] + u * camRight[0] + v * camUpTrue[0]
                      + camDir[0] * focalLength - org_x;
        float dir_y = camPosCorner[1] + u * camRight[1] + v * camUpTrue[1]
                      + camDir[1] * focalLength - org_y;
        float dir_z = camPosCorner[2] + u * camRight[2] + v * camUpTrue[2]
                      + camDir[2] * focalLength - org_z;

        float ray_dir[3] = {dir_x, dir_y, dir_z};

        Ray ray(ray_dir, ray_origin);
        bool wasHit = false;
        int triIdx;
        float distance;
        wasHit = bvh.intersect(ray, connectivityList, vertexCoords, &triIdx, &distance);
        if(wasHit) {
          depthBuffer[pixelIndex] = distance;
          primitiveIds[pixelIndex] = triIdx;
          barycentricCoordinates[bcIndex] = ray.u;
          barycentricCoordinates[bcIndex + 1] = ray.v;
        }
        pixelIndex++;
        bcIndex += 2;
      }
    }
  }
  this->printMsg("Rendering Image ("
                   + std::string(orthographicProjection ? "O" : "P") + "|"
                   + std::to_string(resX) + "x" + std::to_string(resY) + ")",

                 1, timer.getElapsedTime(), this->threadNumber_);

  return 1;
};
