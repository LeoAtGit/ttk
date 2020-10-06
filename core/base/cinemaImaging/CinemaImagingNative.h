/// \ingroup base
/// \class ttk::EmbreeRenderer
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.05.2020
///
/// \brief TTK %EmbreeRenderer processing package.
///
/// %EmbreeRenderer is a TTK processing package that

#pragma once

#include <Debug.h>

namespace ttk {
  class CinemaImagingNative : virtual public Debug {
  public:
    CinemaImagingNative(){
        this->setDebugMsgPrefix("CinemaImaging(Native)");
    }
    ~CinemaImagingNative(){};

    template<typename IT>
    int renderImage(
      float* depthBuffer,
      unsigned int* primitiveIds,
      float* barycentricCoordinates,

      const size_t& nVertices,
      const float* vertexCoords,
      const size_t& nTriangles,
      const IT* connectivityList,

      const double resolution[2],
      const double camCenter[3],
      const double camDir[3],
      const double camUp[3],
      const double& camHeight
    ) const;
  };
};

template<typename IT>
int CinemaImagingNative::renderImage(
  float* depthBuffer,
  unsigned int* primitiveIds,
  float* barycentricCoordinates,

  const size_t& nVertices,
  const float* vertexCoords,
  const size_t& nTriangles,
  const IT* connectivityList,

  const double resolution[2],
  const double camCenter[3],
  const double camDir[3],
  const double camUp[3],
  const double& camHeight
) const {
  ttk::Timer timer;
  this->printMsg(
      "Rendering Image ("+std::string(orthographicProjection ? "O" : "P")+ "|"+std::to_string(resolution[0])+"x"+std::to_string(resolution[1])+")",
      0,0,this->threadNumber_,
      ttk::debug::LineMode::REPLACE
  );

  struct RTCIntersectContext context;
  rtcInitIntersectContext(&context);

  // Compute camera size
  const double aspect = resolution[0] / resolution[1];
  const double camSize[2] = { aspect * camHeight, camHeight};

  const auto normalize = [](double out[3], const double in[3]){
      double temp = sqrt(in[0] * in[0] + in[1] * in[1]
                         + in[2] * in[2]);
      out[0] = in[0]/temp;
      out[1] = in[1]/temp;
      out[2] = in[2]/temp;
  };

  double camDir[3]{0,0,0};
  normalize(camDir,camDirRaw);

  // Compute camRight = camDir x CamUp
  double camRight[3]{camDir[1] * camUp[2] - camDir[2] * camUp[1],
                        camDir[2] * camUp[0] - camDir[0] * camUp[2],
                        camDir[0] * camUp[1] - camDir[1] * camUp[0]};
  normalize(camRight,camRight);

  // Compute true up std::vector
  double camUpTrue[3]{camDir[1] * (-camRight[2]) - camDir[2] * (-camRight[1]),
       camDir[2] * (-camRight[0]) - camDir[0] * (-camRight[2]),
       camDir[0] * (-camRight[1]) - camDir[1] * (-camRight[0])};
  normalize(camUpTrue,camUpTrue);

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
  int resX = resolution[0];
  int resY = resolution[1];
  size_t pixelIndex = 0;
  size_t bcIndex = 0;

  for(int y = 0; y < resY; y++) {
      double v = ((double)y) * pixelHeightWorld;

      for(int x = 0; x < resX; x++) {
        double u = ((double)x) * pixelWidthWorld;

        // set origin
        // rayhit.ray.org_x = camPosCorner[0] + u * camRight[0] + v * camUpTrue[0];
        // rayhit.ray.org_y = camPosCorner[1] + u * camRight[1] + v * camUpTrue[1];
        // rayhit.ray.org_z = camPosCorner[2] + u * camRight[2] + v * camUpTrue[2];

        // set dir
        // rayhit.ray.dir_x = camDir[0];
        // rayhit.ray.dir_y = camDir[1];
        // rayhit.ray.dir_z = camDir[2];

        // compute ray hit:
        //  * iterate over every triangle
        bool wasHit = false;

        // if hit
        if(wasHit){
          // depthBuffer[pixelIndex] = std::max(0.0f,rayhit.ray.tfar);
          // primitiveIds[pixelIndex] = rayhit.hit.primID;
          // barycentricCoordinates[bcIndex++] = rayhit.hit.u;
          // barycentricCoordinates[bcIndex++] = rayhit.hit.v;
        } else {
          // depthBuffer[pixelIndex] = nan;
          // primitiveIds[pixelIndex] = CinemaImagingEmbree::INVALID_ID;
          // barycentricCoordinates[bcIndex++] = nan;
          // barycentricCoordinates[bcIndex++] = nan;
        }
        pixelIndex++;
      }
    }
};