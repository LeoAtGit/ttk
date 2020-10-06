/// \ingroup base
/// \class ttk::CinemaImagingEmbree
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.05.2020
///
/// \brief TTK %CinemaImagingEmbree processing package.
///
/// %CinemaImagingEmbree is a TTK processing package that

#pragma once

#if TTK_ENABLE_EMBREE

// base code includes
#include <Debug.h>
#include <string>
#include <vector>

#include <embree3/rtcore.h>

#include <limits>

namespace ttk {
  class CinemaImagingEmbree : virtual public Debug {
  public:
    static const unsigned int INVALID_ID{std::numeric_limits<unsigned int>::max()};

    CinemaImagingEmbree(){
        this->setDebugMsgPrefix("CinemaImaging(Embree)");
    }
    ~CinemaImagingEmbree(){};

    int initializeDevice(
      RTCDevice& device
    ) const;

    template<typename IT>
    int initializeScene(
      RTCScene& scene,

      const RTCDevice& device,
      const size_t& nVertices,
      const float* vertexCoords,
      const size_t& nTriangles,
      const IT* connectivityList
    ) const;

    int renderImage(
      float* depthBuffer,
      unsigned int* primitiveIds,
      float* barycentricCoordinates,

      const RTCScene& scene,
      const double resolution[2],
      const double camCenter[3],
      const double camDir[3],
      const double camUp[3],
      const double& camHeight,
      const bool& orthographicProjection = true,
      const double& focalLength = 1
    ) const;

    template<typename DT, typename IT>
    int interpolateArray(
      DT* outputArray,

      const unsigned int* primitiveIds,
      const float* barycentricCoordinates,
      const IT* connectivityList,

      const DT* inputArray,
      const size_t& nTuples,
      const size_t& nComponents=1,
      const DT& missingValue= std::numeric_limits<DT>::has_quiet_NaN ? std::numeric_limits<DT>::quiet_NaN() : std::numeric_limits<DT>::max()
    ) const;

    template<typename DT>
    int lookupArray(
      DT* outputArray,

      const unsigned int* primitiveIds,

      const DT* inputArray,
      const size_t& nTuples,
      const size_t& nComponents=1,
      const DT& missingValue= std::numeric_limits<DT>::has_quiet_NaN ? std::numeric_limits<DT>::quiet_NaN() : std::numeric_limits<DT>::max()
    ) const;

    int deallocateScene(
      RTCDevice& device,
      RTCScene& scene
    ) const;
  };
};

template<typename DT, typename IT>
int ttk::CinemaImagingEmbree::interpolateArray(
  DT* outputArray,

  const unsigned int* primitiveIds,
  const float* barycentricCoordinates,
  const IT* connectivityList,

  const DT* inputArray,
  const size_t& nTuples,
  const size_t& nComponents,
  const DT& missingValue
) const {

  if(nComponents!=1){
    this->printErr("Current implementation only supports interpolation of scalars.");
    return 0;
  }

  #ifdef TTK_ENABLE_OPENMP
  #pragma omp parallel for num_threads(this->threadNumber_)
  #endif
  for(size_t i=0; i<nTuples; i++){
    const unsigned int& cellId = primitiveIds[i];
    if(cellId==CinemaImagingEmbree::INVALID_ID){
      outputArray[i] = missingValue;
      continue;
    }

    const size_t cellIndex = cellId*3;
    const IT& v0 = connectivityList[cellIndex+0];
    const IT& v1 = connectivityList[cellIndex+1];
    const IT& v2 = connectivityList[cellIndex+2];

    const size_t bcIndex = i*2;
    const float& u = barycentricCoordinates[bcIndex+0];
    const float& v = barycentricCoordinates[bcIndex+1];
    const float w = 1 - u - v;

    outputArray[i] = w*inputArray[v0] + u*inputArray[v1] + v*inputArray[v2];
  }

  return 1;
};

template<typename DT>
int ttk::CinemaImagingEmbree::lookupArray(
  DT* outputArray,

  const unsigned int* primitiveIds,

  const DT* inputArray,
  const size_t& nTuples,
  const size_t& nComponents,
  const DT& missingValue
) const {

  for(size_t i=0; i<nTuples; i++){
    size_t outputOffset = i*nComponents;
    const unsigned int& cellId = primitiveIds[i];

    if(cellId==CinemaImagingEmbree::INVALID_ID){
      for(size_t j=0; j<nComponents; j++)
        outputArray[ outputOffset++ ] = missingValue;
      continue;
    } else {
      size_t inputOffset = cellId*nComponents;
      for(size_t j=0; j<nComponents; j++)
        outputArray[ outputOffset++ ] = inputArray[ inputOffset++ ];
    }
  }

  return 1;
};

int ttk::CinemaImagingEmbree::deallocateScene(
  RTCDevice& device,
  RTCScene& scene
) const {
  ttk::Timer timer;
  this->printMsg("Deallocating Scene",0,0,ttk::debug::LineMode::REPLACE);

  rtcReleaseScene(scene);
  rtcReleaseDevice(device);

  this->printMsg("Deallocating Scene",1,timer.getElapsedTime());

  return 1;
};

int ttk::CinemaImagingEmbree::initializeDevice(RTCDevice& device) const {
  ttk::Timer timer;
  this->printMsg("Initializing Device",0,0,ttk::debug::LineMode::REPLACE);

  device = rtcNewDevice("hugepages=1,threads=1");

  if(!device){
    this->printErr("Unable to create device");
    this->printErr( std::to_string(rtcGetDeviceError(NULL)) );
    return 0;
  }

  auto errorFunction = [](void* userPtr, enum RTCError error, const char* str){
      printf("error %d: %s\n", error, str);
  };

  rtcSetDeviceErrorFunction(device, errorFunction, NULL);

  this->printMsg("Initializing Device",1,timer.getElapsedTime());

  return 1;
};

template<typename IT>
int ttk::CinemaImagingEmbree::initializeScene(
  RTCScene& scene,

  const RTCDevice& device,
  const size_t& nVertices,
  const float* vertexCoords,
  const size_t& nTriangles,
  const IT* connectivityList
) const {
  ttk::Timer timer;
  this->printMsg("Initializing Scene (#v:"+std::to_string(nVertices)+"|#t:"+std::to_string(nTriangles)+")",0,0,ttk::debug::LineMode::REPLACE);

    scene = rtcNewScene(device);

    RTCGeometry mesh = rtcNewGeometry(
        device,
        RTC_GEOMETRY_TYPE_TRIANGLE
    );

    // vertices
    {
      rtcSetSharedGeometryBuffer(
        mesh,
        RTC_BUFFER_TYPE_VERTEX,
        0,
        RTC_FORMAT_FLOAT3,
        (const void*) vertexCoords,
        0,
        3*sizeof(float),
        nVertices
      );
    }

    // triangles
    {
        // unfortunately embree does not support signed integer based indexing
    //   rtcSetSharedGeometryBuffer(
    //       mesh,
    //       RTC_BUFFER_TYPE_INDEX,
    //       0,
    //       RTC_FORMAT_LLONG3,
    //       static_cast<const void*>(connectivityList),
    //       0,
    //       3*sizeof(long long),
    //       nTriangles
    //   );

      unsigned int* indices = (unsigned int*) rtcSetNewGeometryBuffer(
        mesh,
        RTC_BUFFER_TYPE_INDEX,
        0,
        RTC_FORMAT_UINT3,
        3*sizeof(unsigned int),
        nTriangles
      );

      for(size_t t=0, tn=nTriangles*3; t<tn; t++)
        indices[t] = connectivityList[t];

      rtcCommitGeometry(mesh);
      rtcAttachGeometry(scene, mesh);
      rtcReleaseGeometry(mesh);
    }

    rtcCommitScene(scene);

    this->printMsg("Initializing Scene (#v:"+std::to_string(nVertices)+"|#t:"+std::to_string(nTriangles)+")",1,timer.getElapsedTime());
    return 1;
};

int ttk::CinemaImagingEmbree::renderImage(
  float* depthBuffer,
  unsigned int* primitiveIds,
  float* barycentricCoordinates,

  const RTCScene& scene,
  const double resolution[2],
  const double camPos[3],
  const double camDirRaw[3],
  const double camUp[3],
  const double& camHeight,
  const bool& orthographicProjection,
  const double& viewAngle
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

  struct RTCRayHit rayhit;
  rayhit.ray.dir_x = -1;
  rayhit.ray.dir_y = -1;
  rayhit.ray.dir_z = -1;

  size_t pixelIndex = 0;
  size_t bcIndex = 0;
  float nan = std::numeric_limits<float>::quiet_NaN();


  int resX = resolution[0];
  int resY = resolution[1];
  if(orthographicProjection){

    for(int y = 0; y < resY; y++) {
        double v = ((double)y) * pixelHeightWorld;

        for(int x = 0; x < resX; x++) {
          double u = ((double)x) * pixelWidthWorld;

          // set origin
          rayhit.ray.org_x = camPosCorner[0] + u * camRight[0] + v * camUpTrue[0];
          rayhit.ray.org_y = camPosCorner[1] + u * camRight[1] + v * camUpTrue[1];
          rayhit.ray.org_z = camPosCorner[2] + u * camRight[2] + v * camUpTrue[2];

          // set dir
          rayhit.ray.dir_x = camDir[0];
          rayhit.ray.dir_y = camDir[1];
          rayhit.ray.dir_z = camDir[2];

          // compute hit
          rayhit.ray.tnear = 0.01;
          rayhit.ray.tfar = INFINITY;
          rayhit.ray.mask = 0;
          rayhit.ray.flags = 0;
          rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
          rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
          rtcIntersect1(scene, &context, &rayhit);

          // write depth
          bool hitPrimitive = rayhit.hit.geomID!=RTC_INVALID_GEOMETRY_ID;
          if(hitPrimitive){
            depthBuffer[pixelIndex] = std::max(0.0f,rayhit.ray.tfar);
            primitiveIds[pixelIndex] = rayhit.hit.primID;
            barycentricCoordinates[bcIndex++] = rayhit.hit.u;
            barycentricCoordinates[bcIndex++] = rayhit.hit.v;
          } else {
            depthBuffer[pixelIndex] = nan;
            primitiveIds[pixelIndex] = CinemaImagingEmbree::INVALID_ID;
            barycentricCoordinates[bcIndex++] = nan;
            barycentricCoordinates[bcIndex++] = nan;
          }
          pixelIndex++;
        }
      }
  } else {

      double focalLength = camSize[0]/2 / tan(viewAngle/360.0 * 3.141592653589793);

      for(int y = 0; y < resY; y++) {
        double v = ((double)y) * pixelHeightWorld;

        for(int x = 0; x < resX; x++) {
          double u = ((double)x) * pixelWidthWorld;

          // set origin
          rayhit.ray.org_x = camPos[0];
          rayhit.ray.org_y = camPos[1];
          rayhit.ray.org_z = camPos[2];

          // set dir
          rayhit.ray.dir_x = camPosCorner[0] + u * camRight[0] + v * camUpTrue[0] + camDir[0]*focalLength - rayhit.ray.org_x;
          rayhit.ray.dir_y = camPosCorner[1] + u * camRight[1] + v * camUpTrue[1] + camDir[1]*focalLength - rayhit.ray.org_y;
          rayhit.ray.dir_z = camPosCorner[2] + u * camRight[2] + v * camUpTrue[2] + camDir[2]*focalLength - rayhit.ray.org_z;

          // compute hit
          rayhit.ray.tnear = 0.01;
          rayhit.ray.tfar = INFINITY;
          rayhit.ray.mask = 0;
          rayhit.ray.flags = 0;
          rayhit.hit.geomID = RTC_INVALID_GEOMETRY_ID;
          rayhit.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
          rtcIntersect1(scene, &context, &rayhit);

          // write depth
          bool hitPrimitive = rayhit.hit.geomID!=RTC_INVALID_GEOMETRY_ID;
          if(hitPrimitive){
            depthBuffer[pixelIndex] = std::max(0.0f,rayhit.ray.tfar);
            primitiveIds[pixelIndex] = rayhit.hit.primID;
            barycentricCoordinates[bcIndex++] = rayhit.hit.u;
            barycentricCoordinates[bcIndex++] = rayhit.hit.v;
          } else {
            depthBuffer[pixelIndex] = nan;
            primitiveIds[pixelIndex] = CinemaImagingEmbree::INVALID_ID;
            barycentricCoordinates[bcIndex++] = nan;
            barycentricCoordinates[bcIndex++] = nan;
          }
          pixelIndex++;
        }
      }
  }

  this->printMsg(
      "Rendering Image ("+std::string(orthographicProjection ? "O" : "P")+ "|"+std::to_string(resolution[0])+"x"+std::to_string(resolution[1])+")",
      1,timer.getElapsedTime(),this->threadNumber_
  );

  return 1;
};

#endif