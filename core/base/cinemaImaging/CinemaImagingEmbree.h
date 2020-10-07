/// \ingroup base
/// \class ttk::CinemaImagingEmbree
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.05.2020
///
/// \brief TTK %CinemaImagingEmbree processing package.
///
/// %CinemaImagingEmbree is a TTK processing package that

#pragma once

#include <Debug.h>

#if TTK_ENABLE_EMBREE
#include <string>
#include <embree3/rtcore.h>
#include <limits>
#endif

namespace ttk {

  class CinemaImagingEmbree : virtual public Debug {
  public:

    CinemaImagingEmbree();
    ~CinemaImagingEmbree();

#if TTK_ENABLE_EMBREE
    static const unsigned int INVALID_ID{std::numeric_limits<unsigned int>::max()};

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

#endif
  };
};


#if TTK_ENABLE_EMBREE

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
        indices[t] = (unsigned int)connectivityList[t];

      rtcCommitGeometry(mesh);
      rtcAttachGeometry(scene, mesh);
      rtcReleaseGeometry(mesh);
    }

    rtcCommitScene(scene);

    this->printMsg("Initializing Scene (#v:"+std::to_string(nVertices)+"|#t:"+std::to_string(nTriangles)+")",1,timer.getElapsedTime());
    return 1;
};

#endif