#pragma once

#include <Debug.h>

#include <limits>
#include <vector>

namespace ttk {
  class CinemaImaging : virtual public Debug {
  public:
    static const unsigned int INVALID_ID{
      std::numeric_limits<unsigned int>::max()};
    template <typename DT, typename IT>
    int interpolateArray(float *outputArray,

                         const unsigned int *primitiveIds,
                         const float *barycentricCoordinates,
                         const IT *connectivityList,

                         const DT *inputArray,
                         const size_t &nPixels,
                         const size_t &nComponents = 1) const;

    template <typename DT>
    int lookupArray(float *outputArray,

                    const unsigned int *primitiveIds,

                    const DT *inputArray,
                    const size_t &nPixels,
                    const size_t &nComponents = 1) const;
  };
} // namespace ttk

template <typename DT, typename IT>
int ttk::CinemaImaging::interpolateArray(float *outputArray,

                                         const unsigned int *primitiveIds,
                                         const float *barycentricCoordinates,
                                         const IT *connectivityList,

                                         const DT *inputArray,
                                         const size_t &nPixels,
                                         const size_t &nComponents) const {

  const auto missingValue = std::numeric_limits<float>::quiet_NaN();

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
  for(size_t i = 0; i < nPixels; i++) {
    const unsigned int &cellId = primitiveIds[i];
    const auto offset = i * nComponents;
    if(cellId == CinemaImaging::INVALID_ID) {
      for(size_t c = 0; c < nComponents; c++)
        outputArray[offset + c] = missingValue;
      continue;
    }

    const size_t cellIndex = cellId * 3;
    const auto o0 = nComponents * connectivityList[cellIndex + 0];
    const auto o1 = nComponents * connectivityList[cellIndex + 1];
    const auto o2 = nComponents * connectivityList[cellIndex + 2];

    const size_t bcIndex = i * 2;
    const float &u = barycentricCoordinates[bcIndex + 0];
    const float &v = barycentricCoordinates[bcIndex + 1];
    const float w = 1 - u - v;

    for(size_t c = 0; c < nComponents; c++)
      outputArray[offset + c] = w * static_cast<float>(inputArray[o0 + c])
                                + u * static_cast<float>(inputArray[o1 + c])
                                + v * static_cast<float>(inputArray[o2 + c]);
  }

  return 1;
};

template <typename DT>
int ttk::CinemaImaging::lookupArray(float *outputArray,

                                    const unsigned int *primitiveIds,

                                    const DT *inputArray,
                                    const size_t &nPixels,
                                    const size_t &nComponents) const {

  const auto missingValue = std::numeric_limits<float>::quiet_NaN();

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
  for(size_t i = 0; i < nPixels; i++) {
    size_t outputOffset = i * nComponents;
    const unsigned int &cellId = primitiveIds[i];

    if(cellId == CinemaImaging::INVALID_ID) {
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
