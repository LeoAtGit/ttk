/// TODO 1: Provide your information
///
/// \ingroup base
/// \class ttk::CorrespondenceByDistance
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// This module defines the %CorrespondenceByDistance class that computes for
/// each vertex of a triangulation the average scalar value of itself and its
/// direct neighbors.
///
/// \b Related \b publication: \n
/// 'CorrespondenceByDistance'
/// Jonas Lukasczyk and Julien Tierny.
/// TTK Publications.
/// 2020.
///

#pragma once

#include <Debug.h>

namespace ttk {

  class CorrespondenceByDistance : virtual public Debug {

  public:
    CorrespondenceByDistance() {
      this->setDebugMsgPrefix("CorrespondenceByDistance");
    };
    ~CorrespondenceByDistance(){};

    template <class DT>
    int computeDistanceMatrix(DT *distanceMatrix,
                              const DT *coords0,
                              const DT *coords1,
                              const int nPoints0,
                              const int nPoints1) const {

      ttk::Timer timer;

      const std::string msg = "Computing Distance Matrix ("
                              + std::to_string(nPoints0) + "x"
                              + std::to_string(nPoints1) + ")";
      this->printMsg(
        msg, 0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
      for(int i = 0; i < nPoints0; i++) {
        for(int j = 0; j < nPoints1; j++) {

          const int i3 = i * 3;
          const int j3 = j * 3;

          const DT dx = coords0[i3 + 0] - coords1[j3 + 0];
          const DT dy = coords0[i3 + 1] - coords1[j3 + 1];
          const DT dz = coords0[i3 + 2] - coords1[j3 + 2];

          distanceMatrix[j * nPoints0 + i]
            = std::sqrt(dx * dx + dy * dy + dz * dz);
        }
      }

      this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);

      return 1;
    }
  };
} // namespace ttk
