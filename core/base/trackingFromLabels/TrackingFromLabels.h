/// TODO 1: Provide your information
///
/// \ingroup base
/// \class ttk::TrackingFromLabels
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// This module defines the %TrackingFromLabels class that computes for each
/// vertex of a triangulation the average scalar value of itself and its direct
/// neighbors.
///
/// \b Related \b publication: \n
/// 'TrackingFromLabels'
/// Jonas Lukasczyk and Julien Tierny.
/// TTK Publications.
/// 2020.
///

#pragma once

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>

namespace ttk {

  class TrackingFromLabels : virtual public Debug {

  public:
    TrackingFromLabels() {
      this->setDebugMsgPrefix("TrackingFromLabels");
    };
    ~TrackingFromLabels(){};

    template <class DT>
    int computeEdges(std::vector<int> &adjacencyMatrix,
                     const DT *labels0,
                     const DT *labels1,
                     const size_t nVertices,
                     const size_t nNodes0,
                     const size_t nNodes1) const {

      ttk::Timer timer;
      this->printMsg("Computing Overlap", 0, 0, this->threadNumber_,
                     ttk::debug::LineMode::REPLACE);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
      for(size_t i = 0; i < nVertices; i++) {
        const auto &l0 = labels0[i];
        const auto &l1 = labels1[i];
        if(l0 >= 0 && l1 >= 0) {
#pragma omp atomic update
          adjacencyMatrix[l1 * nNodes0 + l0]++;
        }
      }

      this->printMsg(
        "Computing Overlap", 1, timer.getElapsedTime(), this->threadNumber_);

      return 1; // return success
    }
  };
} // namespace ttk
