/// TODO 1: Provide your information
///
/// \ingroup base
/// \class ttk::TrackingFromOverlap
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// This module defines the %TrackingFromOverlap class that computes for each
/// vertex of a triangulation the average scalar value of itself and its direct
/// neighbors.
///
/// \b Related \b publication: \n
/// 'TrackingFromOverlap'
/// Jonas Lukasczyk and Julien Tierny.
/// TTK Publications.
/// 2020.
///

#pragma once

#include <Debug.h>

namespace ttk {

  class TrackingFromOverlap : virtual public Debug {

  public:
    TrackingFromOverlap() {
      this->setDebugMsgPrefix("TrackingFromOverlap");
    };
    ~TrackingFromOverlap(){};

    template <class DT>
    int computeAdjacencyMatrix(int* adjacencyMatrix,
                     const DT *labels0,
                     const DT *labels1,
                     const int nVertices,
                     const int nLabels0,
                     const int nLabels1) const {

      ttk::Timer timer;

      if(nLabels0<1)
        return this->printWrn("Number of first labels smaller than 1.");
      if(nLabels1<1)
        return this->printWrn("Number of second labels smaller than 1.");

      this->printMsg("Computing Overlap", 0, 0, this->threadNumber_,
                     ttk::debug::LineMode::REPLACE);

      int nLables = nLabels0*nLabels1;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
      for(int i = 0; i < nLables; i++) {
        adjacencyMatrix[i] = 0;
      }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
      for(int i = 0; i < nVertices; i++) {
        const int l0 = labels0[i];
        const int l1 = labels1[i];
        if(l0 >= 0 && l1 >= 0) {
#pragma omp atomic update
          adjacencyMatrix[l1 * nLabels0 + l0]++;
        }
      }

      this->printMsg(
        "Computing Overlap", 1, timer.getElapsedTime(), this->threadNumber_);

      return 1; // return success
    }
  };
} // namespace ttk
