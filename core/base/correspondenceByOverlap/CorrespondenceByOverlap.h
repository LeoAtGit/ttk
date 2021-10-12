/// TODO 1: Provide your information
///
/// \ingroup base
/// \class ttk::CorrespondenceByOverlap
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// This module defines the %CorrespondenceByOverlap class that computes for
/// each vertex of a triangulation the average scalar value of itself and its
/// direct neighbors.
///
/// \b Related \b publication: \n
/// 'CorrespondenceByOverlap'
/// Jonas Lukasczyk and Julien Tierny.
/// TTK Publications.
/// 2020.
///

#pragma once

#include <Debug.h>

#include <unordered_map>

namespace ttk {

  class CorrespondenceByOverlap : virtual public Debug {

  public:
    CorrespondenceByOverlap() {
      this->setDebugMsgPrefix("CorrespondenceByOverlap");
    };
    ~CorrespondenceByOverlap(){};

    template <typename DT, typename IT>
    int computeLabelIndexMap(std::unordered_map<IT, IT> &labelIndexMap,
                             const DT *labels,
                             const IT nLabels) const {

      IT labelIndex = 0;
      for(IT i = 0; i < nLabels; i++) {
        auto l = static_cast<const IT>(labels[i]);
        if(l >= 0 && labelIndexMap.find(l) == labelIndexMap.end())
          labelIndexMap.insert({l, labelIndex++});
      }

      return 1;
    };

    template <typename DT, typename IT>
    int computeAdjacencyMatrix(
      int *adjacencyMatrix,
      const DT *labels0,
      const DT *labels1,
      const int nVertices,
      const std::unordered_map<IT, IT> &labelIndexMap0,
      const std::unordered_map<IT, IT> &labelIndexMap1) const {

      ttk::Timer timer;

      const IT nLabels0 = labelIndexMap0.size();
      const IT nLabels1 = labelIndexMap1.size();

      if(nLabels0 < 1)
        return this->printWrn("Number of first labels smaller than 1.");
      if(nLabels1 < 1)
        return this->printWrn("Number of second labels smaller than 1.");

      const std::string msg = "Computing Overlap " + std::to_string(nLabels0)
                              + "x" + std::to_string(nLabels1);
      this->printMsg(msg, 0, 0, this->threadNumber_, debug::LineMode::REPLACE);

      const IT nLables = nLabels0 * nLabels1;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
      for(IT i = 0; i < nLables; i++) {
        adjacencyMatrix[i] = 0;
      }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
      for(IT i = 0; i < nVertices; i++) {
        auto l0 = static_cast<const IT>(labels0[i]);
        auto l1 = static_cast<const IT>(labels1[i]);
        if(l0 >= 0 && l1 >= 0) {
#pragma omp atomic update
          adjacencyMatrix[labelIndexMap1.at(l1) * nLabels0
                          + labelIndexMap0.at(l0)]++;
        }
      }

      this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);

      return 1; // return success
    }
  };
} // namespace ttk
