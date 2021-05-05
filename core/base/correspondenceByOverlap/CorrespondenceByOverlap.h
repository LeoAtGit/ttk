/// TODO 1: Provide your information
///
/// \ingroup base
/// \class ttk::CorrespondenceByOverlap
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// This module defines the %CorrespondenceByOverlap class that computes for each
/// vertex of a triangulation the average scalar value of itself and its direct
/// neighbors.
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

    typedef std::unordered_map<long long,long long> LabelIndexMap;

    CorrespondenceByOverlap() {
      this->setDebugMsgPrefix("CorrespondenceByOverlap");
    };
    ~CorrespondenceByOverlap(){};

    template <class DT>
    int computeLabelIndexMap(
      LabelIndexMap& labelIndexMap,
      const DT* labels,
      const int nLabels
    ) const {

      long long labelIndex = 0;
      for(int i=0; i<nLabels; i++){
        const auto& l = labels[i];
        if(l>=0 && labelIndexMap.find(l)==labelIndexMap.end())
          labelIndexMap.insert({l,labelIndex++});
      }

      return 1;
    };

    template <class DT>
    int computeAdjacencyMatrix(int *adjacencyMatrix,
                               const DT *labels0,
                               const DT *labels1,
                               const int nVertices,
                               const LabelIndexMap& labelIndexMap0,
                               const LabelIndexMap& labelIndexMap1
    ) const {

      ttk::Timer timer;

      const auto nLabels0 = labelIndexMap0.size();
      const auto nLabels1 = labelIndexMap1.size();

      if(nLabels0 < 1)
        return this->printWrn("Number of first labels smaller than 1.");
      if(nLabels1 < 1)
        return this->printWrn("Number of second labels smaller than 1.");

      const std::string msg = "Computing Overlap "+std::to_string(nLabels0)+"x"+std::to_string(nLabels1);
      this->printMsg(msg, 0, 0, this->threadNumber_, debug::LineMode::REPLACE);

      int nLables = nLabels0 * nLabels1;

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
          adjacencyMatrix[labelIndexMap1.at(l1) * nLabels0 + labelIndexMap0.at(l0)]++;
        }
      }

      this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);

      return 1; // return success
    }
  };
} // namespace ttk
