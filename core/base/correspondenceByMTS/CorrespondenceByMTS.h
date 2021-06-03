#pragma once

#include <Debug.h>

namespace ttk {

  class CorrespondenceByMTS : virtual public Debug {

  public:
    CorrespondenceByMTS() {
      this->setDebugMsgPrefix("CorrespondenceByMTS");
    };
    ~CorrespondenceByMTS(){};

    template <typename IT, typename DT>
    int getBaseRepresentative(IT &base,

                              const DT baselevel,
                              const DT *scalars,
                              const IT *next) const {
      IT nextN = next[base];
      while(nextN >= 0 && scalars[nextN] > baselevel) {
        base = nextN;
        nextN = next[base];
      }

      return 1;
    }

    template <typename IT, typename DT>
    int computeSegmentationOverlap(IT *correlationsArray,

                                   const IT *seg0,
                                   const IT *seg1,
                                   const IT nVertices,
                                   const IT *next0,
                                   const IT *next1,
                                   const DT *scalars0,
                                   const DT *scalars1,
                                   const IT nEdges0,
                                   const IT nEdges1) const {

      ttk::Timer timer;

      const std::string msg
        = "Computing Segmentation Overlap (" + std::to_string(nVertices) + ")";
      this->printMsg(
        msg, 0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);

      // clear
      const IT nEdgePairs = nEdges0 * nEdges1;
      for(IT i = 0; i < nEdgePairs; i++) {
        correlationsArray[i] = 0;
      }

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(this->threadNumber_)
#endif
      for(IT v = 0; v < nVertices; v++) {
        const auto &n0 = seg0[v];
        const auto &n1 = seg1[v];

        const auto &bl = scalars0[next0[n0]] < scalars1[next1[n1]]
                           ? scalars0[next0[n0]]
                           : scalars1[next1[n1]];

        IT base0 = n0;
        IT base1 = n1;
        this->getBaseRepresentative<IT, DT>(base0, bl, scalars0, next0);
        this->getBaseRepresentative<IT, DT>(base1, bl, scalars1, next1);

        correlationsArray[base1 * nEdges0 + base0]++;
//         {
//           IT c0 = base0;
//           IT c1 = base1;

//           while(true) {
// #ifdef TTK_ENABLE_OPENMP
// #pragma omp atomic update
// #endif
//             correlationsArray[c1 * nEdges0 + c0]++;

//             const auto &cNext0 = next0[c0];
//             const auto &cNext1 = next1[c1];

//             if(cNext0 < 0 && cNext1 < 0)
//               break;

//             if(cNext0 < 0) {
//               c1 = cNext1;
//             } else if(cNext1 < 0) {
//               c0 = cNext0;
//             } else {
//               const auto &cNextScalar0 = scalars0[cNext0];
//               const auto &cNextScalar1 = scalars1[cNext1];

//               if(cNextScalar0 == cNextScalar1) {
//                 c0 = cNext0;
//                 c1 = cNext1;
//               } else if(cNextScalar0 < cNextScalar1) {
//                 c1 = cNext1;
//               } else {
//                 c0 = cNext0;
//               }
//             }
//           };
//         }
      }

      this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);

      return 1;
    }
  };
} // namespace ttk
