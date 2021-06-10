///
/// \ingroup base
/// \class ttk::CorrespondenceByPersistencePairs
/// \author Maxime Soler
/// \date 09.06.2021
///
/// \brief Computes a correspondence matrix from persistence diagram matchings.
/// 
/// This module defines the CorrespondenceByPersistencePairs class that computes 
/// the correspondance matrix from the output of a Wasserstein-based matching
/// between persistence diagrams. 1 row = 1 persistence pair in diagram 1, 
/// 1 column = 1 persistence pair in diagram 2, 
/// value = 0 => no matching
/// value > 0 => matching.
///
/// Maxime Soler, Jonas Lukasczyk, Julien Tierny, 2020.
///

#pragma once

#include <Debug.h>
#include <BottleneckDistance.h>

namespace ttk {

  class CorrespondenceByPersistencePairs : virtual public Debug {

  public:
    CorrespondenceByPersistencePairs() {
      this->setDebugMsgPrefix("CorrespondenceByPersistencePairs");
    };
    ~CorrespondenceByPersistencePairs(){};

    template <class DT>
    int computeDistanceMatrix(
        //DT *distanceMatrix,
        //const DT *coords0,
        std::vector< std::tuple<int, ttk::CriticalType, int, ttk::CriticalType, DT,
            int, DT, float, float, float, DT, float, float, float> > &CTDiagram0,
        std::vector< std::tuple<int, ttk::CriticalType, int, ttk::CriticalType, DT,
            int, DT, float, float, float, DT, float, float, float> > &CTDiagram1,
        std::vector< std::tuple<int, int, double> > &matchings,
        //const DT *coords1,
        //const int nFeatures0,
        //const int nFeatures1,
        const double PX, const double PY, const double PZ, 
        const double PS, const double PE,
        const std::string algorithm,
        const std::string wasserstein,
        const double alpha,
        const int pvAlgorithm
        ) const
    {
      ttk::Timer timer;

      //const int nPoints0 = nPairs0 * 2;
      //const int nPoints1 = nPairs1 * 2;
      //const std::string msg = "Computing Distance Matrix ("
      //                        + std::to_string(nPoints0) + "x"
      //                        + std::to_string(nPoints1) + ")";
      //this->printMsg(
      //  msg, 0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);

      ttk::BottleneckDistance bottleneckDistance_;
      bottleneckDistance_.setPersistencePercentThreshold(0);
      // (tolerance can be set earlier, when features are defined)
      bottleneckDistance_.setPX(PX);
      bottleneckDistance_.setPY(PY);
      bottleneckDistance_.setPZ(PZ);
      bottleneckDistance_.setPS(PS);
      bottleneckDistance_.setPE(PE);
      bottleneckDistance_.setAlgorithm(algorithm);
      bottleneckDistance_.setPVAlgorithm(pvAlgorithm);
      bottleneckDistance_.setWasserstein(wasserstein);

      bottleneckDistance_.setCTDiagram1(&CTDiagram0);
      bottleneckDistance_.setCTDiagram2(&CTDiagram1);
      bottleneckDistance_.setOutputMatchings(&matchings);
      int status = bottleneckDistance_.execute<double>(false);
      if (status < 0) return -1;

//#ifdef TTK_ENABLE_OPENMP
//#pragma omp parallel for num_threads(this->threadNumber_)
//#endif
      //for(int i = 0; i < nPoints0; i++) {
      //  for(int j = 0; j < nPoints1; j++) {
      //    const int i3 = i * 3;
      //    const int j3 = j * 3;
      //    const DT dx = coords0[i3 + 0] - coords1[j3 + 0];
      //    const DT dy = coords0[i3 + 1] - coords1[j3 + 1];
      //    const DT dz = coords0[i3 + 2] - coords1[j3 + 2];
      //    distanceMatrix[j * nPoints0 + i]
      //      = std::sqrt(dx * dx + dy * dy + dz * dz);
      //  }
      //}

      //this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);

      return 1;
    }
  };
} // namespace ttk
