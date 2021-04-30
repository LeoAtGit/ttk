/// TODO 1: Provide your information
///
/// \ingroup base
/// \class ttk::CalculatePorosity
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// This module defines the %CalculatePorosity class that computes for each
/// vertex of a triangulation the average scalar value of itself and its direct
/// neighbors.
///
/// \b Related \b publication: \n
/// 'CalculatePorosity'
/// Jonas Lukasczyk and Julien Tierny.
/// TTK Publications.
/// 2021.
///

#pragma once

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>
#include <set>
#include <unordered_map>
#include <cmath>

namespace ttk {

  /**
  * The CalculatePorosity class provides methods to compute for each vertex of
  * a triangulation the average scalar value of itself and its direct
  * neighbors.
  */
  class CalculatePorosity: virtual public Debug {

    public:
    CalculatePorosity() {
      this->setDebugMsgPrefix(
        "CalculatePorosity"); // inherited from Debug: prefix will be printed at
      // the
      // beginning of every msg
    };
    ~CalculatePorosity() {};
    // /**
    // * TODO 2: This method preconditions the triangulation for all operations
    // *         the algorithm of this module requires. For instance,
    // *         preconditionVertexNeighbors, preconditionBoundaryEdges, ...
    // *
    // *         Note: If the algorithm does not require a triangulation then
    // *               this method can be deleted.
    // */
    // int preconditionTriangulation(
    //   ttk::AbstractTriangulation *triangulation) const {
    //   return triangulation->preconditionVertexNeighbors();
    // };

    /**
    * TODO 3: Implmentation of the algorithm.
    *
    *         Note: If the algorithm requires a triangulation then this
    *               method must be called after the triangulation has been
    *               preconditioned for the upcoming operations.
    */
    template <typename dataType>
    int computePorosity(float *outputData,
      const size_t& nVertices,
      const dataType *inputData,
      const float *gradientData,
      const float *divergence,
      const float &distance,
      const float &threshold,
      const float &margin,
      const float &maxThreshold,
      const float &gradientThreshold,
      const int* dim) const {
      // start global timer
      ttk::Timer globalTimer;
      // print horizontal separator
      this->printMsg(ttk::debug::Separator::L1); // L1 is the '=' separator
      // print input parameters in table format
      this->printMsg({
        {"#Threads", std::to_string(this->threadNumber_)},
        {"#Vertices", std::to_string(nVertices)},
      });
      this->printMsg(ttk::debug::Separator::L1);

      {
        // start a local timer for this subprocedure
        ttk::Timer localTimer;

        // print the progress of the current subprocedure (currently 0%)
        this->printMsg("Computing Porosity",
          0, // progress form 0-1
          0, // elapsed time so far
          this->threadNumber_, ttk::debug::LineMode::REPLACE);

	int dimCount = 0;
	for(int i = 0; i < 3; i++) {
	  if(dim[i] > 1) {
	    dimCount++;
	  }
	}

	float sqrDistance = (dimCount == 0) ? 1.0f : 0.0f;
	
	for(int i = 0; i < dimCount; i++) {
	  sqrDistance += distance * distance;
	}
	
        const float marginThreshold = (1.0f + margin) * threshold;
        const float fractionalThreshold = margin * threshold;

	int dimX = dim[0];
        int dimY = dim[1];
        int dimZ = dim[2];

        // init masks
        std::vector < std::tuple < int,
        int,
        int,
        float>> neighbor_mask;

        for(int i = -int(distance); i < int(distance); i++) {
          for(int j = -int(distance); j < int(distance); j++) {
            for(int k = -int(distance); k < int(distance); k++) {
	      float d = i * i + j * j + k * k;
	      //if(d > 0.0f) { 
		neighbor_mask.push_back(std::tuple < int, int, int, float > (i, j, k, 1.0f - (d / sqrDistance)));
		//}
            }
          }
        }


        #ifdef TTK_ENABLE_OPENMP
        #pragma omp parallel for num_threads(this->threadNumber_)
        #endif
        for(int i = 0; i < dimX; i++) {
          for(int j = 0; j < dimY; j++) {
            for(int k = 0; k < dimZ; k++) {
              int idx = k * dimX * dimY + j * dimX + i;

              outputData[idx] = 0;

              if(inputData[idx] > maxThreshold) {
                outputData[idx] = -1;
                continue;
              }

              //if(gradientData[idx] > gradientThreshold)
	      //continue;

              for(auto& it: neighbor_mask) {
                int x = i + std::get < 0 > (it);
                int y = j + std::get < 1 > (it);
                int z = k + std::get < 2 > (it);

                if(x < 0 || x >= dimX || y < 0 || y >= dimY || z < 0 || z >= dimZ) continue;
                int idxN = z * dimX * dimY + y * dimX + x;

                // do the threshold check here
                // marginThreshold is (1 + margin) * threshold - basically the larger
                // threshold fractionalThreshold is threshold * margin that is the
                // maximum distance allowed to be counted as a fraction
                if(inputData[idxN] < threshold) {
                  outputData[idx] += inputData[idxN] * std::get < 3 > (it);// * divergence[idxN];
                } else if(inputData[idxN] < marginThreshold) {
                  outputData[idx] += ((marginThreshold - inputData[idxN]) / (fractionalThreshold)) * std::get <3> (it);// * divergence[idxN];
                }
              }
            }
          }
        }
        // print the progress of the current subprocedure with elapsed time
        this->printMsg("Computing Porosity",
          1, // progress
          localTimer.getElapsedTime(), this->threadNumber_);
      }

      // ---------------------------------------------------------------------
      // print global performance
      // ---------------------------------------------------------------------

      {
        this->printMsg(ttk::debug::Separator::L2); // horizontal '-' separator
        this->printMsg(
          "Complete", 1, globalTimer.getElapsedTime() // global progress, time
        );
        this->printMsg(ttk::debug::Separator::L1); // horizontal '=' separator
      }

      return 1; // return success
    }
  }; // CalculatePorosity class
} // namespace ttk
