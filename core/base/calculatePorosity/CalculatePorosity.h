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
#include <unordered_map>
#include <set>
#include <valarray>

namespace ttk {

  /**
   * The CalculatePorosity class provides methods to compute for each vertex of
   * a triangulation the average scalar value of itself and its direct
   * neighbors.
   */
  class CalculatePorosity : virtual public Debug {

  public:
    CalculatePorosity() {
      this->setDebugMsgPrefix(
        "CalculatePorosity"); // inherited from Debug: prefix will be printed at
                              // the
      // beginning of every msg
    };
    ~CalculatePorosity(){};
    /**
     * TODO 2: This method preconditions the triangulation for all operations
     *         the algorithm of this module requires. For instance,
     *         preconditionVertexNeighbors, preconditionBoundaryEdges, ...
     *
     *         Note: If the algorithm does not require a triangulation then
     *               this method can be deleted.
     */
    int preconditionTriangulation(
      ttk::AbstractTriangulation *triangulation) const {
      return triangulation->preconditionVertexNeighbors();
    };

    /**
     * TODO 3: Implmentation of the algorithm.
     *
     *         Note: If the algorithm requires a triangulation then this
     *               method must be called after the triangulation has been
     *               preconditioned for the upcoming operations.
     */
    template <class dataType,
              class triangulationType = ttk::AbstractTriangulation>
    int computePorosity(float *outputData,
                        const dataType *inputData,
			const float *gradientData,
			const float *divergence,
                        const triangulationType *triangulation,
			const float& distance,
			const float& threshold,
			const float& margin,
			const float& maxThreshold,
			const float& gradientThreshold) const {
      // start global timer
      ttk::Timer globalTimer;
      // print horizontal separator
      this->printMsg(ttk::debug::Separator::L1); // L1 is the '=' separator	             
      // print input parameters in table format
      this->printMsg({
        {"#Threads", std::to_string(this->threadNumber_)},
        {"#Vertices", std::to_string(triangulation->getNumberOfVertices())},
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

        size_t nVertices = triangulation->getNumberOfVertices();
        
	const float sqrDistance = distance * distance;
	const float marginThreshold = (1.0f + margin) * threshold;
	const float fractionalThreshold = margin * threshold;

	std::vector<std::valarray<bool>> mask;

	//init masks
	for(int i = 0; i < this->threadNumber_; i++) {
	  mask.push_back(std::valarray<bool> (nVertices));
	}
	
	
	#ifdef TTK_ENABLE_OPENMP
        #pragma omp parallel for num_threads(this->threadNumber_)
        #endif
	for(size_t i = 0; i < nVertices; i++) {

	    bool safe = false;

	    //check if other threads haven't marked the vertex
	    for(int j = 0; j < this->threadNumber_; j++) {
	      if(!mask[j][i]) {
		safe = true;
		//mark it so that this thread claims ownership
		mask[omp_get_thread_num()][i] = true;
	      }
	    }

	    if(!safe) continue;
	    
	    outputData[i] = 0;    
	    std::vector<ttk::SimplexId> candidates;
	    std::unordered_map<ttk::SimplexId,float> distances;
	    
	    float vX,vY,vZ; //coords of center vertex
	    triangulation->getVertexPoint(i,vX,vY,vZ);
	    
	    candidates.push_back(i);

	    while(!candidates.empty()) {
	     
	      ttk::SimplexId currentVertex = candidates.back();
	      candidates.pop_back();

	      //ignore anything that might be a light artifact
	      if(inputData[currentVertex] > maxThreshold) {
		outputData[currentVertex] = -1;
		continue;
	      }
	      //check if in the list and within distance (std::find != means found)
	      if(distances.find(currentVertex) != distances.end()) continue;

	      if(gradientData[currentVertex] > gradientThreshold) continue;
	      
	      float uX,uY,uZ;
	      triangulation->getVertexPoint(currentVertex,uX,uY,uZ);
	      
	      float curr_distance = (vX - uX) * (vX - uX) + (vY - uY) * (vY - uY) + (vZ - uZ) * (vZ - uZ);
	      
	      //save a sqrt by squaring distance
	      if(curr_distance > sqrDistance) continue;
	      
	      distances.insert(std::pair<ttk::SimplexId,float>(currentVertex, curr_distance));
	      
	      size_t nNeighbors = triangulation->getVertexNeighborNumber(currentVertex);
              ttk::SimplexId neighborId;
	      
              for(size_t j = 0; j < nNeighbors; j++) {
                triangulation->getVertexNeighbor(currentVertex, j, neighborId);
		
		if(distances.find(neighborId) == distances.end()) {
		  candidates.push_back(neighborId);
		}
	      }
	  }

	    //do the threshold check here
	    //marginThreshold is (1 + margin) * threshold - basically the larger threshold
	    //fractionalThreshold is threshold * margin that is the maximum distance allowed to be counted as a fraction
	    	    
	    for(const auto& pair: distances) {

	      if(inputData[pair.first] < threshold) {
		outputData[i] += ((1 - (pair.second / sqrDistance)) * inputData[pair.first]) * divergence[pair.first];
	      } else if(inputData[pair.first] < marginThreshold) {
		outputData[i] += (marginThreshold - inputData[pair.first]) / (fractionalThreshold)  * ((1-(pair.second / sqrDistance))) * divergence[pair.first];
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
