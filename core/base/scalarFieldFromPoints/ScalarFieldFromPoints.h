/// TODO 1: Provide your information
///
/// \ingroup base
/// \class ttk::ScalarFieldFromPoints
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// This module defines the %ScalarFieldFromPoints class that computes for each vertex of a
/// triangulation the average scalar value of itself and its direct neighbors.
///
/// \b Related \b publication: \n
/// 'ScalarFieldFromPoints'
/// Jonas Lukasczyk and Julien Tierny.
/// TTK Publications.
/// 2020.
///

#pragma once

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>

#include <math.h>

namespace ttk {

  /**
   * The ScalarFieldFromPoints class provides methods to compute for each vertex of a
   * triangulation the average scalar value of itself and its direct neighbors.
   */

  class ScalarFieldFromPoints : virtual public Debug {

  public:

    typedef float(*KERNEL)(const float&);

    static float Linear(const float& u) {
      return u>=1 ? 0 : 1-u;
    };

    static float Epanechnikov(const float& u) {
      return u>=1 ? 0 : 0.75 - 0.75*u*u;
    };

    static float Gaussian(const float& u) {
      constexpr float c = (1.0/std::sqrt(2.0*3.14159265359));
      return c * exp(-0.5*u*u);
    };


    // static float Gaussian(const float& x, const float& y) {
    //   constexpr float c = (1.0/std::sqrt(2.0*3.14159265359));
    //   return c * exp(-0.5*u);
    // };

    ScalarFieldFromPoints() {
      this->setDebugMsgPrefix(
        "ScalarFieldFromPoints");
    };
    ~ScalarFieldFromPoints(){};


    // float Epanechnikov(const float& u) const {
    //   return u>=1 ? 0 : 0.75 - 0.75*u*u;
    // };

    template<typename TT, KERNEL k>
    int computeKDE(
      float* outputData,
      const float* pointCoordiantes,
      const size_t nPoints,
      const float bandwidth,
      const TT* triangulation
    ) const {
      ttk::Timer timer;

      this->printMsg("Computing KDE",
                       0, // progress form 0-1
                       0, // elapsed time so far
                       this->threadNumber_, ttk::debug::LineMode::REPLACE);

      float c = 1.0 / (((float)nPoints)*bandwidth*bandwidth);

      // compute the average of each vertex in parallel
      size_t nVertices = triangulation->getNumberOfVertices();
      #ifdef TTK_ENABLE_OPENMP
      #pragma omp parallel for num_threads(this->threadNumber_)
      #endif
      for(size_t i = 0; i < nVertices; i++) {
        float x,y,z;
        triangulation->getVertexPoint(i, x,y,z);

        float& f = outputData[i];
        f = 0;

        for(size_t j=0; j<nPoints; j++){
          const size_t& j3 = j*3;

          float dx = x - pointCoordiantes[j3+0];
          float dy = y - pointCoordiantes[j3+1];
          float dz = z - pointCoordiantes[j3+2];

          float u = std::sqrt(dx*dx + dy*dy + dz*dz)/bandwidth;

          f += k(u);
        }
        f *= c;
      }

      // print the progress of the current subprocedure with elapsed time
      this->printMsg("Computing KDE",
                     1, // progress
                     timer.getElapsedTime(), this->threadNumber_);

      return 1; // return success
    }

    int computeTotalKDE(
      float* outputData,
      const float* kdeByType,

      const size_t& nPixels,
      const size_t& nTypes
    ) const {
      ttk::Timer timer;

      this->printMsg("Computing Total KDE",
                       0, // progress form 0-1
                       0, // elapsed time so far
                       this->threadNumber_, ttk::debug::LineMode::REPLACE);

      #ifdef TTK_ENABLE_OPENMP
      #pragma omp parallel for num_threads(this->threadNumber_)
      #endif
      for(size_t i=0; i<nPixels; i++){
        float sum = 0;
        for(size_t t=0; t<nTypes; t++)
          sum+=kdeByType[i*nTypes+t];
        outputData[i] = sum;
      }

      this->printMsg("Computing KDEs by Category",
             1, // progress
             timer.getElapsedTime(), this->threadNumber_);

      return 1;
    }


    template<KERNEL k>
    int computeKDE2(
      float* outputData,

      const int* count,
      const size_t& nTypes,
      const float bandwidth,
      const double* bounds,
      const double* res
    ) const {
      ttk::Timer timer;

      this->printMsg("Computing KDEs by Category",
                       0, // progress form 0-1
                       0, // elapsed time so far
                       this->threadNumber_, ttk::debug::LineMode::REPLACE);

      const int nPixels = res[0]*res[1];

      const float width = bounds[1]-bounds[0];
      const float height = bounds[3]-bounds[2];

      const float resXm1 = res[0]-1;
      const float resYm1 = res[1]-1;
      const int iResX = res[0];
      const int iResY = res[1];
      const int iResXm1 = iResX-1;
      const int iResYm1 = iResY-1;

      const float dx = width / resXm1;
      const float dy = height / resYm1;

      const int kdx = floor(bandwidth/dx+0.5);
      const int kdy = floor(bandwidth/dy+0.5);

      // clear
      for(int i=0, j=nPixels*nTypes; i<j; i++){
        outputData[i] = 0;
      }

      // compute kde
      #ifdef TTK_ENABLE_OPENMP
      #pragma omp parallel for num_threads(this->threadNumber_)
      #endif
      for(size_t t=0; t<nTypes; t++){

        for(int i = 0; i < nPixels; i++) {

          float count_it = count[i*nTypes+t];
          if(count_it<1.0)
            continue;

          const int xi = i%iResX;
          const int yi = i/iResX;

          const int x0 = std::max(0, std::min(iResXm1, xi - kdx));
          const int x1 = std::max(0, std::min(iResXm1, xi + kdx));
          const int y0 = std::max(0, std::min(iResYm1, yi - kdy));
          const int y1 = std::max(0, std::min(iResYm1, yi + kdy));

          const int xSize = iResX*nTypes;

          for(int x = x0; x<=x1; x++){
            for(int y = y0; y<=y1; y++){

              float xxx = (x-xi)*dx;
              float yyy = (y-yi)*dy;
              const float u = std::sqrt(xxx*xxx + yyy*yyy)/bandwidth;
              const float ku = k(u);

              outputData[ y*xSize + x*nTypes +t ] += count_it*ku;
            }
          }
        }
      }

      // normalize
      // sum up total number of events
      int nTotalEvents = 0;
      for(int i=0, j=nPixels*nTypes; i<j; i++){
        nTotalEvents += count[i];
      }

      float fac = nTotalEvents*bandwidth;
      for(int i=0, j=nPixels*nTypes; i<j; i++){
        outputData[i] /= fac;
      }

      // print the progress of the current subprocedure with elapsed time
      this->printMsg("Computing KDEs by Category",
                     1, // progress
                     timer.getElapsedTime(), this->threadNumber_);

      return 1; // return success
    }

    int computeCounts(
      int* counts,

      const float* pointCoordiantes,
      const unsigned char* types,
      const size_t& nTypes,
      const double* bounds,
      const double* res,
      const size_t& nPoints,
      const size_t& nPixels
    ) const {
      ttk::Timer timer;

      this->printMsg("Computing Counts", 0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);

      // clear
      size_t nTotal = nPixels*nTypes;
      for(size_t i=0; i<nTotal; i++){
        counts[i] = 0;
      }

      const float width = bounds[1]-bounds[0];
      const float height = bounds[3]-bounds[2];
      const float resX = res[0];
      const float resXm1 = res[0]-1;
      const float resYm1 = res[1]-1;

      const int iResX = resX;
      const float dx = width / resXm1;
      const float dy = height / resYm1;
      const float dx2 = dx/2.0;
      const float dy2 = dy/2.0;

      const float x0 = bounds[0] - dx2;
      const float y0 = bounds[2] - dy2;

      #ifdef TTK_ENABLE_OPENMP
      #pragma omp parallel for num_threads(this->threadNumber_)
      #endif
      for(size_t i=0; i<nPoints; i++){

        const float& x = pointCoordiantes[i*3+0];
        const float& y = pointCoordiantes[i*3+1];

        const int xi = std::min(
          resXm1,
          std::max(
            0.0f,
            floor( (x - x0)/dx )
          )
        );
        const int yi = std::min(
          resYm1,
          std::max(
            0.0f,
            floor( (y - y0)/dy )
          )
        );

        const int pixelIndex = yi * iResX + xi;
        const int countIndex = pixelIndex*nTypes + ((int)types[i]);

        #ifdef TTK_ENABLE_OPENMP
        #pragma omp atomic update
        #endif
        counts[countIndex]++;
      }

      // print the progress of the current subprocedure with elapsed time
      this->printMsg("Computing Counts",
                     1, // progress
                     timer.getElapsedTime(), this->threadNumber_);

      return 1; // return success
    }

  }; // ScalarFieldFromPoints class

} // namespace ttk
