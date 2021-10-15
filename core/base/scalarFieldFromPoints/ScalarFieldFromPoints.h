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

    typedef double(*KERNEL)(const double&);

    static double Linear(const double& u) {
      return u>=1 ? 0 : 1-u;
    };

    static double Epanechnikov(const double& u) {
      return u>=1 ? 0 : 0.75 - 0.75*u*u;
    };

    static double Gaussian(const double& u) {
      constexpr double c = (1.0/std::sqrt(pow(2.0*3.14159265359, 3)));
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
    int computeScalarField(
      double* outputData,
      const double* pointCoordiantes,
      const size_t nPoints,
      const float bandwidth,
      const TT* triangulation
    ) const {
      ttk::Timer timer;

      this->printMsg("Computing scalar field",
                       0, // progress form 0-1
                       0, // elapsed time so far
                       this->threadNumber_, ttk::debug::LineMode::REPLACE);

      // compute the average of each vertex in parallel
      size_t nVertices = triangulation->getNumberOfVertices();
      #ifdef TTK_ENABLE_OPENMP
      #pragma omp parallel for num_threads(this->threadNumber_)
      #endif
      for(size_t i = 0; i < nVertices; i++) {
        float x,y,z;
        triangulation->getVertexPoint(i, x,y,z);

        double& f = outputData[i];
        f = 0;

        for(size_t j=0; j<nPoints; j++){
          const size_t& j3 = j*3;

          double dx = x - pointCoordiantes[j3+0];
          double dy = y - pointCoordiantes[j3+1];
          double dz = z - pointCoordiantes[j3+2];

          // Calculate norm of vector to be used in kernel
          // Divide with bandwidth bc. KDE
          double u = std::sqrt(dx*dx + dy*dy + dz*dz)/bandwidth;

          f += k(u);
        }
      }

      // print the progress of the current subprocedure with elapsed time
      this->printMsg("Computing scalar field",
                     1, // progress
                     timer.getElapsedTime(), this->threadNumber_);

      return 1; // return success
    }

    int computeTotalKDE(
      double* outputData,
      const double* kdeByType,

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
      double* outputData,

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

      const double width = bounds[1]-bounds[0];
      const double height = bounds[3]-bounds[2];

      const double resXm1 = res[0]-1;
      const double resYm1 = res[1]-1;
      const int iResX = res[0];
      const int iResY = res[1];
      const int iResXm1 = iResX-1;
      const int iResYm1 = iResY-1;

      // Calculate width and height for the output kernel, and the number of points
      // in each direction that should be given scalar values
      const double dx = width / resXm1;
      const double dy = height / resYm1;

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

          double count_it = count[i*nTypes+t];
          if(count_it<1.0)
            continue;

          const int xi = i%iResX;
          const int yi = i/iResX;

          const int x0 = std::max(0, std::min(iResXm1, xi - kdx));
          const int x1 = std::max(0, std::min(iResXm1, xi + kdx));
          const int y0 = std::max(0, std::min(iResYm1, yi - kdy));
          const int y1 = std::max(0, std::min(iResYm1, yi + kdy));

          const int xSize = iResX*nTypes;
          
          // for all points in the bandwidth interval, calculate scalar value
          for(int x = x0; x<=x1; x++){
            for(int y = y0; y<=y1; y++){

              double xxx = (x-xi)*dx;
              double yyy = (y-yi)*dy;
              const double u = std::sqrt(xxx*xxx + yyy*yyy)/bandwidth;
              const double ku = k(u);

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

      double fac = nTotalEvents*bandwidth;
      for(int i=0, j=nPixels*nTypes; i<j; i++){
        outputData[i] /= fac;
      }

      // print the progress of the current subprocedure with elapsed time
      this->printMsg("Computing KDEs by Category",
                     1, // progress
                     timer.getElapsedTime(), this->threadNumber_);

      return 1; // return success
    }

    template<KERNEL k>
    int computeScalarField3D(
      double* outputData,
      const int* count,
      const double bandwidth,
      const double* bounds,
      const double* res
    ) const {
      ttk::Timer timer;

      this->printMsg("Computing Scalar Field 3D",
                       0, // progress form 0-1
                       0, // elapsed time so far
                       this->threadNumber_, ttk::debug::LineMode::REPLACE);

      const int nPixels = res[0]*res[1]*res[2];

      const double width = bounds[1]-bounds[0];
      const double height = bounds[3]-bounds[2];
      const double depth = bounds[5]-bounds[4];
      const double resXm1 = res[0]-1;
      const double resYm1 = res[1]-1;
      const double resZm1 = res[2]-1;
      const int iResX = res[0];
      const int iResY = res[1];
      const int iResZ = res[2];
      const int iResXm1 = iResX-1;
      const int iResYm1 = iResY-1;
      const int iResZm1 = iResZ-1;

      // Calculate width, height and depth for the output kernel, and the number of points
      // in each direction that should be given scalar values
      const double dx = width / resXm1;
      const double dy = height / resYm1;
      const double dz = depth / resZm1;

      const int kdx = 4 * floor(bandwidth/dx+0.5);
      const int kdy = 4 * floor(bandwidth/dy+0.5);
      const int kdz = 4 * floor(bandwidth/dz+0.5);

      // bandwidth is the same for each dimension --> same variance
      const double variance = bandwidth;
      const double stdev = sqrt(variance);

      // const int kdx = 8 * stdev;
      // const int kdy = 8 * stdev;
      // const int kdz = 8 * stdev;

      this->printMsg(std::to_string(kdx));
      this->printMsg(std::to_string(dx));
      this->printMsg(std::to_string(variance));
      this->printMsg("-------");

      // clear
      for(int i=0, j=nPixels; i<j; i++){
        outputData[i] = 0;
      }

      // compute scalar field
      #ifdef TTK_ENABLE_OPENMP
      #pragma omp parallel for num_threads(this->threadNumber_)
      #endif
      for(int i = 0; i < nPixels; i++) {

        double count_it = count[i];
        if(count_it<1.0)
          continue;

        const int xi = i % iResX;
        const int yi = (i / iResX) % iResY;
        const int zi = (i / (iResX * iResY));

        const int x0 = std::max(0, std::min(iResXm1, xi - kdx));
        const int x1 = std::max(0, std::min(iResXm1, xi + kdx));
        const int y0 = std::max(0, std::min(iResYm1, yi - kdy));
        const int y1 = std::max(0, std::min(iResYm1, yi + kdy));
        const int z0 = std::max(0, std::min(iResZm1, zi - kdz));
        const int z1 = std::max(0, std::min(iResZm1, zi + kdz));

        if (i == 1977138) {
          this->printMsg(std::to_string(x0));
          this->printMsg(std::to_string(x1));
          this->printMsg(std::to_string(xi));
          this->printMsg("-------");
        }
        
        // for all points in the bandwidth interval, calculate scalar value
        for(int x = x0; x <= x1; x++){
          for(int y = y0; y <= y1; y++){
            for (int z = z0; z <= z1; z++) {
              double xxx = (x - xi) * dx;
              double yyy = (y - yi) * dy;
              double zzz = (z - zi) * dz;
              const double u = std::sqrt(xxx * xxx + yyy * yyy + zzz * zzz)/variance;
              const double ku = k(u);// / pow(stdev, 3);
              if (x == x0 && y == y0 && z == z0 && i == 1977138) {
                this->printMsg(std::to_string(xxx));
                this->printMsg(std::to_string(u));
                this->printMsg(std::to_string(ku));
                this->printMsg("-------");
              }

              outputData[z * iResY * iResX + y * iResX + x ] += ku;
            }
          }
        }
      }

      // print the progress of the current subprocedure with elapsed time
      this->printMsg("Computing Scalar Field 3D",
                     1, // progress
                     timer.getElapsedTime(), this->threadNumber_);

      return 1; // return success
    }

    int computeCounts(
      int* counts,
      const double* pointCoordiantes,
      const double* bounds,
      const double* res,
      const size_t& nPoints,
      const size_t& nPixels
    ) const {
      ttk::Timer timer;

      this->printMsg("Computing Counts", 0, 0, this->threadNumber_, ttk::debug::LineMode::REPLACE);

      // clear
      size_t nTotal = nPixels;
      for(size_t i=0; i<nTotal; i++){
        counts[i] = 0;
      }

      const float width = bounds[1]-bounds[0];
      const float height = bounds[3]-bounds[2];
      const float depth = (bounds[5] > 1) ? bounds[5]-bounds[4] : 0;
      const float resX = res[0];
      const float resXm1 = res[0]-1;
      const float resYm1 = res[1]-1;
      const float resZm1 = (res[2] > 1) ? res[2]-1 : 1;

      const int iResX = resX;
      const int iResY = res[1];
      const float dx = width / resXm1;
      const float dy = height / resYm1;
      const float dz = depth / resZm1;
      const float dx2 = dx/2.0;
      const float dy2 = dy/2.0;
      const float dz2 = dz/2.0;

      const float x0 = bounds[0] - dx2;
      const float y0 = bounds[2] - dy2;
      const float z0 = bounds[4] - dz2;


      #ifdef TTK_ENABLE_OPENMP
      #pragma omp parallel for num_threads(this->threadNumber_)
      #endif
      for(size_t i=0; i<nPoints; i++){

        const float& x = pointCoordiantes[i*3+0];
        const float& y = pointCoordiantes[i*3+1];
        const float& z = pointCoordiantes[i*3+2];

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
        const int zi = std::min(
          resZm1,
          std::max(
            0.0f,
            floor((z - z0)/dz)
          )
        );

        const int pixelIndex = zi * iResX * iResY + yi * iResX + xi;
        const int countIndex = pixelIndex;

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
