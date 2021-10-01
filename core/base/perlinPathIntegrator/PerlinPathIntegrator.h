/// \ingroup base
/// \class ttk::PerlinPathIntegrator
/// \author Emma Nilsson <emma.nilsson@liu.se>
/// \date 2021-09-20.

#pragma once

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>
#include <PerlinNoise.h>

namespace ttk {

  /**
   * The PerlinPathIntegrator class provides methods to compute for each vertex of a
   * triangulation the average scalar value of itself and its direct neighbors.
   */
  class PerlinPathIntegrator : virtual public Debug {

  public:
    PerlinPathIntegrator();

    struct Point {
      int pointId{-1};
      int timestep{-1};
      double x{0.0};
      double y{0.0};
      double z{0.0};
      double v[3]{0.0, 0.0, 0.0};

      Point() {

      }

      Point(const Point& p) {
        pointId = p.pointId;
        timestep = p.timestep;
        x = p.x;
        y = p.y;
        z = p.z;
        v[0] = p.v[0];
        v[1] = p.v[1];
        v[2] = p.v[2];
      }

      Point operator+(const Point& a) const {
          Point p;
          p.pointId = pointId;
          p.timestep = timestep;
          p.x = x + a.x;
          p.y = y + a.y;
          p.z = z + a.z;

          return p;
      }

      Point operator*(double k) {
          Point p;
          p.pointId = pointId;
          p.timestep = timestep;
          p.x = k * x;
          p.y = k * y;
          p.z = k * z;
          return p;

      }

      void setVelocity(const double vel[3]) {
        v[0] = vel[0];
        v[1] = vel[1];
        v[2] = vel[2];
      }

    };

    // Set functions for class variables
    void setDomainDimension(const int dims[3]) {
      dimensions_[0] = dims[0];
      dimensions_[1] = dims[1];
      dimensions_[2] = dims[2];
    }

    void setStepLength(const double stepLength) {
      h_ = stepLength;
    }

    void setPerlinScaleFactor(const double psf) {
      psf_[0] = {(dimensions_[0] - 1) / psf};
      psf_[1] = {(dimensions_[1] - 1) / psf};
      psf_[2] = {(dimensions_[2] - 1) / psf};
    }    

    Point sampleVectorField(const Point& p, double t) {
      Point dv(p);
      pn.perlin4D<double>(p.x/psf_[0], p.y/psf_[1], p.z/psf_[2], t, dv.x);
      pn.perlin4D<double>(p.z/psf_[2], p.x/psf_[0], p.y/psf_[1], t, dv.y);
      pn.perlin4D<double>(p.y/psf_[1], p.z/psf_[2], p.x/psf_[0], t, dv.z);

      return dv;
    }

    int RK4(Point& prevP, Point& newP, double time) {
      // Call perlin function to execute vector field
      // add const point&
      Point q1 = sampleVectorField(prevP, time) * h_;
      Point q2 = sampleVectorField(prevP + (q1 * 0.5), time) * h_;
      Point q3 = sampleVectorField(prevP + (q2 * 0.5), time) * h_;
      Point q4 = sampleVectorField(prevP + q3, time) * h_; 

      Point vel = (q1 + q2 * 2 + q3 * 2 + q4) * (1.0/6);
      newP = prevP + vel;

      // Set velocity of previous point
      double v[3] = {vel.x, vel.y, vel.z};
      prevP.setVelocity(v);

      return 1;
    }

    template <class dataType>
    int integrate(
      const std::vector<Point>& initPoints,
      std::vector<std::vector<Point>>& outPoints,
      const int nTimesteps,
      const double timeInterval,
      const int dims[3],
      const double stepLength,
      const double psf
    ) {
      // Set class variables
      setDomainDimension(dims);
      setStepLength(stepLength);
      setPerlinScaleFactor(psf);


      // Integrate the paths of the initial points by moving the points
      // along the vector field for all timesteps
      for (int i = 0; i < nTimesteps - 1; i++) {
        std::vector<Point>& curPoints = outPoints[i];
        double time = i * timeInterval;
        for (size_t j = 0; j < outPoints[i].size(); j++) {
          Point newP;
          auto& curPoint = curPoints[j];

          // Integrate using RK4
          RK4(curPoint, newP, time);

          // Check if points is within the domain dimensions
          if ((newP.x > 0.0 && newP.x < dimensions_[0]) && 
            (newP.y > 0.0 && newP.y < dimensions_[1]) &&
            (newP.z > 0.0 && newP.z < dimensions_[2])
          ) {
            // Add point to the next time-step
            newP.timestep = i+1;
            newP.pointId = curPoint.pointId;
            outPoints[i+1].push_back(newP);

          }
        }
      }
      return 1;
    }

  protected:
  int dimensions_[3]{};
  double h_{};
  double psf_[3]{};
  PerlinNoise pn;

  }; // PerlinPathIntegrator class

} // namespace ttk
