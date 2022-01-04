/// \ingroup base
/// \class ttk::PathIntegrator
/// \author Emma Nilsson <emma.nilsson@liu.se>
/// \date 2021-09-20.

#pragma once

// ttk common includes
#include <Debug.h>
#include <PerlinNoise.h>
#include <Triangulation.h>

namespace ttk {

  /**
   * The PathIntegrator class provides methods to compute for each vertex of a
   * triangulation the average scalar value of itself and its direct neighbors.
   */
  class PathIntegrator : virtual public Debug {

  public:
    PathIntegrator();

    struct Point {
      int pointId{-1};
      int timestep{-1};
      double x{0.0};
      double y{0.0};
      double z{0.0};
      double v[3]{0.0, 0.0, 0.0};
      double amplitude{0.0};
      double spread{0.0};
      double rate{0.0};
      bool outsideDomain{false};

      Point() {
      }

      Point(const double &px, const double &py, const double &pz) {
        x = px;
        y = py;
        z = pz;
      }

      Point(const Point &p) {
        pointId = p.pointId;
        timestep = p.timestep;
        x = p.x;
        y = p.y;
        z = p.z;
        v[0] = p.v[0];
        v[1] = p.v[1];
        v[2] = p.v[2];
        amplitude = p.amplitude;
        spread = p.spread;
        rate = p.rate;
        outsideDomain = p.outsideDomain;
      }

      Point &operator=(const Point &p) {
        pointId = p.pointId;
        timestep = p.timestep;
        x = p.x;
        y = p.y;
        z = p.z;
        v[0] = p.v[0];
        v[1] = p.v[1];
        v[2] = p.v[2];
        amplitude = p.amplitude;
        spread = p.spread;
        rate = p.rate;
        outsideDomain = p.outsideDomain;
        return *this;
      }

      Point operator+(const Point &a) const {
        Point p;
        p.pointId = pointId;
        p.timestep = timestep;
        p.amplitude = amplitude;
        p.spread = spread;
        p.rate = rate;
        p.outsideDomain = outsideDomain;
        p.x = x + a.x;
        p.y = y + a.y;
        p.z = z + a.z;

        return p;
      }

      Point operator*(double k) {
        Point p;
        p.pointId = pointId;
        p.timestep = timestep;
        p.amplitude = amplitude;
        p.spread = spread;
        p.rate = rate;
        p.outsideDomain = outsideDomain;
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

    enum class VectorField {
      PerlinPerturbed,
      PerlinGradient,
      PosDiagonal,
      PosX
    };

    // Set functions for class variables
    void setDomainDimension(const int dims[3]) {
      dimensions_[0] = dims[0];
      dimensions_[1] = dims[1];
      dimensions_[2] = dims[2];
    }

    void setVirtualBoxDimension(const double &maxSpread) {
      int boxLimit = ceil(3 * std::sqrt(maxSpread));

      // X min and max
      vbDimensions_[0] = 0 - boxLimit;
      vbDimensions_[1] = dimensions_[0] + boxLimit;
      // Y min and max
      vbDimensions_[2] = 0 - boxLimit;
      vbDimensions_[3] = dimensions_[1] + boxLimit;
      // Z min and max
      vbDimensions_[4] = 0 - boxLimit;
      vbDimensions_[5] = dimensions_[2] + boxLimit;
    }

    void setStepLength(const double stepLength) {
      h_ = stepLength;
    }

    void setPerlinScaleFactor(const double psf) {
      psf_[0] = {(dimensions_[0] - 1) / psf};
      psf_[1] = {(dimensions_[1] - 1) / psf};
      psf_[2] = {(dimensions_[2] - 1) / psf};
    }

    void setVectorField(const PathIntegrator::VectorField &vf) {
      vf_ = vf;
    }

    Point sampleVectorField(const Point &p, double t) {
      Point dv(p);
      switch(vf_) {
        case PathIntegrator::VectorField::PerlinPerturbed: {
          // Use perturbation in the 3D domain for the vector field
          pn.perlin4D<double>(
            p.x / psf_[0], p.y / psf_[1], p.z / psf_[2], t, dv.x);
          pn.perlin4D<double>(
            p.z / psf_[2], p.x / psf_[0], p.y / psf_[1], t, dv.y);
          pn.perlin4D<double>(
            p.y / psf_[1], p.z / psf_[2], p.x / psf_[0], t, dv.z);
          break;
        }
        case PathIntegrator::VectorField::PerlinGradient: {
          // Calculate gradient using finite differences
          double s1, s2 = 0.0;
          pn.perlin4D<double>(
            (p.x - h_) / psf_[0], p.y / psf_[1], p.z / psf_[2], t, s1);
          pn.perlin4D<double>(
            (p.x + h_) / psf_[0], p.y / psf_[1], p.z / psf_[2], t, s2);
          dv.x = (s2 - s1) / (2 * h_);

          pn.perlin4D<double>(
            p.x / psf_[0], (p.y - h_) / psf_[1], p.z / psf_[2], t, s1);
          pn.perlin4D<double>(
            p.x / psf_[0], (p.y + h_) / psf_[1], p.z / psf_[2], t, s2);
          dv.y = (s2 - s1) / (2 * h_);

          pn.perlin4D<double>(
            p.x / psf_[0], p.y / psf_[1], (p.z - h_) / psf_[2], t, s1);
          pn.perlin4D<double>(
            p.x / psf_[0], p.y / psf_[1], (p.z + h_) / psf_[2], t, s2);
          dv.z = (s2 - s1) / (2 * h_);
          break;
        }
        case PathIntegrator::VectorField::PosDiagonal: {
          // Go (1, 1, 1) along positive diagonal
          dv.x = 1.0;
          dv.y = 1.0;
          dv.z = 1.0;
          break;
        }
        case PathIntegrator::VectorField::PosX: {
          // Go (1, 0, 0)
          dv.x = 1.0;
          dv.y = 0.0;
          dv.z = 0.0;
        }
      }

      return dv;
    }

    int RK4(Point &prevP, Point &newP, double time) {
      // Call perlin function to execute vector field
      // add const point&
      if((prevP.x > 0.0 && prevP.x < dimensions_[0])
         && (prevP.y > 0.0 && prevP.y < dimensions_[1])
         && (prevP.z > 0.0 && prevP.z < dimensions_[2])) {
        // Check if point has been outside domain and set that it has re-entered
        if(prevP.outsideDomain) {
          prevP.outsideDomain = false;
        }

        Point q1 = sampleVectorField(prevP, time) * h_;
        Point q2 = sampleVectorField(prevP + (q1 * 0.5), time) * h_;
        Point q3 = sampleVectorField(prevP + (q2 * 0.5), time) * h_;
        Point q4 = sampleVectorField(prevP + q3, time) * h_;

        Point vel = (q1 + q2 * 2 + q3 * 2 + q4) * (1.0 / 6);
        newP = prevP + vel;

        // Set velocity of previous point
        double v[3] = {vel.x, vel.y, vel.z};
        prevP.setVelocity(v);
      } else {
        // Check if domain has been left for this point
        if(!prevP.outsideDomain) {
          // Set new velocity direction depending on what face the point has
          // exited
          prevP.v[0] = (prevP.x < 0.0)              ? -1
                       : (prevP.x > dimensions_[0]) ? 1
                                                    : 0;
          prevP.v[1] = (prevP.y < 0.0)              ? -1
                       : (prevP.y > dimensions_[1]) ? 1
                                                    : 0;
          prevP.v[2] = (prevP.z < 0.0)              ? -1
                       : (prevP.z > dimensions_[2]) ? 1
                                                    : 0;

          // Set point to be outside domain
          prevP.outsideDomain = true;
        }

        newP = prevP + Point(prevP.v[0], prevP.v[1], prevP.v[2]) * (1.0 / 6);

        // Transform coordinates into virtual box coordinates and make periodic,
        // then transform back to domain coordinates
        double Xvb = (newP.x - vbDimensions_[0])
                     - floor((newP.x - vbDimensions_[0])
                             / double(vbDimensions_[1] - vbDimensions_[0]))
                         * double(vbDimensions_[1] - vbDimensions_[0]);
        newP.x = Xvb - abs(vbDimensions_[0]);

        double Yvb = (newP.y - vbDimensions_[2])
                     - floor((newP.y - vbDimensions_[2])
                             / double(vbDimensions_[3] - vbDimensions_[2]))
                         * double(vbDimensions_[3] - vbDimensions_[2]);
        newP.y = Yvb - abs(vbDimensions_[2]);

        double Zvb = (newP.z - vbDimensions_[4])
                     - floor((newP.z - vbDimensions_[4])
                             / double(vbDimensions_[5] - vbDimensions_[4]))
                         * double(vbDimensions_[5] - vbDimensions_[4]);
        newP.z = Zvb - abs(vbDimensions_[4]);

        // Set velocity of new point
        double v[3] = {prevP.v[0], prevP.v[1], prevP.v[2]};
        newP.setVelocity(v);
      }

      return 1;
    }

    template <class dataType>
    int integrate(std::vector<std::vector<Point>> &outPoints,
                  const int nTimesteps,
                  const double timeInterval,
                  const int dims[3],
                  const double stepLength,
                  const double psf,
                  const double maxSpread,
                  const VectorField &vf) {
      // Set class variables
      setDomainDimension(dims);
      setVirtualBoxDimension(maxSpread);
      setStepLength(stepLength);
      setPerlinScaleFactor(psf);
      setVectorField(vf);

      int maxPointId = outPoints[0].size() - 1;

      // Integrate the paths of the initial points by moving the points
      // along the vector field for all timesteps
      ttk::Timer timer;
      this->printMsg("Integrating " + std::to_string(maxPointId + 1)
                       + " particles along vector field",
                     0, 0, this->threadNumber_, debug::LineMode::REPLACE);

      for(int i = 0; i < nTimesteps - 1; i++) {
        std::vector<Point> &curPoints = outPoints[i];
        double time = i * timeInterval;
        for(size_t j = 0; j < outPoints[i].size(); j++) {
          Point newP;
          auto &curPoint = curPoints[j];

          // Integrate using RK4
          RK4(curPoint, newP, time);

          // Add point to the next time-step
          newP.timestep = i + 1;
          newP.pointId = curPoint.pointId;
          outPoints[i + 1].push_back(newP);
        }
      }

      this->printMsg("Integrating " + std::to_string(maxPointId + 1)
                       + " particles along vector field",
                     1, timer.getElapsedTime(), this->threadNumber_);

      return 1;
    }

  protected:
    int dimensions_[3]{};
    int vbDimensions_[6]{};
    double h_{};
    double psf_[3]{};
    VectorField vf_{};
    PerlinNoise pn;

  }; // PathIntegrator class

} // namespace ttk
