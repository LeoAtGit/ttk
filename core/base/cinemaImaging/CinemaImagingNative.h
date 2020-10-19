/// \ingroup base
/// \class ttk::EmbreeRenderer
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.05.2020
///
/// \brief TTK %EmbreeRenderer processing package.
///
/// %EmbreeRenderer is a TTK processing package that

#pragma once

#include <Debug.h>

struct Ray {
  float* dir;
  float* origin;
  float distance;
  float u;
  float v;
};
  
namespace ttk {  
  class CinemaImagingNative : virtual public Debug {
  public:
    CinemaImagingNative(){
        this->setDebugMsgPrefix("CinemaImaging(Native)");
    }
    ~CinemaImagingNative(){};

    void multiplyByScalar(float* out, const float* a, const float& scalar) const;
    void addVectors(float* out, const float* a, const float* b) const;
    void subVectors(float* out, const float* a, const float* b) const;
    void cross(float* out, const float* a, const float* b) const;
    float dot(const float* a, const float* b) const;
    
    template<typename IT>
    bool MollerTrumbore(
      Ray &ray,
      const IT v0,
      const IT v1,
      const IT v2,
      const float* vertexCoords) const;
    
    template<typename IT>
    int renderImage(
      float* depthBuffer,
      unsigned int* primitiveIds,
      float* barycentricCoordinates,

      const size_t& nVertices,
      const float* vertexCoords,
      const size_t& nTriangles,
      const IT* connectivityList,

      const double resolution[2],
      const double camPos[3],
      const double camDirRaw[3],
      const double camUp[3],
      const double& camHeight,
      const bool& orthographicProjection
    ) const;
  };
};

template <typename IT>
int ttk::CinemaImagingNative::renderImage(
  float* depthBuffer,
  unsigned int* primitiveIds,
  float* barycentricCoordinates,

  const size_t& nVertices,
  const float* vertexCoords,
  const size_t& nTriangles,
  const IT* connectivityList,

  const double resolution[2],
  const double camPos[3],
  const double camDirRaw[3],
  const double camUp[3],
  const double& camHeight,
  const bool& orthographicProjection
) const {
  ttk::Timer timer;
  int resX = resolution[0];
  int resY = resolution[1];

  this->printMsg(
		 "Rendering Image ("+std::string(orthographicProjection ? "O" : "P")+ "|"+std::to_string(resX)+"x"+std::to_string(resY)+")",
      0,0,this->threadNumber_,
      ttk::debug::LineMode::REPLACE
  );

  // Compute camera size
  const double aspect = resolution[0] / resolution[1];
  const double camSize[2] = { aspect * camHeight, camHeight};

  const auto normalize = [](double out[3], const double in[3]){
      double temp = sqrt(in[0] * in[0] + in[1] * in[1]
                         + in[2] * in[2]);
      out[0] = in[0]/temp;
      out[1] = in[1]/temp;
      out[2] = in[2]/temp;
  };

  double camDir[3]{0,0,0};
  normalize(camDir,camDirRaw);

  // Compute camRight = camDir x CamUp
  double camRight[3]{camDir[1] * camUp[2] - camDir[2] * camUp[1],
                        camDir[2] * camUp[0] - camDir[0] * camUp[2],
                        camDir[0] * camUp[1] - camDir[1] * camUp[0]};
  normalize(camRight,camRight);

  // Compute true up std::vector
  double camUpTrue[3]{camDir[1] * (-camRight[2]) - camDir[2] * (-camRight[1]),
       camDir[2] * (-camRight[0]) - camDir[0] * (-camRight[2]),
       camDir[0] * (-camRight[1]) - camDir[1] * (-camRight[0])};
  normalize(camUpTrue,camUpTrue);

  // Compute pixel size in world coordinates
  double pixelWidthWorld = camSize[0] / resolution[0];
  double pixelHeightWorld = camSize[1] / resolution[1];

  // Optimization: precompute half of the camera size to reduce the number of
  // operations in the for loop. Include a half pixel offset (-0.5) to center
  // vertices at pixel centers
  double camWidthWorldHalf = 0.5 * camSize[0] - 0.5 * pixelWidthWorld;
  double camHeightWorldHalf = 0.5 * camSize[1] - 0.5 * pixelHeightWorld;

  // Optimization: reorient camera model to bottom left corner to reduce
  // operations in for loop
  double camPosCorner[3] = {camPos[0] - camRight[0] * camWidthWorldHalf
                              - camUpTrue[0] * camHeightWorldHalf,
                            camPos[1] - camRight[1] * camWidthWorldHalf
                              - camUpTrue[1] * camHeightWorldHalf,
                            camPos[2] - camRight[2] * camWidthWorldHalf
                              - camUpTrue[2] * camHeightWorldHalf};

  float nan = std::numeric_limits<float>::quiet_NaN();
  size_t pixelIndex = 0;
  size_t bcIndex = 0;
  
  for(int y = 0; y < resY; y++) {
      double v = ((double)y) * pixelHeightWorld;

      for(int x = 0; x < resX; x++) {
        double u = ((double)x) * pixelWidthWorld;

	depthBuffer[pixelIndex] = nan;
	primitiveIds[pixelIndex] = nan;
	barycentricCoordinates[bcIndex] = nan;
	barycentricCoordinates[bcIndex+1] = nan;
	
	
        // set origin
        float org_x = camPosCorner[0] + u * camRight[0] + v * camUpTrue[0];
        float org_y = camPosCorner[1] + u * camRight[1] + v * camUpTrue[1];
        float org_z = camPosCorner[2] + u * camRight[2] + v * camUpTrue[2];

	float ray_origin[3] = { org_x,  org_y, org_z };
	
        // set dir
        float dir_x = camDir[0];
        float dir_y = camDir[1];
        float dir_z = camDir[2];

	float ray_dir[3] = {dir_x, dir_y, dir_z};
	
	float nearestTriangle = std::numeric_limits<float>::max();
	Ray ray = {ray_dir, ray_origin, 0, 0, 0};

	
	for(int i = 0; i < nTriangles; i++)
	{
	  bcIndex = pixelIndex;
	  //go through the list of all triangles

	  //ensure all vertices are correct
	  bool wasHit = false;
	  IT v0 = connectivityList[i*3+0];
	  IT v1 = connectivityList[i*3+1];
	  IT v2 = connectivityList[i*3+2];
	  v0 *= 3; 
	  v1 *= 3; 
	  v2 *= 3; 
     	  wasHit = MollerTrumbore(ray,v0,v1,v2,vertexCoords);
    	  if(wasHit && nearestTriangle > ray.distance)
	    { 		  
		  depthBuffer[pixelIndex] = ray.distance;
		  primitiveIds[pixelIndex] = i;
		  barycentricCoordinates[bcIndex] = ray.u;
		  barycentricCoordinates[bcIndex+1] = ray.v;
		  nearestTriangle = ray.distance;
	    }
     	}
        pixelIndex++;
	bcIndex += 2;
      }
  }

  this->printMsg(
      "Rendering Image ("+std::string(orthographicProjection ? "O" : "P")+ "|"+std::to_string(resX)+"x"+std::to_string(resY)+")",
      1,timer.getElapsedTime(),this->threadNumber_
  );

  return 1;
};

template <typename IT>
bool ttk::CinemaImagingNative::MollerTrumbore(  
      Ray &ray,
      const IT v0,
      const IT v1,
      const IT v2,
      const float* vertexCoords) const
{

  constexpr float kEpsilon = 1e-8;
  
  float v0v1[3], v0v2[3], pvec[3], tvec[3], qvec[3];
  
  subVectors(v0v1,&vertexCoords[v0],&vertexCoords[v1]);
  subVectors(v0v2,&vertexCoords[v0],&vertexCoords[v2]);
  
  cross(pvec,ray.dir,v0v2);
  float det = dot(v0v1,pvec);
  if (det > -kEpsilon && det < kEpsilon) return false;
  
  float invDet = 1.0f / det;

  subVectors(tvec,&vertexCoords[v0],ray.origin);
  float u = dot(tvec, pvec) * invDet;
  if (u < 0.0 || u > 1.0) return false;
  
  cross(qvec,tvec,v0v1);
  float v = dot(ray.dir,qvec) * invDet;
  if (v < 0.0 || u + v > 1.0) return false;

  float t = dot(v0v2,qvec) * invDet;
  ray.distance = t;

  ray.u = u;
  ray.v = v;
  
  return true; 
};



