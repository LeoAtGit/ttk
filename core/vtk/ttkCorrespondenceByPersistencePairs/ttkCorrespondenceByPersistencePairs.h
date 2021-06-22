#pragma once

// VTK Module
#include <ttkCorrespondenceByPersistencePairsModule.h>

// VTK Includes
#include <ttkCorrespondenceAlgorithm.h>
#include <vtkUnstructuredGrid.h>

// TTK Base Includes
#include <CorrespondenceByPersistencePairs.h>

class TTKCORRESPONDENCEBYPERSISTENCEPAIRS_EXPORT
  ttkCorrespondenceByPersistencePairs
  : public ttkCorrespondenceAlgorithm,
    protected ttk::CorrespondenceByPersistencePairs {

public:
  static ttkCorrespondenceByPersistencePairs *New();
  vtkTypeMacro(ttkCorrespondenceByPersistencePairs, ttkCorrespondenceAlgorithm);

  vtkSetMacro(PX, double);
  vtkGetMacro(PX, double);

  vtkSetMacro(PY, double);
  vtkGetMacro(PY, double);

  vtkSetMacro(PZ, double);
  vtkGetMacro(PZ, double);

  vtkSetMacro(PE, double);
  vtkGetMacro(PE, double);

  vtkSetMacro(PS, double);
  vtkGetMacro(PS, double);

  vtkSetMacro(Alpha, double);
  vtkGetMacro(Alpha, double);

  vtkSetMacro(PVAlgorithm, int);
  vtkGetMacro(PVAlgorithm, int);

  vtkSetMacro(WassersteinMetric, const std::string &);
  vtkGetMacro(WassersteinMetric, std::string);

  vtkSetMacro(DistanceAlgorithm, const std::string &);
  vtkGetMacro(DistanceAlgorithm, std::string);

  template <typename dataType>
  int getDiagram(std::vector<std::tuple<int,
                                        ttk::CriticalType,
                                        int,
                                        ttk::CriticalType,
                                        dataType,
                                        int,
                                        dataType,
                                        float,
                                        float,
                                        float,
                                        dataType,
                                        float,
                                        float,
                                        float>> &diagram,
                 vtkUnstructuredGrid *CTPersistenceDiagram_,
                 const double spacing,
                 const int diagramNumber);

protected:
  ttkCorrespondenceByPersistencePairs();
  ~ttkCorrespondenceByPersistencePairs();

  int ComputeCorrespondences(vtkImageData *correspondenceMatrix,
                             vtkDataObject *inputDataObjects0,
                             vtkDataObject *inputDataObjects1) override;

private:
  // Metric weights
  double PX{1};
  double PY{1};
  double PZ{1};
  double PE{1}; // extrema
  double PS{1}; // saddles
  // Metric config
  double Alpha{1.0}; // power
  std::string DistanceAlgorithm{"ttk"}; // distance between PPs
  std::string WassersteinMetric{"2"}; // wass vs inf (bottleneck)
  int PVAlgorithm{-1};
};
