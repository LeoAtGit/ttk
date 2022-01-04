/// \ingroup vtk
/// \class ttkPathIntegrator
/// \author Emma Nilsson <emma.nilsson@liu.se>
/// \date 2021-09-20.
///
/// \brief TTK VTK-filter that wraps the ttk::PathIntegrator module.
///
///
/// \param Input vtkPolyDataset.
/// \param Output vtkPolyDataset.

#pragma once

// VTK Module
#include <ttkPathIntegratorModule.h>

// VTK Includes
#include <ttkAlgorithm.h>

// TTK Base Includes
#include <PathIntegrator.h>

class TTKPATHINTEGRATOR_EXPORT ttkPathIntegrator
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
  ,
    protected ttk::PathIntegrator // and we inherit from the base class
{
private:
  std::vector<std::vector<ttk::PathIntegrator::Point>> pointsPerTimestep;

  int nTimesteps{0};
  double TimeInterval{0.0};
  double StepLength{0.0};
  double PerlinScaleFactor{0.0};

  int VecField{0};

public:
  // Properties macros
  vtkSetMacro(nTimesteps, int);
  vtkGetMacro(nTimesteps, int);

  vtkSetMacro(TimeInterval, double);
  vtkGetMacro(TimeInterval, double);

  vtkSetMacro(StepLength, double);
  vtkGetMacro(StepLength, double);

  vtkSetMacro(PerlinScaleFactor, double);
  vtkGetMacro(PerlinScaleFactor, double);

  vtkSetMacro(VecField, int);
  vtkGetMacro(VecField, int);

  static ttkPathIntegrator *New();
  vtkTypeMacro(ttkPathIntegrator, ttkAlgorithm);

protected:
  ttkPathIntegrator();
  ~ttkPathIntegrator() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
