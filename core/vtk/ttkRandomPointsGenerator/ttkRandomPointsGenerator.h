/// \ingroup vtk
/// \class ttkRandomPointsGenerator
/// \author Emma Nilsson <emma.nilsson@liu.se>
/// \date 2021-08-31.
///
/// \sa ttk::RandomPointsGenerator
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkRandomPointsGeneratorModule.h>

// VTK Includes
#include <ttkAlgorithm.h>

// TTK Base Includes
#include <RandomPointsGenerator.h>

class TTKRANDOMPOINTSGENERATOR_EXPORT ttkRandomPointsGenerator
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
  ,
    protected ttk::RandomPointsGenerator // and we inherit from the base class
{
private:
  int nPoints{0};
  int PositionDomain[3]{0, 0, 0};
  double Amplitude[2]{0.0, 1.0};
  double Spread[2]{0.0, 1.0};
  int RandomSeed{0};

public:
  // Properties macros
  vtkSetMacro(nPoints, int);
  vtkGetMacro(nPoints, int);

  vtkSetVector3Macro(PositionDomain, int);
  vtkGetVector3Macro(PositionDomain, int);

  vtkSetVector2Macro(Amplitude, double);
  vtkGetVector2Macro(Amplitude, double);

  vtkSetVector2Macro(Spread, double);
  vtkGetVector2Macro(Spread, double);

  vtkSetMacro(RandomSeed, int);
  vtkGetMacro(RandomSeed, int);

  static ttkRandomPointsGenerator *New();
  vtkTypeMacro(ttkRandomPointsGenerator, ttkAlgorithm);

protected:

  ttkRandomPointsGenerator();
  ~ttkRandomPointsGenerator() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
