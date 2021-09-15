/// \ingroup vtk
/// \class ttkPerlinVectorField
/// \author Emma Nilsson <emma.nilsson@liu.se>
/// \date 2021-08-31.
///
/// \sa ttk::PerlinVectorField
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkPerlinVectorFieldModule.h>

// VTK Includes
#include <ttkAlgorithm.h>

// TTK Base Includes
#include <PerlinVectorField.h>
#include<PerlinNoise.h>

class TTKPERLINVECTORFIELD_EXPORT ttkPerlinVectorField
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
  ,
    protected ttk::PerlinVectorField // and we inherit from the base class
{
private:
  int ImageDimension[3]{0, 0, 0};
  int Scale{0};

public:
  // Properties macros
  vtkSetVector3Macro(ImageDimension, int);
  vtkGetVector3Macro(ImageDimension, int);

  vtkSetMacro(Scale, int);
  vtkGetMacro(Scale, int);

  static ttkPerlinVectorField *New();
  vtkTypeMacro(ttkPerlinVectorField, ttkAlgorithm);

protected:

  ttkPerlinVectorField();
  ~ttkPerlinVectorField() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  // Needs to be filled in for vtkImageData
  int RequestInformation(vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector *outputVector) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
