/// \ingroup vtk
/// \class ttkIcosphereFromPoint
/// \author Jonas Lukasczyk (jl@jluk.de)
/// \date 01.09.2019
///
/// This filter creates an icosphere with a specified radius, center, and number
/// of subdivisions at each vertex of an input vtkPointset.
///
/// \sa ttk::IcoSphere
/// \sa ttk::ttkAlgorithm

#pragma once

// VTK Module
#include <ttkIcosphereFromPointModule.h>

// VTK Includes
#include <ttkIcosphere.h>

class TTKICOSPHEREFROMPOINT_EXPORT ttkIcosphereFromPoint : public ttkIcosphere {

private:
  bool CopyPointData{true};

public:
  vtkSetMacro(CopyPointData, bool);
  vtkGetMacro(CopyPointData, bool);

  static ttkIcosphereFromPoint *New();
  vtkTypeMacro(ttkIcosphereFromPoint, ttkIcosphere);

protected:
  ttkIcosphereFromPoint();
  ~ttkIcosphereFromPoint();

  int FillInputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};