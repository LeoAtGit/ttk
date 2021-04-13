/// TODO 4: Provide your information
///
/// \ingroup vtk
/// \class ttkKernelDensityEstimation
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the ttk::KernelDensityEstimation module.
///
/// This VTK filter uses the ttk::KernelDensityEstimation module to compute the bounding box
/// of a vtkDataSet, which is returned as a vtkUnstructuredGrid.
///
/// \param Input vtkDataSet whose bounding box will be computed.
/// \param Output vtkUnstructuredGrid that corresponds to bounding box of the
/// input.
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutputDataObject()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::KernelDensityEstimation
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkKernelDensityEstimationModule.h>

// VTK Includes
#include <ttkAlgorithm.h>

// TTK Base Includes
#include <KernelDensityEstimation.h>

class TTKKERNELDENSITYESTIMATION_EXPORT ttkKernelDensityEstimation
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
  ,
    protected ttk::KernelDensityEstimation // and we inherit from the base class
{
private:
  double ImageBounds[6]{0,1,0,1,0,1};
  double Resolution[3]{1,1,1};
  double Bandwidth{1};
  int Kernel{0};

public:

  vtkSetVector6Macro(ImageBounds, double);
  vtkGetVector6Macro(ImageBounds, double);
  vtkSetVector3Macro(Resolution, double);
  vtkGetVector3Macro(Resolution, double);
  vtkSetMacro(Bandwidth, double);
  vtkGetMacro(Bandwidth, double);
  vtkSetMacro(Kernel, int);
  vtkGetMacro(Kernel, int);

  static ttkKernelDensityEstimation *New();
  vtkTypeMacro(ttkKernelDensityEstimation, ttkAlgorithm);

protected:
  ttkKernelDensityEstimation();
  ~ttkKernelDensityEstimation() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestInformation(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector);
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
