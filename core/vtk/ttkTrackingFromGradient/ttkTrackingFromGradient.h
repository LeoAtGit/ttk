/// TODO 4: Provide your information

#pragma once

//#include <Debug.h>

// VTK Module
#include <ttkTrackingFromGradientModule.h>

// VTK Includes
#include <ttkAlgorithm.h>

// TTK Base Includes
#include <TrackingFromGradient.h>

#include <vtkSmartPointer.h>

class vtkDataSet;
class vtkImageData;

class TTKTRACKINGFROMGRADIENT_EXPORT ttkTrackingFromGradient: public ttkAlgorithm, protected ttk::TrackingFromGradient {

  private:
  vtkSmartPointer<vtkDataSet> PreviousDomain;
  vtkSmartPointer<vtkDataSet> PreviousCriticalPoints;

  public:
  static ttkTrackingFromGradient *New();
  vtkTypeMacro(ttkTrackingFromGradient, ttkAlgorithm);

  protected:

  ttkTrackingFromGradient();
  ~ttkTrackingFromGradient() override;

  int Correlate(vtkImageData* correlations);

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request, vtkInformationVector **inputVector, vtkInformationVector *outputVector) override;
};