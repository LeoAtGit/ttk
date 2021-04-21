/// \ingroup vtk
/// \class ttkTrackingFromOverlap
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TODO
///
/// TODO
///
/// \sa ttk::TrackingFromOverlap
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkTrackingFromOverlapModule.h>

// VTK Includes
#include <ttkAlgorithm.h>
#include <vtkSmartPointer.h>

// TTK Base Includes
#include <TrackingFromOverlap.h>

class vtkDataArray;
class vtkImageData;
class vtkDataSet;

class TTKTRACKINGFROMOVERLAP_EXPORT ttkTrackingFromOverlap
  : public ttkAlgorithm,
    protected ttk::TrackingFromOverlap {
private:
  vtkSmartPointer<vtkDataSet> PreviousInput;

public:
  static ttkTrackingFromOverlap *New();
  vtkTypeMacro(ttkTrackingFromOverlap, ttkAlgorithm);

protected:
  ttkTrackingFromOverlap();
  ~ttkTrackingFromOverlap() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  int Correlate(vtkImageData* correlations, vtkDataArray* l0, vtkDataArray* l1);
};
