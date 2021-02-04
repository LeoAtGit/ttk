/// \ingroup vtk
/// \class ttkTrackingFromLabels
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TODO
///
/// TODO
///
/// \sa ttk::TrackingFromLabels
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkTrackingFromLabelsModule.h>

// VTK Includes
#include <ttkAlgorithm.h>
#include <vtkSmartPointer.h>

// TTK Base Includes
#include <TrackingFromLabels.h>

class vtkPolyData;
class vtkDataSet;

class TTKTRACKINGFROMLABELS_EXPORT ttkTrackingFromLabels
  : public ttkAlgorithm,
    protected ttk::TrackingFromLabels {
private:
  vtkSmartPointer<vtkDataSet> previousSegmentation;
  std::vector<vtkSmartPointer<vtkPolyData>> Nodes;
  std::vector<std::vector<int>> AdjacencyMatrices;

public:
  static ttkTrackingFromLabels *New();
  vtkTypeMacro(ttkTrackingFromLabels, ttkAlgorithm);

protected:
  ttkTrackingFromLabels();
  ~ttkTrackingFromLabels() override;

  int Reset();
  int Finalize(vtkPolyData *trackingGraph);

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
