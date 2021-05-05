/// \ingroup vtk
/// \class ttkCorrespondenceAlgorithm
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TODO
///
/// TODO
///
/// \sa ttk::CorrespondenceAlgorithm
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkCorrespondenceAlgorithmModule.h>

// VTK Includes
#include <ttkAlgorithm.h>
#include <vtkSmartPointer.h>

class vtkImageData;
class vtkMultiBlockDataSet;
class vtkDataSet;

class TTKCORRESPONDENCEALGORITHM_EXPORT ttkCorrespondenceAlgorithm : public ttkAlgorithm {

private:
  vtkSmartPointer<vtkMultiBlockDataSet> PreviousInputs;

public:
  static ttkCorrespondenceAlgorithm *New();
  vtkTypeMacro(ttkCorrespondenceAlgorithm, ttkAlgorithm);

protected:
  ttkCorrespondenceAlgorithm();
  ~ttkCorrespondenceAlgorithm();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  virtual int Correlate(vtkImageData *correlations,
                vtkDataObject *inputDataObjects0,
                vtkDataObject *inputDataObjects1
  ) {
    return 0;
  };
};
