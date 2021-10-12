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

#include <unordered_map>

// VTK Module
#include <ttkCorrespondenceAlgorithmModule.h>

// VTK Includes
#include <ttkAlgorithm.h>
#include <vtkSmartPointer.h>

class vtkImageData;
class vtkMultiBlockDataSet;
class vtkDataArray;
class vtkFieldData;

class TTKCORRESPONDENCEALGORITHM_EXPORT ttkCorrespondenceAlgorithm
  : public ttkAlgorithm {

private:
  vtkSmartPointer<vtkMultiBlockDataSet> PreviousInputs;

public:
  static ttkCorrespondenceAlgorithm *New();
  vtkTypeMacro(ttkCorrespondenceAlgorithm, ttkAlgorithm);

  static int GetIndexLabelMaps(vtkDataArray *&indexLabelMapR,
                               vtkDataArray *&indexLabelMapC,
                               vtkFieldData *fieldData);
  static int AddIndexLabelMaps(vtkImageData *correspondenceMatrix,
                               vtkDataArray *indexLabelMapR,
                               vtkDataArray *indexLabelMapC,
                               const std::string& labelIdentifier = "");
  static int AddIndexLabelMaps(
    vtkImageData *correspondenceMatrix,
    const std::unordered_map<ttk::SimplexId, ttk::SimplexId> &labelIndexMapP,
    const std::unordered_map<ttk::SimplexId, ttk::SimplexId> &labelIndexMapC,
    const std::string& labelIdentifier);

  static int
    BuildLabelIndexMap(std::unordered_map<ttk::SimplexId, ttk::SimplexId> &,
                       const vtkDataArray *indexLabelMap);

protected:
  ttkCorrespondenceAlgorithm();
  ~ttkCorrespondenceAlgorithm();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  virtual int ComputeCorrespondences(vtkImageData *vtkNotUsed(correspondenceMatrix),
                                     vtkDataObject *vtkNotUsed(inputDataObjects0),
                                     vtkDataObject *vtkNotUsed(inputDataObjects1)) {
    return 0;
  };
};
