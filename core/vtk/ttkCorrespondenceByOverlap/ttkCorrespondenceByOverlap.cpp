#include <ttkCorrespondenceByOverlap.h>

#include <vtkObjectFactory.h>

#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkImageData.h>
#include <vtkIntArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointData.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkCorrespondenceByOverlap);

ttkCorrespondenceByOverlap::ttkCorrespondenceByOverlap() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkCorrespondenceByOverlap::~ttkCorrespondenceByOverlap() {
}

int ttkCorrespondenceByOverlap::ComputeCorrespondences(
  vtkImageData *correspondenceMatrix,
  vtkDataObject *inputDataObjects0,
  vtkDataObject *inputDataObjects1) {

  if(this->GetInputArrayAssociation(0, inputDataObjects0) != 0)
    return !this->printErr("Labels must be point data.");

  // get label arrays
  auto labels0 = this->GetInputArrayToProcess(0, inputDataObjects0);
  auto labels1 = this->GetInputArrayToProcess(0, inputDataObjects1);

  // validate arrays
  if(!labels0 || !labels1)
    return !this->printErr("Unable to retrieve labels.");

  if(labels0->GetNumberOfComponents() != 1
     || labels1->GetNumberOfComponents() != 1)
    return !this->printErr("Labels must have exactly one component.");

  if(labels0->GetNumberOfTuples() != labels1->GetNumberOfTuples())
    return !this->printErr("Labels must have same number of values.");

  if(labels0->GetDataType() != labels1->GetDataType())
    return !this->printErr("Labels must have same data type.");

  const int nVertices = labels0->GetNumberOfTuples();

  // extract unique labels from volume
  std::unordered_map<ttk::SimplexId, ttk::SimplexId> labelIndexMap0;
  std::unordered_map<ttk::SimplexId, ttk::SimplexId> labelIndexMap1;
  int status = 0;
  for(auto &it : std::vector<std::pair<
        vtkDataArray *, std::unordered_map<ttk::SimplexId, ttk::SimplexId> *>>(
        {{labels0, &labelIndexMap0}, {labels1, &labelIndexMap1}})) {
    ttkTypeMacroA(
      labels0->GetDataType(),
      (status = this->computeLabelIndexMap<T0, ttk::SimplexId>(
         *it.second, ttkUtils::GetPointer<const T0>(it.first), nVertices)));
    if(!status)
      return 0;
  }

  const int nLabels0 = labelIndexMap0.size();
  const int nLabels1 = labelIndexMap1.size();

  // initialize correspondence matrix
  correspondenceMatrix->SetDimensions(nLabels0, nLabels1, 1);
  correspondenceMatrix->AllocateScalars(VTK_INT, 1);
  auto matrixData = correspondenceMatrix->GetPointData()->GetArray(0);
  matrixData->SetName("Overlap");

  // compute overlaps
  ttkTypeMacroA(labels0->GetDataType(),
                (status = this->computeAdjacencyMatrix<T0, ttk::SimplexId>(
                   ttkUtils::GetPointer<int>(matrixData),
                   ttkUtils::GetPointer<const T0>(labels0),
                   ttkUtils::GetPointer<const T0>(labels1), nVertices,
                   labelIndexMap0, labelIndexMap1)));
  if(!status)
    return 0;

  status = ttkCorrespondenceAlgorithm::AddIndexLabelMaps(
    correspondenceMatrix, labelIndexMap0, labelIndexMap1, labels0->GetName());
  if(!status)
    return 0;

  return 1;
}