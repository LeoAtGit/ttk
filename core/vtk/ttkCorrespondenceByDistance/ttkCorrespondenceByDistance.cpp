#include <ttkCorrespondenceByDistance.h>

#include <vtkInformation.h>
#include <vtkObjectFactory.h>

#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointSet.h>

#include <vtkPointData.h>

#include <ttkUtils.h>

vtkStandardNewMacro(ttkCorrespondenceByDistance);

ttkCorrespondenceByDistance::ttkCorrespondenceByDistance() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkCorrespondenceByDistance::~ttkCorrespondenceByDistance() {
}

int ttkCorrespondenceByDistance::ComputeCorrespondences(
  vtkImageData *correspondenceMatrix,
  vtkDataObject *inputDataObjects0,
  vtkDataObject *inputDataObjects1) {
  // unpack input
  auto p0 = vtkPointSet::SafeDownCast(inputDataObjects0);
  auto p1 = vtkPointSet::SafeDownCast(inputDataObjects1);
  if(!p0 || !p1)
    return !this->printErr("Input data objects need to be vtkPointSets.");

  // get point coordinates
  auto coords0 = p0->GetPoints()->GetData();
  auto coords1 = p1->GetPoints()->GetData();

  if(coords0->GetDataType() != coords1->GetDataType())
    return !this->printErr("Input vtkPointSet need to have same precision.");

  const int nPoints0 = coords0->GetNumberOfTuples();
  const int nPoints1 = coords1->GetNumberOfTuples();

  // initialize correspondence matrix i.e., distance matrix
  correspondenceMatrix->SetDimensions(nPoints0, nPoints1, 1);
  correspondenceMatrix->AllocateScalars(coords0->GetDataType(), 1);
  auto matrixData = correspondenceMatrix->GetPointData()->GetArray(0);
  matrixData->SetName("Distance");

  // compute distance matrix
  int status = 0;
  switch(coords0->GetDataType()) {
    vtkTemplateMacro(status = this->computeDistanceMatrix<VTK_TT>(
                       ttkUtils::GetPointer<VTK_TT>(matrixData),
                       ttkUtils::GetPointer<const VTK_TT>(coords0),
                       ttkUtils::GetPointer<const VTK_TT>(coords1), nPoints0,
                       nPoints1));
  }
  if(!status)
    return 0;

  // add index label maps
  int a = 0;
  auto fd = correspondenceMatrix->GetFieldData();
  for(auto &it : std::vector<vtkPointSet *>({p0, p1})) {
    auto labels = this->GetInputArrayToProcess(0, it);
    if(!labels)
      return !this->printErr("Unable to retrieve labels.");
    auto array = vtkSmartPointer<vtkDataArray>::Take(labels->NewInstance());
    array->ShallowCopy(labels);
    array->SetName(("IndexLabelMap" + std::to_string(a++)).data());
    fd->AddArray(array);
  }

  return 1;
}