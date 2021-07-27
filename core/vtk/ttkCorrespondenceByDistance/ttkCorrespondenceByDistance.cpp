#include <ttkCorrespondenceByDistance.h>

#include <vtkInformation.h>
#include <vtkObjectFactory.h>

#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointSet.h>
#include <vtkFloatArray.h>

#include <vtkPointData.h>
#include <vtkStringArray.h>

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

  const int nPoints0 = p0->GetNumberOfPoints();
  const int nPoints1 = p1->GetNumberOfPoints();

  // get point coordinates

  auto temp = vtkSmartPointer<vtkFloatArray>::New();
  auto coords0 = nPoints0>0 ? p0->GetPoints()->GetData() : temp;
  auto coords1 = nPoints1>0 ? p1->GetPoints()->GetData() : temp;

  if(coords0->GetDataType() != coords1->GetDataType())
    return !this->printErr("Input vtkPointSet need to have same precision.");

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
  std::string labelIdentifier;
  std::string labelType;
  for(auto &it : std::vector<vtkPointSet *>({p0, p1})) {
    auto labels = this->GetInputArrayToProcess(0, it);
    if(!labels)
      return !this->printErr("Unable to retrieve labels.");
    labelIdentifier = labels->GetName();
    auto array = vtkSmartPointer<vtkDataArray>::Take(labels->NewInstance());
    array->ShallowCopy(labels);
    array->SetName(("IndexLabelMap" + std::to_string(a++)).data());
    fd->AddArray(array);
  }

  auto labelIdentifierArray = vtkSmartPointer<vtkStringArray>::New();
  labelIdentifierArray->SetName("LabelIdentifier");
  labelIdentifierArray->InsertNextValue(labelIdentifier);
  fd->AddArray(labelIdentifierArray);

  return 1;
}