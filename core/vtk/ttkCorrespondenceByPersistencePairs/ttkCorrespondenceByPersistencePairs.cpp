#include <ttkCorrespondenceByPersistencePairs.h>

#include <vtkObjectFactory.h>
#include <vtkInformation.h>

#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointSet.h>

#include <vtkPointData.h>

#include <ttkUtils.h>

vtkStandardNewMacro(ttkCorrespondenceByPersistencePairs);

ttkCorrespondenceByPersistencePairs::ttkCorrespondenceByPersistencePairs() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkCorrespondenceByPersistencePairs::~ttkCorrespondenceByPersistencePairs() {
}

int ttkCorrespondenceByPersistencePairs::Correlate(vtkImageData *correspondences,
                                       vtkDataObject* inputDataObjects0,
                                       vtkDataObject* inputDataObjects1
){
  // unpack input
  auto p0 = vtkPointSet::SafeDownCast(inputDataObjects0); // set of persist. pairs from previous timestep
  auto p1 = vtkPointSet::SafeDownCast(inputDataObjects1); // PPs from current timestep
  if(!p0 || !p1)
    return !this->printErr("Input data objects need to be vtkPointSets.");

  // get point coordinates
  auto coords0 = p0->GetPoints()->GetData();
  auto coords1 = p1->GetPoints()->GetData();

  // get birth & death scalars
  auto birth0 = this->GetInputArrayToProcess(0, p0);
  auto birth1 = this->GetInputArrayToProcess(0, p1);
  auto death0 = this->GetInputArrayToProcess(1, p0);
  auto death1 = this->GetInputArrayToProcess(1, p1);

  if(coords0->GetDataType() != coords1->GetDataType())
    return !this->printErr("Input vtkPointSet need to have same precision.");

  const int nPairs0 = coords0->GetNumberOfTuples() / 2;
  const int nPairs1 = coords1->GetNumberOfTuples() / 2;

  // initialize correspondence matrix i.e., distance matrix
  correspondences->SetDimensions(nPairs0, nPairs1, 1);
  correspondences->AllocateScalars(VTK_DOUBLE, 1); // matching output = double

  auto correspondencesArray = correspondences->GetPointData()->GetArray(0);
  correspondencesArray->SetName("LiftedWassersteinDistance");

  // compute distance matrix
  int status = 0;
  switch(coords0->GetDataType()) {
    vtkTemplateMacro(status = this->computeDistanceMatrix<VTK_TT>(
                       ttkUtils::GetPointer<VTK_TT>(correspondencesArray),
                       ttkUtils::GetPointer<const VTK_TT>(coords0),
                       ttkUtils::GetPointer<const VTK_TT>(coords1), nPairs0,
                       nPairs1));
  }
  if(!status)
    return 0;

  // add index label maps
  int a=0;
  auto fd = correspondences->GetFieldData();
  for(auto& it : std::vector<vtkPointSet*>({p0,p1})){
    auto labels = this->GetInputArrayToProcess(0, it);
    if(!labels)
      return !this->printErr("Unable to retrieve labels.");
    auto array = vtkSmartPointer<vtkDataArray>::Take( labels->NewInstance() );
    array->ShallowCopy(labels);
    array->SetName(("IndexLabelMap"+std::to_string(a++)).data());
    fd->AddArray(array);
  }

  return 1;
}