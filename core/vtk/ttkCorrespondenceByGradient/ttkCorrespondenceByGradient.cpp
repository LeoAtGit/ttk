#include <ttkCorrespondenceByGradient.h>

#include <vtkDataSet.h>
#include <vtkInformation.h>
#include <vtkMultiBlockDataSet.h>

#include <vtkImageData.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkCorrespondenceByGradient);

ttkCorrespondenceByGradient::ttkCorrespondenceByGradient() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}
ttkCorrespondenceByGradient::~ttkCorrespondenceByGradient() {
}

vtkIntArray *GetVertexIdArray(vtkDataSet *input) {
  return vtkIntArray::SafeDownCast(
    input->GetPointData()->GetArray("ttkVertexScalarField"));
}

int ttkCorrespondenceByGradient::Correlate(
  vtkImageData *correspondences,
  vtkDataObject *inputDataObjects0,
  vtkDataObject *inputDataObjects1
) {
  auto inputs0AsMB = static_cast<vtkMultiBlockDataSet*>(inputDataObjects0);
  auto inputs1AsMB = static_cast<vtkMultiBlockDataSet*>(inputDataObjects1);

  auto domain0 = vtkDataSet::SafeDownCast(inputs0AsMB->GetBlock(0));
  auto domain1 = vtkDataSet::SafeDownCast(inputs1AsMB->GetBlock(0));
  auto seeds0 = vtkDataSet::SafeDownCast(inputs0AsMB->GetBlock(1));
  auto seeds1 = vtkDataSet::SafeDownCast(inputs1AsMB->GetBlock(1));

  int nFeatures0 = seeds0->GetNumberOfElements(0);
  int nFeatures1 = seeds1->GetNumberOfElements(0);

  // allocate correlation matrices
  correspondences->SetDimensions(nFeatures0, nFeatures1, 1);
  correspondences->AllocateScalars(VTK_INT, 1);

  auto forward = correspondences->GetPointData()->GetArray(0);
  forward->SetName("Forward");

  auto backward = vtkSmartPointer<vtkIntArray>::New();
  backward->DeepCopy(forward);
  backward->SetName("Backward");
  correspondences->GetPointData()->AddArray(backward);

  auto orderArray0 = ttkAlgorithm::GetOrderArray(domain0, 0);
  auto orderArray1 = ttkAlgorithm::GetOrderArray(domain1, 0);

  int status = 0;
  auto triangulation = this->GetTriangulation(domain0);
  this->preconditionTriangulation(triangulation);

  ttkVtkTemplateMacro(
    orderArray1->GetDataType(), triangulation->getType(),
    (status = this->computeCorrespondences<ttk::SimplexId, TTK_TT>(
       ttkUtils::GetPointer<int>(forward),

       ttkUtils::GetPointer<ttk::SimplexId>(orderArray1),
       static_cast<TTK_TT *>(triangulation->getData()),

       ttkUtils::GetPointer<ttk::SimplexId>(GetVertexIdArray(seeds0)),
       ttkUtils::GetPointer<ttk::SimplexId>(GetVertexIdArray(seeds1)),
       nFeatures0, nFeatures1, std::greater<ttk::SimplexId>{},
       [](ttk::SimplexId i, ttk::SimplexId j, ttk::SimplexId n,
          ttk::SimplexId m) { return j * n + i; })));
  if(!status)
    return 0;

  ttkVtkTemplateMacro(
    orderArray0->GetDataType(), triangulation->getType(),
    (status = this->computeCorrespondences<ttk::SimplexId, TTK_TT>(
       ttkUtils::GetPointer<int>(backward),

       ttkUtils::GetPointer<ttk::SimplexId>(orderArray0),
       static_cast<TTK_TT *>(triangulation->getData()),

       ttkUtils::GetPointer<ttk::SimplexId>(GetVertexIdArray(seeds1)),
       ttkUtils::GetPointer<ttk::SimplexId>(GetVertexIdArray(seeds0)),
       nFeatures1, nFeatures0, std::greater<ttk::SimplexId>{},
       [](ttk::SimplexId i, ttk::SimplexId j, ttk::SimplexId n,
          ttk::SimplexId m) { return i * m + j; })));
  if(!status)
    return 0;

  // add index label maps
  int a=0;
  auto fd = correspondences->GetFieldData();
  for(auto& it : std::vector<vtkDataSet*>({seeds0,seeds1})){
    auto labels = this->GetInputArrayToProcess(1, it);
    if(!labels)
      return !this->printErr("Unable to retrieve labels.");
    auto array = vtkSmartPointer<vtkDataArray>::Take( labels->NewInstance() );
    array->ShallowCopy(labels);
    array->SetName(("IndexLabelMap"+std::to_string(a++)).data());
    fd->AddArray(array);
  }

  return 1;
}