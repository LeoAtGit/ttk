#include <ttkCorrespondenceAlgorithm.h>

#include <vtkInformation.h>

#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointSet.h>

#include <vtkDoubleArray.h>
#include <vtkPointData.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkCorrespondenceAlgorithm);

ttkCorrespondenceAlgorithm::ttkCorrespondenceAlgorithm() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);

  this->PreviousInputs = vtkSmartPointer<vtkMultiBlockDataSet>::New();
}

ttkCorrespondenceAlgorithm::~ttkCorrespondenceAlgorithm() {
}

int ttkCorrespondenceAlgorithm::FillInputPortInformation(int port,
                                                         vtkInformation *info) {
  if(port >= 0 && port < this->GetNumberOfInputPorts()) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    info->Append(
      vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

int ttkCorrespondenceAlgorithm::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

template <typename DT>
void BuildLabelIndexMapDT(std::unordered_map<int, int> &map,
                          const int n,
                          const DT *indexLabelMapData) {
  for(int i = 0; i < n; i++)
    map.emplace(std::make_pair(static_cast<int>(indexLabelMapData[i]), i));
}

int ttkCorrespondenceAlgorithm::BuildLabelIndexMap(
  std::unordered_map<int, int> &map, const vtkDataArray *indexLabelMap) {
  if(!indexLabelMap)
    return 0;

  switch(indexLabelMap->GetDataType()) {
    vtkTemplateMacro(BuildLabelIndexMapDT<VTK_TT>(
      map, indexLabelMap->GetNumberOfValues(),
      ttkUtils::GetConstPointer<VTK_TT>(indexLabelMap)));
  }

  return 1;
}

int ttkCorrespondenceAlgorithm::RequestData(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  auto output = vtkMultiBlockDataSet::GetData(outputVector);

  const int nInputs = this->GetNumberOfInputPorts();
  std::vector<vtkDataObject *> inputs(nInputs);
  for(int i = 0; i < nInputs; i++) {
    auto input = vtkDataObject::GetData(inputVector[i]);
    if(!input)
      return !this->printErr("Unable to retrieve input data objects.");
    inputs[i] = input;
  }

  auto firstInput = inputs[0];
  const bool streamingMode = !firstInput->IsA("vtkMultiBlockDataSet");

  std::vector<vtkSmartPointer<vtkMultiBlockDataSet>> sequence;

  if(streamingMode) {
    auto iterationInformation = vtkDoubleArray::SafeDownCast(
      firstInput->GetFieldData()->GetArray("_ttk_IterationInfo"));
    if(this->PreviousInputs->GetNumberOfBlocks() > 0
       && (!iterationInformation || iterationInformation->GetValue(0) > 0)) {
      sequence.resize(2);
      sequence[0] = vtkSmartPointer<vtkMultiBlockDataSet>::New();
      sequence[0]->ShallowCopy(this->PreviousInputs);

      sequence[1] = vtkSmartPointer<vtkMultiBlockDataSet>::New();
      for(int i = 0; i < nInputs; i++)
        sequence[1]->SetBlock(i, inputs[i]);
    }
  } else {
    const int nSteps
      = static_cast<vtkMultiBlockDataSet *>(firstInput)->GetNumberOfBlocks();
    sequence.resize(nSteps);
    for(int s = 0; s < nSteps; s++) {
      sequence[s] = vtkSmartPointer<vtkMultiBlockDataSet>::New();
      for(int i = 0; i < nInputs; i++) {
        sequence[s]->SetBlock(
          i, static_cast<vtkMultiBlockDataSet *>(inputs[i])->GetBlock(s));
      }
    }
  }

  for(size_t s = 1; s < sequence.size(); s++) {
    auto correspondenceMatrix = vtkSmartPointer<vtkImageData>::New();
    auto data0 = sequence[s - 1];
    auto data1 = sequence[s + 0];
    bool singleInput = data0->GetNumberOfBlocks() == 1;

    if(!this->ComputeCorrespondences(correspondenceMatrix,
                                     singleInput ? data0->GetBlock(0) : data0,
                                     singleInput ? data1->GetBlock(0) : data1))
      return 0;

    output->SetBlock(s - 1, correspondenceMatrix);
  }

  if(streamingMode) {
    this->PreviousInputs = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    for(int i = 0; i < nInputs; i++) {
      auto copy = vtkSmartPointer<vtkDataSet>::Take(
        static_cast<vtkDataSet *>(inputs[i])->NewInstance());
      copy->ShallowCopy(inputs[i]);
      this->PreviousInputs->SetBlock(i, copy);
    }
  }

  return 1;
}