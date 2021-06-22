#include <ttkBlockAggregator.h>

#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkSmartPointer.h>

#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkMultiBlockDataSet.h>

#include <unordered_set>
#include <vtkThreshold.h>

vtkStandardNewMacro(ttkBlockAggregator);

ttkBlockAggregator::ttkBlockAggregator() {
  this->setDebugMsgPrefix("BlockAggregator");

  this->Reset();

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkBlockAggregator::~ttkBlockAggregator() {
}

int ttkBlockAggregator::FillInputPortInformation(int port,
                                                 vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataObject");
    info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), 1);
    return 1;
  }
  return 0;
}

int ttkBlockAggregator::FillOutputPortInformation(int port,
                                                  vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

int ttkBlockAggregator::Reset() {
  this->AggregatedMultiBlockDataSet
    = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  return 1;
}

int copyObjects(vtkDataObject *source, vtkDataObject *copy) {
  if(source->IsA("vtkMultiBlockDataSet")) {
    auto sourceAsMB = vtkMultiBlockDataSet::SafeDownCast(source);
    auto copyAsMB = vtkMultiBlockDataSet::SafeDownCast(copy);

    if(sourceAsMB == nullptr || copyAsMB == nullptr) {
      return 0;
    }

    const auto sourceFD = sourceAsMB->GetFieldData();
    auto copyFD = copyAsMB->GetFieldData();

    if(sourceFD == nullptr || copyFD == nullptr) {
      return 0;
    }

    copyFD->ShallowCopy(sourceFD);

    for(size_t i = 0; i < sourceAsMB->GetNumberOfBlocks(); i++) {
      auto block = sourceAsMB->GetBlock(i);
      auto blockCopy
        = vtkSmartPointer<vtkDataObject>::Take(block->NewInstance());

      copyObjects(block, blockCopy);
      copyAsMB->SetBlock(i, blockCopy);
    }
  } else {
    copy->ShallowCopy(source);
  }

  return 1;
}

int ttkBlockAggregator::AggregateBlock(vtkDataObject *dataObject) {
  ttk::Timer t;
  size_t nBlocks = this->AggregatedMultiBlockDataSet->GetNumberOfBlocks();
  this->printMsg("Adding object add index " + std::to_string(nBlocks), 0,
                 ttk::debug::LineMode::REPLACE);

  auto copy = vtkSmartPointer<vtkDataObject>::Take(dataObject->NewInstance());
  copyObjects(dataObject, copy);

  this->AggregatedMultiBlockDataSet->SetBlock(nBlocks, copy);

  this->printMsg("Adding object at block index " + std::to_string(nBlocks), 1,
                 t.getElapsedTime());

  return 1;
}

int ttkBlockAggregator::DecomposeDataSet(vtkDataObject *data,
                                         vtkMultiBlockDataSet *collection) {

  if(data->IsA("vtkMultiBlockDataSet")) {
    auto dataAsMB = (vtkMultiBlockDataSet *)data;
    const size_t nBlocks = dataAsMB->GetNumberOfBlocks();
    for(size_t b = 0; b < nBlocks; b++) {
      auto outputAsMB = vtkSmartPointer<vtkMultiBlockDataSet>::New();
      this->DecomposeDataSet(dataAsMB->GetBlock(b), outputAsMB);
      collection->SetBlock(b, outputAsMB);
    }
  } else if(data->IsA("vtkDataSet")) {
    auto array = this->GetInputArrayToProcess(0, data);
    if(!array)
      return !this->printErr("Unable to retrieve input array.");

    if(array->GetNumberOfComponents() != 1)
      return !this->printErr("Input array must be a scalar array.");

    const size_t nTuples = array->GetNumberOfTuples();

    std::vector<double> uniqueValues;
    {
      std::unordered_set<double> uniqueValuesSet;
      for(size_t i = 0; i < nTuples; i++)
        uniqueValuesSet.insert(array->GetTuple1(i));

      uniqueValues.resize(uniqueValuesSet.size());
      size_t i = 0;
      for(auto it : uniqueValuesSet)
        uniqueValues[i++] = it;

      std::sort(uniqueValues.begin(), uniqueValues.end());
    }

    const size_t nUniqueValues = uniqueValues.size();

    auto inputArrayInfo = this->GetInputArrayInformation(0);

    for(size_t i = 0; i < nUniqueValues; i++) {
      auto threshold = vtkSmartPointer<vtkThreshold>::New();
      threshold->SetAllScalars(true);
      threshold->ThresholdBetween(uniqueValues[i], uniqueValues[i]);
      threshold->SetInputArrayToProcess(0, inputArrayInfo);
      threshold->SetInputDataObject(0, data);
      threshold->Update();

      collection->SetBlock(i, threshold->GetOutputDataObject(0));
    }
  } else {
    return !this->printErr("Unsupported Input Data Object Type.");
  }

  return 1;
}

int flattenMB(vtkMultiBlockDataSet *collection, vtkDataObject *data) {
  if(data->IsA("vtkMultiBlockDataSet")) {
    auto dataAsMB = (vtkMultiBlockDataSet *)data;
    const size_t nBlocks = dataAsMB->GetNumberOfBlocks();
    for(size_t b = 0; b < nBlocks; b++)
      flattenMB(collection, dataAsMB->GetBlock(b));
  } else {
    const size_t nBlocks = collection->GetNumberOfBlocks();
    collection->SetBlock(nBlocks, data);
  }

  return 1;
}

int ttkBlockAggregator::RequestData(vtkInformation *request,
                                    vtkInformationVector **inputVector,
                                    vtkInformationVector *outputVector) {

  if(!this->Decompose) {
    // Get iteration information
    double iterationIndex = 0;
    this->SetInputArrayToProcess(1, 0, 0, 2, "_ttk_IterationInfo");
    auto iterationInformation = vtkDoubleArray::SafeDownCast(
      this->GetInputArrayToProcess(1, inputVector));
    if(iterationInformation) {
      iterationIndex = iterationInformation->GetValue(0);
      this->AggregatedMultiBlockDataSet->GetFieldData()->AddArray(
        iterationInformation);
    }

    // Check if AggregatedMultiBlockDataSet needs to be reset
    if(!iterationInformation || this->GetForceReset() || iterationIndex == 0)
      this->Reset();

    // Add all inputs
    const size_t nInputs = inputVector[0]->GetNumberOfInformationObjects();
    for(size_t i = 0; i < nInputs; i++) {
      auto input = vtkDataObject::GetData(inputVector[0], i);

      if(this->GetFlattenInput() && input->IsA("vtkMultiBlockDataSet")) {
        auto inputAsMB = (vtkMultiBlockDataSet *)input;
        const size_t nBlocks = inputAsMB->GetNumberOfBlocks();
        for(size_t j = 0; j < nBlocks; j++)
          this->AggregateBlock(inputAsMB->GetBlock(j));
      } else
        this->AggregateBlock(input);
    }
  } else {
    this->Reset();

    ttk::Timer timer;
    this->printMsg(
      "Decomposing Input into Blocks", 0, ttk::debug::LineMode::REPLACE);

    auto inputAsMB = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    const size_t nInputs = inputVector[0]->GetNumberOfInformationObjects();

    for(size_t i = 0; i < nInputs; i++) {
      auto input = vtkDataObject::GetData(inputVector[0], i);
      if(nInputs == 1 && input->IsA("vtkMultiBlockDataSet"))
        inputAsMB->ShallowCopy(input);
      else
        inputAsMB->SetBlock(i, input);
    }

    auto outputAsMB = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    this->DecomposeDataSet(inputAsMB, outputAsMB);

    if(this->GetFlattenInput()) {
      flattenMB(this->AggregatedMultiBlockDataSet, outputAsMB);
    } else {
      this->AggregatedMultiBlockDataSet->ShallowCopy(outputAsMB);
    }

    this->printMsg("Decomposing Input into Blocks", 1, timer.getElapsedTime());
  }

  // Prepare output
  auto output = vtkMultiBlockDataSet::GetData(outputVector);
  output->ShallowCopy(this->AggregatedMultiBlockDataSet);

  return 1;
}
