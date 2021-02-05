#include <ttkTrackingFromLabels.h>

#include <vtkInformation.h>

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include <vtkDoubleArray.h>
#include <vtkPointData.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkTrackingFromLabels);

ttkTrackingFromLabels::ttkTrackingFromLabels() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);

  this->Reset();
}

ttkTrackingFromLabels::~ttkTrackingFromLabels() {
}

int ttkTrackingFromLabels::FillInputPortInformation(int port,
                                                    vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    return 1;
  } else if(port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkTrackingFromLabels::FillOutputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
    return 1;
  }
  return 0;
}

int ttkTrackingFromLabels::Reset() {
  this->previousSegmentation = nullptr;
  this->Nodes.clear();
  this->AdjacencyMatrices.clear();
  return 1;
}

template <class ArrayType>
vtkSmartPointer<ArrayType> prepArray(const char *name,
                                     const size_t &nTuples,
                                     const size_t &nComponents) {
  auto array = vtkSmartPointer<ArrayType>::New();
  array->SetName(name);
  array->SetNumberOfComponents(nComponents);
  array->SetNumberOfTuples(nTuples);
  return array;
}

int ttkTrackingFromLabels::Finalize(vtkPolyData *trackingGraph) {

  const size_t nTimesteps = this->Nodes.size();
  std::vector<size_t> offsets(nTimesteps);

  size_t nTGNodes = 0;
  for(size_t t = 0; t < nTimesteps; t++) {
    offsets[t] = nTGNodes;
    nTGNodes += this->Nodes[t]->GetNumberOfPoints();
  }

  // Points
  {
    auto timeArray = prepArray<vtkIntArray>("Time", nTGNodes, 1);
    auto timeArrayData
      = static_cast<int *>(ttkUtils::GetVoidPointer(timeArray));

    auto tgPoints = vtkSmartPointer<vtkPoints>::New();
    tgPoints->SetNumberOfPoints(nTGNodes);
    auto tgPointsData
      = static_cast<float *>(ttkUtils::GetVoidPointer(tgPoints));
    for(size_t t = 0, c = 0, p = 0; t < nTimesteps; t++) {
      const auto &nodes = this->Nodes[t];
      const size_t n = nodes->GetNumberOfPoints();
      auto pointsData = static_cast<const float *>(
        ttkUtils::GetVoidPointer(nodes->GetPoints()));

      for(size_t i = 0; i < n; i++) {
        size_t offset = i * 3;
        tgPointsData[c++] = pointsData[offset++];
        tgPointsData[c++] = pointsData[offset++];
        tgPointsData[c++] = pointsData[offset];

        timeArrayData[p++] = t;
      }
    }

    trackingGraph->SetPoints(tgPoints);
    trackingGraph->GetPointData()->AddArray(timeArray);
  }

  // Cells
  {
    // compute nCells
    size_t nCells = 0;
    for(size_t t = 0; t < nTimesteps - 1; t++) {
      const auto &matrix = this->AdjacencyMatrices[t];
      const size_t nNodes = matrix.size();

      for(size_t i = 0; i < nNodes; i++)
        if(matrix[i] > 0)
          nCells++;
    }

    auto connectivityArray = vtkSmartPointer<vtkIdTypeArray>::New();
    connectivityArray->SetNumberOfTuples(nCells * 2);
    auto connectivityArrayData
      = static_cast<vtkIdType *>(ttkUtils::GetVoidPointer(connectivityArray));
    for(size_t t = 0, q = 0; t < nTimesteps - 1; t++) {
      const auto &matrix = this->AdjacencyMatrices[t];

      const size_t nPNodes = this->Nodes[t]->GetNumberOfPoints();
      const size_t nCNodes = this->Nodes[t + 1]->GetNumberOfPoints();

      for(size_t i = 0; i < nCNodes; i++) {
        for(size_t j = 0; j < nPNodes; j++) {
          if(matrix[i * nPNodes + j] > 0) {
            connectivityArrayData[q++] = offsets[t] + j;
            connectivityArrayData[q++] = offsets[t + 1] + i;
          }
        }
      }
    }

    auto offsetArray = vtkSmartPointer<vtkIdTypeArray>::New();
    offsetArray->SetNumberOfTuples(nCells + 1);
    auto offsetArrayData
      = static_cast<vtkIdType *>(ttkUtils::GetVoidPointer(offsetArray));
    for(size_t i = 0; i <= nCells; i++)
      offsetArrayData[i] = i * 2;

    auto cellArray = vtkSmartPointer<vtkCellArray>::New();
    cellArray->SetData(offsetArray, connectivityArray);

    trackingGraph->SetLines(cellArray);
  }

  return 1;
}

int ttkTrackingFromLabels::RequestData(vtkInformation *request,
                                       vtkInformationVector **inputVector,
                                       vtkInformationVector *outputVector) {

  // Get iteration information
  this->SetInputArrayToProcess(1, 0, 0, 2, "_ttk_IterationInfo");
  auto iterationInformation = vtkDoubleArray::SafeDownCast(
    this->GetInputArrayToProcess(1, inputVector));
  if(!iterationInformation)
    return 0;

  double iterationIndex = iterationInformation->GetValue(0);
  double iterationLimit = iterationInformation->GetValue(1);

  auto currentNodes = vtkPolyData::GetData(inputVector[0]);
  if(!currentNodes)
    return 0;

  auto currentSegmentation = vtkDataSet::GetData(inputVector[1]);
  if(!currentSegmentation)
    return 0;

  // Check if AggregatedMultiBlockDataSet needs to be reset
  if(iterationIndex == 0) {
    this->printMsg("First Iteration");
    this->Reset();
  } else {
    this->printMsg("Iteration " + std::to_string(iterationIndex));

    auto cLabels = this->GetInputArrayToProcess(0, currentSegmentation);
    if(!cLabels) {
      this->printErr("Unable to retrieve input array.");
      return 0;
    }
    if(this->GetInputArrayAssociation(0, inputVector) != 0) {
      this->printErr("Input array needs to be a point data array.");
      return 0;
    }
    if(cLabels->GetNumberOfComponents() != 1) {
      this->printErr("Input array needs to be a scalar array.");
      return 0;
    }

    auto pLabels = this->GetInputArrayToProcess(0, this->previousSegmentation);
    if(!pLabels) {
      this->printErr("Unable to retrieve input array.");
      return 0;
    }
    if(pLabels->GetDataType() != cLabels->GetDataType()) {
      this->printErr("Previous and current label array must be of same type.");
      return 0;
    }
    if(pLabels->GetNumberOfTuples() != cLabels->GetNumberOfTuples()) {
      this->printErr(
        "Previous and current label array must be of same length.");
      return 0;
    }

    this->AdjacencyMatrices.resize(iterationIndex);

    const size_t nPNodes = this->Nodes[iterationIndex - 1]->GetNumberOfPoints();
    const size_t nCNodes = currentNodes->GetNumberOfPoints();

    auto &adjacencyMatrix = this->AdjacencyMatrices[iterationIndex - 1];
    adjacencyMatrix.resize(nPNodes * nCNodes, 0);

    int status = 0;
    switch(cLabels->GetDataType()) {
      vtkTemplateMacro(
        status = this->computeEdges<VTK_TT>(
          adjacencyMatrix,
          static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(pLabels)),
          static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(cLabels)),
          cLabels->GetNumberOfTuples(), nPNodes, nCNodes));
    }
    if(!status)
      return 0;
  }

  this->Nodes.resize(iterationIndex + 1);
  this->Nodes[iterationIndex] = vtkSmartPointer<vtkPolyData>::New();
  this->Nodes[iterationIndex]->ShallowCopy(currentNodes);

  this->previousSegmentation
    = vtkSmartPointer<vtkDataSet>::Take(currentSegmentation->NewInstance());
  this->previousSegmentation->ShallowCopy(currentSegmentation);

  if(iterationIndex == iterationLimit - 1) {
    if(!this->Finalize(vtkPolyData::GetData(outputVector)))
      return 0;
  }

  return 1;
}
