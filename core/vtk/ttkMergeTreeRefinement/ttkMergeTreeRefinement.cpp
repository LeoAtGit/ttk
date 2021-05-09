#include <ttkMergeTreeRefinement.h>

#include <vtkDataSet.h>
#include <vtkInformation.h>
#include <vtkUnstructuredGrid.h>

#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>
#include <vtkIntArray.h>
#include <vtkSignedCharArray.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkMergeTreeRefinement);

ttkMergeTreeRefinement::ttkMergeTreeRefinement() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(2);
}

ttkMergeTreeRefinement::~ttkMergeTreeRefinement() {
}

int ttkMergeTreeRefinement::FillInputPortInformation(int port,
                                                     vtkInformation *info) {
  switch(port) {
    case 0:
      info->Set(
        vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
      return 1;
    case 1:
      info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
      return 1;
    default:
      return 0;
  }
}

int ttkMergeTreeRefinement::FillOutputPortInformation(int port,
                                                      vtkInformation *info) {
  switch(port) {
    case 0:
      info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
      return 1;
    case 1:
      info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 1);
      return 1;
    default:
      return 0;
  }
}

int ttkMergeTreeRefinement::RequestData(vtkInformation *request,
                                        vtkInformationVector **inputVector,
                                        vtkInformationVector *outputVector) {
  ttk::Timer timer;
  this->printMsg(
    "Retrieving Input Data", 0, 0, 1, ttk::debug::LineMode::REPLACE);

  // Get the input
  auto i_MergeTree = vtkUnstructuredGrid::GetData(inputVector[0]);
  auto i_Segmentation = vtkDataSet::GetData(inputVector[1]);
  if(!i_MergeTree || !i_Segmentation) {
    this->printErr("Unable to retrieve input data objects.");
    return 0;
  }
  auto o_MergeTree = vtkUnstructuredGrid::GetData(outputVector, 0);
  auto o_Segmentation = vtkDataSet::GetData(outputVector, 1);
  o_Segmentation->ShallowCopy(i_Segmentation);

  const size_t n_i_mtPoints = i_MergeTree->GetNumberOfPoints();
  const size_t n_i_mtEdges = i_MergeTree->GetNumberOfCells();
  const size_t n_i_sPoints = i_Segmentation->GetNumberOfPoints();

  if(n_i_mtEdges < 1) {
    o_MergeTree->ShallowCopy(i_MergeTree);
    return 1;
  }

  // get arrays
  auto i_mtPointData = i_MergeTree->GetPointData();
  auto i_mtCellData = i_MergeTree->GetCellData();

  auto i_mtPointCoords = i_MergeTree->GetPoints()->GetData();
  auto i_mtConnectivity = i_MergeTree->GetCells()->GetConnectivityArray();

  auto i_mtScalars = this->GetInputArrayToProcess(0, i_MergeTree);
  auto i_sScalars = this->GetInputArrayToProcess(0, i_Segmentation);
  if(!i_mtScalars || !i_sScalars) {
    this->printErr("Unable to retrieve input scalar arrays.");
    return 0;
  }
  if(i_mtScalars->GetDataType() != i_sScalars->GetDataType()) {
    this->printErr("Scalar arrays are not of the same type.");
    return 0;
  }

  bool isSplitTree = i_mtScalars->GetTuple1(i_mtConnectivity->GetTuple1(0))
                     > i_mtScalars->GetTuple1(
                       i_mtConnectivity->GetTuple1(n_i_mtEdges * 2 - 1));

  auto i_sBranchId = this->GetInputArrayToProcess(1, i_Segmentation);
  if(!i_sBranchId || i_sBranchId->GetDataType() != VTK_INT) {
    this->printErr("Unable to retrieve BranchId input array of type int*.");
    return 0;
  }
  auto i_sBranchIdData = ttkUtils::GetPointer<int>(i_sBranchId);

  double interval;
  {
    std::string finalExpressionString;

    std::string errorMsg;
    if(!ttkUtils::replaceVariables(this->GetInterval(),
                                   i_MergeTree->GetFieldData(),
                                   finalExpressionString, errorMsg)) {
      this->printErr(errorMsg);
      return 0;
    }

    std::vector<double> values;
    ttkUtils::stringListToDoubleVector(finalExpressionString, values);
    if(values.size() < 1) {
      this->printErr("Unable to parse 'Interval' parameter.");
      return 0;
    }
    interval = values[0];
  }

  this->printMsg("Retrieving Input Data", 1, timer.getElapsedTime(), 1);
  timer.reStart();

  this->printMsg("Refining Merge Tree", 0, 0, 1, ttk::debug::LineMode::REPLACE);

  // compute number of output points and edges
  std::vector<int> edgeNodesOffset(n_i_mtEdges + 1, 0);
  size_t n_o_mtEdges = 0;
  size_t n_o_mtPoints = n_i_mtPoints;
  {
    edgeNodesOffset[0] = n_i_mtPoints;

    for(size_t i = 0; i < n_i_mtEdges; i++) {
      const vtkIdType v0 = i_mtConnectivity->GetTuple1(i * 2);
      const vtkIdType v1 = i_mtConnectivity->GetTuple1(i * 2 + 1);
      const double s0 = i_mtScalars->GetTuple1(v0);
      const double s1 = i_mtScalars->GetTuple1(v1);

      const double delta = std::abs(s0 - s1);

      const size_t nIntervals = std::max(1.0, std::ceil(delta / interval));

      n_o_mtEdges += nIntervals;
      n_o_mtPoints += nIntervals - 1;
      edgeNodesOffset[i + 1] = n_o_mtPoints;
    }
  }

  const auto prepArray
    = [](vtkDataArray *array, const int nTuples, const int nComponents,
         const std::string &name = "") {
        if(name.length() > 0)
          array->SetName(name.data());
        array->SetNumberOfComponents(nComponents);
        array->SetNumberOfTuples(nTuples);
        return ttkUtils::GetVoidPointer(array);
      };

  auto o_mtLambda = vtkSmartPointer<vtkDoubleArray>::New();
  auto o_mtLambdaData
    = static_cast<double *>(prepArray(o_mtLambda, n_o_mtPoints, 1, "Lambda"));

  auto o_mtParentEdge = vtkSmartPointer<vtkIntArray>::New();
  auto o_mtParentEdgeData
    = static_cast<int *>(prepArray(o_mtParentEdge, n_o_mtPoints, 1, "Parent"));

  auto o_mtParentEdge2 = vtkSmartPointer<vtkIntArray>::New();
  auto o_mtParentEdge2Data
    = static_cast<int *>(prepArray(o_mtParentEdge2, n_o_mtEdges, 1, "Parent"));

  auto o_mtNodeId = vtkSmartPointer<vtkIntArray>::New();
  auto o_mtNodeIdData = static_cast<int *>(prepArray(o_mtNodeId, n_o_mtPoints, 1, "NodeId"));
  for(size_t i=0; i<n_o_mtPoints; i++)
    o_mtNodeIdData[i] = i;

  auto nextId = vtkSmartPointer<vtkIntArray>::New();
  auto nextIdData
    = static_cast<int *>(prepArray(nextId, n_o_mtPoints, 1, "NextId"));

  auto size = vtkSmartPointer<vtkIntArray>::New();
  auto sizeData = static_cast<int *>(prepArray(size, n_o_mtPoints, 1, "Size"));

  auto nodeId = vtkSmartPointer<vtkIntArray>::New();
  auto nodeIdData
    = static_cast<int *>(prepArray(nodeId, n_i_sPoints, 1, "NodeId"));

  // edges
  auto o_mtOffsets = vtkSmartPointer<vtkIdTypeArray>::New();
  auto o_mtOffsetsData
    = static_cast<vtkIdType *>(prepArray(o_mtOffsets, n_o_mtEdges + 1, 1));

  auto o_mtConnectivity = vtkSmartPointer<vtkIdTypeArray>::New();
  auto o_mtConnectivityData
    = static_cast<vtkIdType *>(prepArray(o_mtConnectivity, n_o_mtEdges * 2, 1));

  {
    // compute offsets
    for(size_t i = 0; i <= n_o_mtEdges; i++)
      o_mtOffsetsData[i] = i * 2;

    for(size_t i = 0; i < n_o_mtPoints; i++) {
      sizeData[i] = 0;
      nextIdData[i] = -1;
    }

    // compute connectivity
    for(size_t i = 0, c = 0, e = 0; i < n_i_mtEdges; i++) {
      const vtkIdType v0 = i_mtConnectivity->GetTuple1(i * 2);
      const vtkIdType v1 = i_mtConnectivity->GetTuple1(i * 2 + 1);

      const double s0 = i_mtScalars->GetTuple1(v0);
      const double s1 = i_mtScalars->GetTuple1(v1);

      const double delta = std::abs(s0 - s1);

      const auto nIntervals = std::ceil(delta / interval);

      o_mtParentEdgeData[v0] = -1;
      o_mtParentEdgeData[v1] = -1;

      if(nIntervals < 2) {
        o_mtConnectivityData[c++] = v0;
        o_mtConnectivityData[c++] = v1;
        nextIdData[v0] = v1;

        o_mtParentEdge2Data[e++] = i;
      } else {
        auto offset = edgeNodesOffset[i];

        // first edge
        {
          o_mtConnectivityData[c++] = v0;
          o_mtConnectivityData[c++] = offset;
          nextIdData[v0] = offset;

          o_mtLambdaData[offset] = std::abs(delta - interval) / delta;
          o_mtParentEdgeData[offset] = i;

          o_mtParentEdge2Data[e++] = i;
        }

        for(size_t j = 1; j < nIntervals - 1; j++) {
          o_mtConnectivityData[c++] = offset;
          o_mtConnectivityData[c++] = offset + 1;
          nextIdData[offset] = offset + 1;

          offset++;

          o_mtLambdaData[offset] = std::abs(delta - (j + 1) * interval) / delta;
          o_mtParentEdgeData[offset] = i;

          o_mtParentEdge2Data[e++] = i;
        }

        // last edge
        {
          o_mtConnectivityData[c++] = offset;
          o_mtConnectivityData[c++] = v1;
          nextIdData[offset] = v1;

          o_mtParentEdge2Data[e++] = i;
        }
      }
    }
  }

  // compute output point arrays
  std::vector<vtkSmartPointer<vtkDataArray>> outPointArrays;
  {
    std::vector<vtkDataArray *> inArrays;
    inArrays.push_back(i_mtPointCoords);
    for(int i = 0; i < i_mtPointData->GetNumberOfArrays(); i++) {
      auto array = i_mtPointData->GetArray(i);
      if(array && std::string("NextId").compare(array->GetName()) != 0)
        inArrays.push_back(array);
    }
    const size_t nArrays = inArrays.size();
    outPointArrays.resize(nArrays);

    for(size_t a = 0; a < nArrays; a++) {
      const auto inArray = inArrays[a];

      outPointArrays[a]
        = vtkSmartPointer<vtkDataArray>::Take(inArray->NewInstance());
      auto &outArray = outPointArrays[a];
      outArray->SetName(inArray->GetName());
      outArray->SetNumberOfComponents(inArray->GetNumberOfComponents());
      outArray->SetNumberOfTuples(n_o_mtPoints);

      // copy old data
      for(size_t i = 0; i < n_i_mtPoints; i++)
        outArray->SetTuple(i, i, inArray);

      if(std::string("VertexId").compare(outArray->GetName()) == 0){
        for(size_t i = n_i_mtPoints; i < n_o_mtPoints; i++)
          outArray->SetTuple1(i, -1);
        continue;
      } else if(std::string("Type").compare(outArray->GetName()) == 0) {
        for(size_t i = n_i_mtPoints; i < n_o_mtPoints; i++)
          outArray->SetTuple1(i, 2.0);
      } else if(outArray->IsA("vtkDoubleArray")
                || outArray->IsA("vtkFloatArray")) {
        // interpolate new data
        for(size_t i = n_i_mtPoints; i < n_o_mtPoints; i++) {
          const int parentEdge = o_mtParentEdgeData[i];
          const vtkIdType v0 = i_mtConnectivity->GetTuple1(parentEdge * 2);
          const vtkIdType v1 = i_mtConnectivity->GetTuple1(parentEdge * 2 + 1);

          outArray->InterpolateTuple(
            i, v0, inArray, v1, inArray, 1. - o_mtLambdaData[i]);
        }
      } else {
        // snap new data
        for(size_t i = n_i_mtPoints; i < n_o_mtPoints; i++) {
          const int parentEdge = o_mtParentEdgeData[i];
          const vtkIdType v0 = i_mtConnectivity->GetTuple1(parentEdge * 2);
          outArray->SetTuple(i, v0, inArray);
        }
      }
    }
  }

  // compute output cell arrays
  std::vector<vtkSmartPointer<vtkDataArray>> outCellArrays;
  {
    std::vector<vtkDataArray *> inArrays;
    for(int i = 0; i < i_mtCellData->GetNumberOfArrays(); i++) {
      if(auto array = i_mtCellData->GetArray(i))
        inArrays.push_back(array);
    }
    const size_t nArrays = inArrays.size();
    outCellArrays.resize(nArrays);

    for(size_t a = 0; a < nArrays; a++) {
      const auto inArray = inArrays[a];

      outCellArrays[a]
        = vtkSmartPointer<vtkDataArray>::Take(inArray->NewInstance());
      auto &outArray = outCellArrays[a];
      outArray->SetName(inArray->GetName());
      outArray->SetNumberOfComponents(inArray->GetNumberOfComponents());
      outArray->SetNumberOfTuples(n_o_mtEdges);

      for(size_t i = 0; i < n_o_mtEdges; i++)
        outArray->SetTuple(i, o_mtParentEdge2Data[i], inArray);
    }
  }

  // creating refined merge tree output
  {
    auto cells = vtkSmartPointer<vtkCellArray>::New();
    cells->SetData(o_mtOffsets, o_mtConnectivity);
    o_MergeTree->SetCells(VTK_LINE, cells);

    auto points = vtkSmartPointer<vtkPoints>::New();
    points->SetData(outPointArrays[0]);
    o_MergeTree->SetPoints(points);

    auto o_MergeTreePD = o_MergeTree->GetPointData();
    for(size_t a = 1; a < outPointArrays.size(); a++)
      o_MergeTreePD->AddArray(outPointArrays[a]);
    o_MergeTreePD->AddArray(nextId);
    o_MergeTreePD->AddArray(size);
    o_MergeTreePD->AddArray(o_mtNodeId);

    auto o_MergeTreeCD = o_MergeTree->GetCellData();
    for(size_t a = 0; a < outCellArrays.size(); a++)
      o_MergeTreeCD->AddArray(outCellArrays[a]);
  }
  this->printMsg("Refining Merge Tree", 1, timer.getElapsedTime(), 1);

  // compute refined segmentation
  {
    auto o_mtScalars = this->GetInputArrayToProcess(0, o_MergeTree);
    int status = 0;
    ttkTypeMacroA(
      i_sScalars->GetDataType(),
      (
        status = this->computeRefinedSegmentation<T0, ttk::SimplexId>(
        nodeIdData,
        sizeData,

        isSplitTree,
        i_sBranchIdData,
        ttkUtils::GetPointer<const T0>(i_sScalars),
        ttkUtils::GetPointer<const T0>(o_mtScalars),
        nextIdData,
        n_o_mtPoints,
        n_i_sPoints
        )
      )
    );
    if(!status)
      return 0;

    o_Segmentation->GetPointData()->AddArray(nodeId);
  }

  return 1;
}