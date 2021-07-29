#include <ttkTrackingGraph.h>

#include <vtkInformation.h>

#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkStringArray.h>

#include <ttkCorrespondenceAlgorithm.h>
#include <ttkMacros.h>
#include <ttkUtils.h>

typedef std::unordered_map<std::string, vtkAbstractArray *> ArrayMap;
typedef std::vector<std::tuple<vtkFieldData *, ArrayMap *, int>> ArrayMapSet;

vtkStandardNewMacro(ttkTrackingGraph);

ttkTrackingGraph::ttkTrackingGraph() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

ttkTrackingGraph::~ttkTrackingGraph() {
}

int ttkTrackingGraph::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0 || port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    if(port == 1)
      info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
    return 1;
  }
  return 0;
}

int ttkTrackingGraph::FillOutputPortInformation(int port,
                                                vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
    return 1;
  }
  return 0;
}

template <typename CT, typename IT>
int countEdges(int &nEdges,
               const int *dim,
               const CT *correspondences,
               const int nLabels0 = 0,
               const int nLables1 = 0,
               const IT *labels0 = nullptr,
               const IT *labels1 = nullptr,
               const vtkDataArray *indexLabelMap0 = nullptr,
               const vtkDataArray *indexLabelMap1 = nullptr) {

  if(nLabels0 > 0 && nLables1 > 0) {
    std::unordered_map<int, int> labelIndexMap0;
    std::unordered_map<int, int> labelIndexMap1;
    ttkCorrespondenceAlgorithm::BuildLabelIndexMap(
      labelIndexMap0, indexLabelMap0);
    ttkCorrespondenceAlgorithm::BuildLabelIndexMap(
      labelIndexMap1, indexLabelMap1);

    for(int i = 0; i < nLabels0; i++) {
      for(int j = 0; j < nLables1; j++) {
        const auto &iLabel = labels0[i];
        const auto &jLabel = labels1[j];
        const auto &iIt = labelIndexMap0.find(iLabel);
        const auto &jIt = labelIndexMap1.find(jLabel);
        if(iIt == labelIndexMap0.end() || jIt == labelIndexMap1.end())
          continue;

        const auto &iIdx = iIt->second;
        const auto &jIdx = jIt->second;
        const int index = jIdx * dim[0] + iIdx;
        if(correspondences[index] > 0)
          nEdges++;
      }
    }
  } else {
    for(int i = 0; i < dim[0]; i++)
      for(int j = 0; j < dim[1]; j++) {
        const auto &iIdx = i;
        const auto &jIdx = j;
        const int index = jIdx * dim[0] + iIdx;
        if(correspondences[index] > 0)
          nEdges++;
      }
  }

  return 1;
}

template <typename CT, typename IT>
int generateEdges(vtkPolyData *output,
                  int &iEdge,
                  vtkFieldData *trackingGraphCD,
                  vtkFieldData *correspondencesPD,
                  const CT *correspondences,
                  const int &offset0,
                  const int &offset1,
                  const int *dim,
                  const int nLabels0 = 0,
                  const int nLables1 = 0,
                  const IT *labels0 = nullptr,
                  const IT *labels1 = nullptr,
                  const vtkDataArray *indexLabelMap0 = nullptr,
                  const vtkDataArray *indexLabelMap1 = nullptr) {

  std::vector<std::pair<vtkAbstractArray *, vtkAbstractArray *>> arrayMap;
  for(int a = 0; a < trackingGraphCD->GetNumberOfArrays(); a++) {
    auto oArray = trackingGraphCD->GetAbstractArray(a);
    arrayMap.push_back(
      {oArray, correspondencesPD->GetArray(oArray->GetName())});
  }

  if(nLabels0 > 0 && nLables1 > 0) {
    std::unordered_map<int, int> labelIndexMap0;
    std::unordered_map<int, int> labelIndexMap1;
    ttkCorrespondenceAlgorithm::BuildLabelIndexMap(
      labelIndexMap0, indexLabelMap0);
    ttkCorrespondenceAlgorithm::BuildLabelIndexMap(
      labelIndexMap1, indexLabelMap1);

    for(int i = 0; i < nLabels0; i++) {
      for(int j = 0; j < nLables1; j++) {
        const auto &iLabel = labels0[i];
        const auto &jLabel = labels1[j];
        const auto &iIt = labelIndexMap0.find(iLabel);
        const auto &jIt = labelIndexMap1.find(jLabel);
        if(iIt == labelIndexMap0.end() || jIt == labelIndexMap1.end())
          continue;

        const auto &iIdx = iIt->second;
        const auto &jIdx = jIt->second;
        const int index = jIdx * dim[0] + iIdx;
        if(correspondences[index] <= 0)
          continue;

        vtkIdType points[2]{offset0 + i, offset1 + j};
        output->InsertNextCell(VTK_LINE, 2, points);

        for(auto &it : arrayMap)
          it.first->SetTuple(iEdge, index, it.second);

        iEdge++;
      }
    }
  } else {
    for(int i = 0; i < dim[0]; i++) {
      for(int j = 0; j < dim[1]; j++) {
        const auto &iIdx = i;
        const auto &jIdx = j;
        const int index = jIdx * dim[0] + iIdx;
        if(correspondences[index] <= 0)
          continue;

        const vtkIdType points[2]{offset0 + i, offset1 + j};
        output->InsertNextCell(VTK_LINE, 2, points);

        for(auto &it : arrayMap)
          it.first->SetTuple(iEdge, index, it.second);

        iEdge++;
      }
    }
  }

  return 1;
}

int ttkTrackingGraph::CountNodesAndEdges(int &nNodes,
                                         int &nEdges,
                                         std::vector<int> &nodeIdxOffsets,
                                         vtkMultiBlockDataSet *correspondences,
                                         vtkMultiBlockDataSet *features
                                         = nullptr) {

  const int nSteps = correspondences->GetNumberOfBlocks();

  nNodes = 0;
  nEdges = 0;
  nodeIdxOffsets.resize(nSteps + 1, 0);
  nodeIdxOffsets[0] = 0;

  if(features) {
    // nodes
    for(int t = 0; t <= nSteps; t++) {
      auto f = static_cast<vtkPointSet *>(features->GetBlock(t));
      const auto n = f->GetNumberOfPoints();
      nNodes += n;
      if(t < nSteps)
        nodeIdxOffsets[t + 1] = n + nodeIdxOffsets[t];
    }

    // edges
    for(int t = 0; t < nSteps; t++) {
      auto f0 = static_cast<vtkPointSet *>(features->GetBlock(t));
      auto f1 = static_cast<vtkPointSet *>(features->GetBlock(t + 1));
      auto l0 = this->GetInputArrayToProcess(1, f0);
      auto l1 = this->GetInputArrayToProcess(1, f1);

      auto c = static_cast<vtkImageData *>(correspondences->GetBlock(t));
      int dim[3];
      c->GetDimensions(dim);
      auto matrix = this->GetInputArrayToProcess(0, c);
      auto indexLabelMap0 = c->GetFieldData()->GetArray("IndexLabelMap0");
      auto indexLabelMap1 = c->GetFieldData()->GetArray("IndexLabelMap1");

      ttkTypeMacroAI(
        matrix->GetDataType(), l0->GetDataType(),
        countEdges(nEdges, dim, ttkUtils::GetPointer<T0>(matrix),
                   l0->GetNumberOfValues(), l1->GetNumberOfValues(),
                   ttkUtils::GetPointer<T1>(l0), ttkUtils::GetPointer<T1>(l1),
                   indexLabelMap0, indexLabelMap1));
    }
  } else {
    for(int t = 0; t < nSteps; t++) {
      auto c = static_cast<vtkImageData *>(correspondences->GetBlock(t));
      int dim[3];
      c->GetDimensions(dim);

      // nodes
      if(t == 0)
        nNodes += dim[0];
      nNodes += dim[1];

      nodeIdxOffsets[t + 1] = dim[0] + nodeIdxOffsets[t];

      // edges
      auto matrix = this->GetInputArrayToProcess(0, c);
      switch(matrix->GetDataType()) {
        vtkTemplateMacro((countEdges<VTK_TT, int>(
          nEdges, dim, ttkUtils::GetPointer<VTK_TT>(matrix))));
      }
    }
  }

  return 1;
};

int collectArrays(ArrayMap &arrayMap,
                  vtkMultiBlockDataSet *correspondences,
                  int attribute) {
  const int nSteps = correspondences->GetNumberOfBlocks();
  for(int t = 0; t < nSteps; t++) {
    auto c = static_cast<vtkImageData *>(correspondences->GetBlock(t));
    auto data = c->GetAttributesAsFieldData(attribute);
    for(int a = 0; a < data->GetNumberOfArrays(); a++) {
      auto array = data->GetAbstractArray(a);
      arrayMap.insert({array->GetName(), array});
    }
  }

  return 1;
}

// int ttkTrackingGraph::GenerateSpatialTrackingGraph(
//   vtkPolyData *output,
//   vtkMultiBlockDataSet *correspondences,
//   vtkMultiBlockDataSet *features) {
//   ttk::Timer timer;
//   std::string msg = "Initializing Output";
//   this->printMsg(msg, 0, 0, ttk::debug::LineMode::REPLACE);

//   const int nSteps = features->GetNumberOfBlocks();
//   if(nSteps != (int)correspondences->GetNumberOfBlocks() + 1)
//     return !this->printErr("Incompatible vtkMultiBlockDataSet structures for
//     "
//                           "features and correspondences.");

//   auto outputPD = output->GetPointData();
//   auto trackingGraphCD = output->GetCellData();

//   typedef std::unordered_map<std::string, vtkAbstractArray *> ArrayMap;
//   typedef std::vector<std::tuple<vtkFieldData *, ArrayMap *, int>>
//   ArrayMapSet;

//   int nNodes = 0;
//   int nEdges = 0;
//   std::vector<int> pointIndexOffset(1, 0);
//   {
//     ArrayMap iPointDataMap;
//     ArrayMap iCellDataMap;
//     ArrayMap iFieldDataMap;

//     for(int t = 0; t < nSteps; t++) {
//       auto featuresC = vtkPointSet::SafeDownCast(features->GetBlock(t));
//       if(!featuresC)
//         return !this->printErr("Features must be a list of vtkPointSets.");

//       const int nFeaturesC = featuresC->GetNumberOfPoints();
//       pointIndexOffset.push_back(
//         t < 1 ? nFeaturesC
//               : pointIndexOffset[pointIndexOffset.size() - 1] + nFeaturesC);
//       nNodes += nFeaturesC;

//       ArrayMapSet iSet;
//       for(auto it :
//           ArrayMapSet({{featuresC->GetPointData(), &iPointDataMap, 0},
//                       {featuresC->GetFieldData(), &iFieldDataMap, 0}})) {
//         for(int a = 0; a < std::get<0>(it)->GetNumberOfArrays(); a++) {
//           auto array = std::get<0>(it)->GetArray(a);
//           if(array)
//             std::get<1>(it)->insert(std::make_pair(array->GetName(), array));
//         }
//       }

//       if(t == 0)
//         continue;

//       auto correspondencesPC
//         = vtkImageData::SafeDownCast(correspondences->GetBlock(t - 1));
//       if(!correspondencesPC)
//         return !this->printErr(
//           "Correspondences must be a list of vtkImageData objects.");

//       int dim[3];
//       correspondencesPC->GetDimensions(dim);
//       if(dim[0] < 1 || dim[1] < 1)
//         continue;

//       auto matrix = this->GetInputArrayToProcess(0, correspondencesPC);
//       if(!matrix)
//         return !this->printErr("Unable to retrieve correspondence matrix.");

//       auto featuresP = vtkPointSet::SafeDownCast(features->GetBlock(t - 1));
//       const int nFeaturesP = featuresP->GetNumberOfPoints();
//       if(nFeaturesP < 1)
//         continue;

//       auto labelsP = this->GetInputArrayToProcess(1, featuresP);
//       auto labelsC = this->GetInputArrayToProcess(1, featuresC);
//       if(!labelsP || !labelsC)
//         return !this->printErr("Label lookup requires label arrays.");

//       if(labelsP->GetDataType() != labelsC->GetDataType())
//         return !this->printErr("Labels must be of same data type.");

//       auto cFD = correspondencesPC->GetFieldData();
//       auto indexLabelMap0 = cFD->GetArray("IndexLabelMap0");
//       auto indexLabelMap1 = cFD->GetArray("IndexLabelMap1");
//       if(!indexLabelMap0 || !indexLabelMap1)
//         return !this->printErr(
//           "Label lookup requires IndexLabelMaps for Correspondence
//           Matrices.");

//       switch(
//         vtkTemplate2PackMacro(labelsC->GetDataType(), matrix->GetDataType()))
//         { ttkTemplate2IdMacro(countEdges(
//           nEdges, ttkUtils::GetPointer<VTK_T2>(matrix), dim, nFeaturesP,
//           nFeaturesC, ttkUtils::GetPointer<VTK_T1>(labelsP),
//           ttkUtils::GetPointer<VTK_T1>(labelsC),
//           ttkUtils::GetPointer<int>(indexLabelMap0),
//           ttkUtils::GetPointer<int>(indexLabelMap1)));
//       }

//       auto data = correspondencesPC->GetPointData();
//       for(int a = 0; a < data->GetNumberOfArrays(); a++) {
//         auto iArray = data->GetArray(a);
//         if(iArray)
//           iCellDataMap.insert(std::make_pair(iArray->GetName(), iArray));
//       }
//     }

//     ArrayMapSet oSet{{outputPD, &iPointDataMap, nNodes},
//                     {outputPD, &iFieldDataMap, nNodes},
//                     {trackingGraphCD, &iCellDataMap, nEdges}};
//     for(auto it : oSet) {
//       for(auto ait : *(std::get<1>(it))) {
//         auto array
//           =
//           vtkSmartPointer<vtkAbstractArray>::Take(ait.second->NewInstance());
//         array->SetName(ait.second->GetName());
//         array->SetNumberOfComponents(ait.second->GetNumberOfComponents());
//         array->SetNumberOfTuples(std::get<2>(it));
//         std::get<0>(it)->AddArray(array);
//       }
//     }
//   }

//   auto points = vtkSmartPointer<vtkPoints>::New();
//   points->SetDataTypeToFloat();
//   points->SetNumberOfPoints(nNodes);
//   auto pointCoords = points->GetData();

//   output->SetPoints(points);
//   output->AllocateExact(0, 0, nEdges, 2, 0, 0, 0, 0);

//   this->printMsg(msg, 1, timer.getElapsedTime());
//   timer.reStart();
//   msg = "Generating Spatial Tracking Graph";
//   this->printMsg(msg, 0, 0, ttk::debug::LineMode::REPLACE);

//   const int nPointArrays = outputPD->GetNumberOfArrays();
//   int iEdge = 0;
//   for(int t = 0, q = 0; t < nSteps; t++) {
//     auto featuresC = vtkPointSet::SafeDownCast(features->GetBlock(t));
//     const int nFeaturesC = featuresC->GetNumberOfPoints();
//     if(nFeaturesC < 1)
//       continue;

//     auto pointData = featuresC->GetPointData();
//     auto fieldData = featuresC->GetFieldData();

//     pointCoords->InsertTuples(
//       q, nFeaturesC, 0, featuresC->GetPoints()->GetData());

//     for(int a = 0; a < nPointArrays; a++) {
//       auto oArray = outputPD->GetAbstractArray(a);
//       auto iArray = pointData->GetAbstractArray(oArray->GetName());
//       if(iArray) {
//         oArray->InsertTuples(q, nFeaturesC, 0, iArray);
//       } else {
//         iArray = fieldData->GetAbstractArray(oArray->GetName());
//         for(int i = 0; i < nFeaturesC; i++)
//           oArray->InsertTuples(q + i, 1, 0, iArray);
//       }
//     }

//     q += nFeaturesC;

//     if(t == 0)
//       continue;

//     auto correspondencesPC
//       = vtkImageData::SafeDownCast(correspondences->GetBlock(t - 1));
//     int dim[3];
//     correspondencesPC->GetDimensions(dim);

//     auto featuresP = vtkPointSet::SafeDownCast(features->GetBlock(t - 1));
//     const int nFeaturesP = featuresP->GetNumberOfPoints();
//     if(nFeaturesP < 1 || nFeaturesC < 1)
//       continue;

//     auto matrix = this->GetInputArrayToProcess(0, correspondencesPC);
//     auto labelsP = this->GetInputArrayToProcess(1, featuresP);
//     auto labelsC = this->GetInputArrayToProcess(1, featuresC);
//     auto cFD = correspondencesPC->GetFieldData();
//     auto indexLabelMap0 = cFD->GetArray("IndexLabelMap0");
//     auto indexLabelMap1 = cFD->GetArray("IndexLabelMap1");

//     switch(
//       vtkTemplate2PackMacro(labelsC->GetDataType(), matrix->GetDataType())) {
//       ttkTemplate2IdMacro((generateEdges<VTK_T2, VTK_T1>(
//         output, iEdge, trackingGraphCD, correspondencesPC->GetPointData(),
//         ttkUtils::GetPointer<VTK_T2>(matrix), pointIndexOffset[t - 1],
//         pointIndexOffset[t + 0], dim, nFeaturesP, nFeaturesC,
//         ttkUtils::GetPointer<VTK_T1>(labelsP),
//         ttkUtils::GetPointer<VTK_T1>(labelsC),
//         ttkUtils::GetPointer<int>(indexLabelMap0),
//         ttkUtils::GetPointer<int>(indexLabelMap1))));
//     }
//   }

//   this->printMsg(msg, 1, timer.getElapsedTime());
//   return 1;
// }

int ttkTrackingGraph::Validate(vtkMultiBlockDataSet *correspondences,
                               vtkMultiBlockDataSet *features) {
  ttk::Timer timer;
  const std::string msg = "Validating Input";
  this->printMsg(msg, 0, 0, ttk::debug::LineMode::REPLACE);

  const std::string errmsg0
    = "Correspondences input must be a flat vtkMultiBlockDataSet that contains "
      "only vtkImageData objects.";
  const std::string errmsg1 = "Features must be vtkPointSet where each feature "
                              "is represented by a single point.";

  if(!correspondences)
    return !this->printErr(errmsg0);

  const int nSteps = correspondences->GetNumberOfBlocks();
  for(int t = 0; t < nSteps; t++) {
    auto c = correspondences->GetBlock(t);
    if(!c || !c->IsA("vtkImageData"))
      return !this->printErr(errmsg0);
    auto matrix = this->GetInputArrayToProcess(0, c);
    if(!matrix)
      return !this->printErr("Unable to retrieve correspondence matrix.");
  }

  if(features) {
    if(features->GetNumberOfBlocks() - 1
       != correspondences->GetNumberOfBlocks())
      return !this->printErr(
        "Number of feature sets (F) and number of correspondence matrices (C) "
        "must satisfy F-1 = C.");

    auto ct = vtkSmartPointer<vtkCellTypes>::New();

    for(int t = 0; t <= nSteps; t++) {
      auto f = vtkPointSet::SafeDownCast(features->GetBlock(t));
      if(!f)
        return !this->printErr(errmsg1);
      f->GetCellTypes(ct);
      if(ct->GetNumberOfTypes() != 1 || !ct->IsType(VTK_VERTEX))
        return !this->printErr(errmsg1);
      auto labels = this->GetInputArrayToProcess(1, f);
      if(!labels)
        return !this->printErr("Unable to retrieve feature labels.");
    }
  }

  this->printMsg(msg, 1, timer.getElapsedTime());

  return 1;
}

int ttkTrackingGraph::GenerateTrackingGraphFromFeatures(
  vtkPolyData *trackingGraph,
  vtkMultiBlockDataSet *correspondences,
  vtkMultiBlockDataSet *features) {

  // initialize trackingGraph
  ttk::Timer timer;
  std::string msg = "Initializing Output";
  this->printMsg(msg, 0, 0, ttk::debug::LineMode::REPLACE);

  const int nSteps = correspondences->GetNumberOfBlocks();

  int nNodes, nEdges;
  std::vector<int> nodeIdxOffsets;
  if(!this->CountNodesAndEdges(
       nNodes, nEdges, nodeIdxOffsets, correspondences, features))
    return 0;

  ArrayMap featuresPD;
  collectArrays(featuresPD, features, 0);
  ArrayMap featuresFD;
  collectArrays(featuresFD, features, 2);
  ArrayMap correspondencesPD;
  collectArrays(correspondencesPD, correspondences, 0);

  // allocating memory
  auto trackingGraphPD = trackingGraph->GetPointData();
  auto trackingGraphCD = trackingGraph->GetCellData();

  auto points = vtkSmartPointer<vtkPoints>::New();
  points->SetDataTypeToFloat();
  points->SetNumberOfPoints(nNodes);
  trackingGraph->SetPoints(points);
  trackingGraph->AllocateExact(0, 0, nEdges, 2, 0, 0, 0, 0);

  {
    ArrayMapSet oSet{{trackingGraphPD, &featuresPD, nNodes},
                     {trackingGraphPD, &featuresFD, nNodes},
                     {trackingGraphCD, &correspondencesPD, nEdges}};
    for(auto it : oSet) {
      for(auto ait : *(std::get<1>(it))) {
        auto array
          = vtkSmartPointer<vtkAbstractArray>::Take(ait.second->NewInstance());
        array->SetName(ait.second->GetName());
        array->SetNumberOfComponents(ait.second->GetNumberOfComponents());
        array->SetNumberOfTuples(std::get<2>(it));
        std::get<0>(it)->AddArray(array);
      }
    }
  }

  this->printMsg(msg, 1, timer.getElapsedTime());
  timer.reStart();

  // ---------------------------------------------------------------------------
  msg = "Generating Spatial Tracking Graph";
  this->printMsg(msg, 0, 0, ttk::debug::LineMode::REPLACE);

  // Nodes
  {
    auto pointCoords = points->GetData();
    for(int t = 0, q = 0; t <= nSteps; t++) {
      auto f = vtkPointSet::SafeDownCast(features->GetBlock(t));
      const int n = f->GetNumberOfPoints();
      if(n < 1)
        continue;
      pointCoords->InsertTuples(q, n, 0, f->GetPoints()->GetData());

      auto fPD = f->GetPointData();
      auto fFD = f->GetFieldData();
      for(int a = 0; a < trackingGraphPD->GetNumberOfArrays(); a++) {
        auto oArray = trackingGraphPD->GetAbstractArray(a);
        auto iArray = fPD->GetAbstractArray(oArray->GetName());
        if(iArray) {
          oArray->InsertTuples(q, n, 0, iArray);
        } else {
          iArray = fFD->GetAbstractArray(oArray->GetName());
          for(int i = 0; i < n; i++)
            oArray->InsertTuples(q + i, 1, 0, iArray);
        }
      }

      q += n;
    }
  }

  // Edges
  {
    for(int t = 0, q = 0; t < nSteps; t++) {
      auto f0 = static_cast<vtkPointSet *>(features->GetBlock(t));
      auto f1 = static_cast<vtkPointSet *>(features->GetBlock(t + 1));
      auto l0 = this->GetInputArrayToProcess(1, f0);
      auto l1 = this->GetInputArrayToProcess(1, f1);

      auto c = static_cast<vtkImageData *>(correspondences->GetBlock(t));
      int dim[3];
      c->GetDimensions(dim);
      auto matrix = this->GetInputArrayToProcess(0, c);
      auto indexLabelMap0 = c->GetFieldData()->GetArray("IndexLabelMap0");
      auto indexLabelMap1 = c->GetFieldData()->GetArray("IndexLabelMap1");

      ttkTypeMacroAI(
        matrix->GetDataType(), l0->GetDataType(),
        (generateEdges<T0, T1>(
          trackingGraph, q, trackingGraphCD, c->GetPointData(),
          ttkUtils::GetPointer<T0>(matrix), nodeIdxOffsets[t],
          nodeIdxOffsets[t + 1], dim, f0->GetNumberOfPoints(),
          f1->GetNumberOfPoints(), ttkUtils::GetPointer<T1>(l0),
          ttkUtils::GetPointer<T1>(l1), indexLabelMap0, indexLabelMap1)));
    }
  }

  this->printMsg(msg, 1, timer.getElapsedTime());

  return 1;
}

int ttkTrackingGraph::GenerateTrackingGraphFromMatrix(
  vtkPolyData *trackingGraph, vtkMultiBlockDataSet *correspondences) {

  // initialize output
  ttk::Timer timer;
  std::string msg = "Initializing Output";
  this->printMsg(msg, 0, 0, ttk::debug::LineMode::REPLACE);

  const int nSteps = correspondences->GetNumberOfBlocks();

  int nNodes, nEdges;
  std::vector<int> nodeIdxOffsets;
  if(!this->CountNodesAndEdges(nNodes, nEdges, nodeIdxOffsets, correspondences))
    return 0;

  ArrayMap correspondencesPD;
  collectArrays(correspondencesPD, correspondences, 0);

  ArrayMap correspondencesFD;
  collectArrays(correspondencesFD, correspondences, 2);

  // allocate memory
  auto points = vtkSmartPointer<vtkPoints>::New();
  points->SetDataTypeToFloat();
  points->SetNumberOfPoints(nNodes);

  auto trackingGraphPD = trackingGraph->GetPointData();
  auto trackingGraphCD = trackingGraph->GetCellData();
  for(auto ait : correspondencesPD) {
    auto array
      = vtkSmartPointer<vtkAbstractArray>::Take(ait.second->NewInstance());
    array->SetName(ait.second->GetName());
    array->SetNumberOfComponents(ait.second->GetNumberOfComponents());
    array->SetNumberOfTuples(nEdges);
    trackingGraphCD->AddArray(array);
  }

  vtkSmartPointer<vtkDataArray> labels;
  {
    auto someIndexLabelMap = correspondencesFD.find("IndexLabelMap0");
    labels = vtkSmartPointer<vtkDataArray>::Take(
      vtkDataArray::CreateDataArray(someIndexLabelMap->second->GetDataType()));
    labels->SetName(static_cast<vtkStringArray *>(
                      correspondencesFD.find("LabelIdentifier")->second)
                      ->GetValue(0)
                      .data());
    labels->SetNumberOfTuples(nNodes);
    trackingGraphPD->AddArray(labels);
  }

  auto timeIdx = vtkSmartPointer<vtkIntArray>::New();
  {
    timeIdx->SetName("TimeIdx");
    timeIdx->SetNumberOfTuples(nNodes);
    trackingGraphPD->AddArray(timeIdx);
  }

  trackingGraph->SetPoints(points);
  trackingGraph->AllocateExact(0, 0, nEdges, 2, 0, 0, 0, 0);

  this->printMsg(msg, 1, timer.getElapsedTime());

  // Generate Tracking Graph
  timer.reStart();
  msg = "Generating Tracking Graph";
  this->printMsg(msg, 0, 0, ttk::debug::LineMode::REPLACE);

  auto pointCoordsData = ttkUtils::GetPointer<float>(points->GetData());
  auto timeIdxData = ttkUtils::GetPointer<int>(timeIdx);

  for(int t = 0, nodeIdx = 0, nodeIdx3 = 0, edgeIdx = 0; t < nSteps; t++) {
    auto c = static_cast<vtkImageData *>(correspondences->GetBlock(t));
    int dim[3];
    c->GetDimensions(dim);

    auto cFD = c->GetFieldData();
    auto indexLabelMap0 = cFD->GetArray("IndexLabelMap0");
    auto indexLabelMap1 = cFD->GetArray("IndexLabelMap1");

    // nodes
    if(t == 0) {
      labels->InsertTuples(0, dim[0], 0, indexLabelMap0);
      for(int i = 0; i < dim[0]; i++) {
        pointCoordsData[nodeIdx3++] = 0;
        pointCoordsData[nodeIdx3++] = i;
        pointCoordsData[nodeIdx3++] = 0;

        timeIdxData[nodeIdx++] = 0;
      }
    }

    labels->InsertTuples(nodeIdx, dim[1], 0, indexLabelMap1);
    for(int i = 0; i < dim[1]; i++) {
      pointCoordsData[nodeIdx3++] = t + 1;
      pointCoordsData[nodeIdx3++] = i;
      pointCoordsData[nodeIdx3++] = 0;

      timeIdxData[nodeIdx++] = t + 1;
    }

    // edges
    if(dim[0] < 1 || dim[1] < 1)
      continue;

    auto cMatrix = this->GetInputArrayToProcess(0, c);
    if(!cMatrix)
      return !this->printErr("Unable to retrieve correspondence matrix.");

    switch(cMatrix->GetDataType()) {
      vtkTemplateMacro((generateEdges<VTK_TT, int>(
        trackingGraph, edgeIdx, trackingGraphCD, c->GetPointData(),
        ttkUtils::GetPointer<VTK_TT>(cMatrix), nodeIdxOffsets[t + 0],
        nodeIdxOffsets[t + 1], dim)));
    }
  }

  this->printMsg(msg, 1, timer.getElapsedTime());

  return 1;
}

int ttkTrackingGraph::RequestData(vtkInformation *request,
                                  vtkInformationVector **inputVector,
                                  vtkInformationVector *outputVector) {
  // get input / output
  auto correspondences = vtkMultiBlockDataSet::GetData(inputVector[0]);
  auto features = vtkMultiBlockDataSet::GetData(inputVector[1]);
  auto trackingGraph = vtkPolyData::GetData(outputVector);

  if(!this->Validate(correspondences, features))
    return 0;

  if(features) {
    if(!this->GenerateTrackingGraphFromFeatures(
         trackingGraph, correspondences, features))
      return 0;
  } else {
    if(!this->GenerateTrackingGraphFromMatrix(trackingGraph, correspondences))
      return 0;
  }

  return 1;
}