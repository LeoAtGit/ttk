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

#include <ttkMacros.h>
#include <ttkUtils.h>

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
    if(port==1)
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
int countEdges(
  int &nEdges,
  const CT *correspondences,
  const int* dim,
  const int nLabels0 = 0,
  const int nLables1 = 0,
  const IT* labels0 = nullptr,
  const IT* labels1 = nullptr,
  const int* indexLabelMap0 = nullptr,
  const int* indexLabelMap1 = nullptr
) {

  if(nLabels0>0 && nLables1>0){
    std::unordered_map<int,int> labelIndexMap0;
    std::unordered_map<int,int> labelIndexMap1;

    for(int i=0; i<dim[0]; i++)
      labelIndexMap0.emplace(std::make_pair(indexLabelMap0[i],i));
    for(int i=0; i<dim[1]; i++)
      labelIndexMap1.emplace(std::make_pair(indexLabelMap1[i],i));

    for(int i=0; i<nLabels0; i++){
      for(int j=0; j<nLables1; j++){
        const auto& iLabel = labels0[i];
        const auto& jLabel = labels1[j];
        const auto& iIt = labelIndexMap0.find(iLabel);
        const auto& jIt = labelIndexMap1.find(jLabel);
        if(iIt==labelIndexMap0.end() || jIt==labelIndexMap1.end())
          continue;

        const auto& iIdx = iIt->second;
        const auto& jIdx = jIt->second;
        const int index = jIdx * dim[0] + iIdx;
        if(correspondences[index] > 0)
          nEdges++;
      }
    }
  } else {
    for(int i=0; i<dim[0]; i++)
      for(int j=0; j<dim[1]; j++){
        const auto& iIdx = i;
        const auto& jIdx = j;
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
                  vtkFieldData *outputCD,
                  vtkFieldData *inputPD,
                  const CT *correspondences,
                  const int &offset0,
                  const int &offset1,
                  const int *dim,
                  const int nLabels0=0,
                  const int nLables1=0,
                  const IT* labels0=nullptr,
                  const IT* labels1=nullptr,
                  const int* indexLabelMap0=nullptr,
                  const int* indexLabelMap1=nullptr
  ) {

  std::vector<std::pair<vtkAbstractArray *, vtkAbstractArray *>> arrayMap;
  for(int a = 0; a < outputCD->GetNumberOfArrays(); a++) {
    auto oArray = outputCD->GetArray(a);
    if(oArray)
      arrayMap.push_back({oArray, inputPD->GetArray(oArray->GetName())});
  }

  if(nLabels0>0 && nLables1>0){
    std::unordered_map<int,int> labelIndexMap0;
    std::unordered_map<int,int> labelIndexMap1;

    for(int i=0; i<dim[0]; i++)
      labelIndexMap0.emplace(std::make_pair(indexLabelMap0[i],i));
    for(int i=0; i<dim[1]; i++)
      labelIndexMap1.emplace(std::make_pair(indexLabelMap1[i],i));

    for(int i=0; i<nLabels0; i++){
      for(int j=0; j<nLables1; j++){
        const auto& iLabel = labels0[i];
        const auto& jLabel = labels1[j];
        const auto& iIt = labelIndexMap0.find(iLabel);
        const auto& jIt = labelIndexMap1.find(jLabel);
        if(iIt==labelIndexMap0.end() || jIt==labelIndexMap1.end())
          continue;

        const auto& iIdx = iIt->second;
        const auto& jIdx = jIt->second;
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
    for(int i=0; i<dim[0]; i++)
      for(int j=0; j<dim[1]; j++){
        const auto& iIdx = i;
        const auto& jIdx = j;
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

  return 1;
}

int ttkTrackingGraph::GeneratePlanarTrackingGraph(
  vtkPolyData* output,
  vtkMultiBlockDataSet* correspondences
) {
  ttk::Timer timer;
  std::string msg = "Initializing Output";
  this->printMsg(msg, 0, 0, ttk::debug::LineMode::REPLACE);

  const int nSteps = correspondences->GetNumberOfBlocks();

  int nNodes = 0;
  int nEdges = 0;
  std::vector<int> pointIndexOffset(1, 0);
  typedef std::unordered_map<std::string, vtkAbstractArray *> ArrayMap;
  ArrayMap iPointData;

  for(int t=0; t<nSteps; t++) {
    auto correspondencesPC = vtkImageData::SafeDownCast(correspondences->GetBlock(t));
    if(!correspondencesPC)
      return !this->printErr("Correspondences must be a list of vtkImageData objects.");

    int dim[3];
    correspondencesPC->GetDimensions(dim);

    if(t==0)
      nNodes+=dim[0];
    nNodes+=dim[1];

    pointIndexOffset.push_back(dim[0] + pointIndexOffset[t]);

    if(dim[0]<1 || dim[1]<1)
      continue;

    auto pd = correspondencesPC->GetPointData();
    for(int a=0; a<pd->GetNumberOfArrays(); a++){
      auto array = pd->GetArray(a);
      if(array)
        iPointData.insert({array->GetName(),array});
    }

    auto matrix = this->GetInputArrayToProcess(0, correspondencesPC);
    if(!matrix)
      return !this->printErr("Unable to retrieve correspondence matrix.");

    switch(matrix->GetDataType()) {
      vtkTemplateMacro(
        (countEdges<VTK_TT,int>(
          nEdges,
          ttkUtils::GetPointer<VTK_TT>(matrix),
          dim
        ))
      );
    }
  }

  auto outputCD = output->GetCellData();
  for(auto ait : iPointData) {
    auto array
      = vtkSmartPointer<vtkAbstractArray>::Take(ait.second->NewInstance());
    array->SetName(ait.second->GetName());
    array->SetNumberOfComponents(ait.second->GetNumberOfComponents());
    array->SetNumberOfTuples(nEdges);
    outputCD->AddArray(array);
  }

  auto points = vtkSmartPointer<vtkPoints>::New();
  points->SetDataTypeToFloat();
  points->SetNumberOfPoints(nNodes);
  auto pointCoordsData = ttkUtils::GetPointer<float>(points->GetData());

  auto sequenceId = vtkSmartPointer<vtkIntArray>::New();
  sequenceId->SetName("SequenceId");
  sequenceId->SetNumberOfTuples(nNodes);
  auto sequenceIdData = ttkUtils::GetPointer<int>(sequenceId);
  output->GetPointData()->AddArray(sequenceId);

  output->SetPoints(points);
  output->AllocateExact(0, 0, nEdges, 2, 0, 0, 0, 0);

  this->printMsg(msg, 1, timer.getElapsedTime());
  timer.reStart();
  msg = "Generating Planar Tracking Graph";
  this->printMsg(msg, 0, 0, ttk::debug::LineMode::REPLACE);

  int iNode = 0;
  int iNodeC = 0;
  int iEdge = 0;
  for(int t=0; t<nSteps; t++) {
    auto correspondencesPC = vtkImageData::SafeDownCast(correspondences->GetBlock(t));
    if(!correspondencesPC)
      return !this->printErr("Correspondences must be a list of vtkImageData objects.");

    int dim[3];
    correspondencesPC->GetDimensions(dim);

    if(t==0){
      for(int i=0; i<dim[0]; i++){
        pointCoordsData[iNodeC++] = 0;
        pointCoordsData[iNodeC++] = i;
        pointCoordsData[iNodeC++] = 0;

        sequenceIdData[iNode++] = 0;
      }
    }

    for(int i=0; i<dim[1]; i++){
      pointCoordsData[iNodeC++] = t+1;
      pointCoordsData[iNodeC++] = i;
      pointCoordsData[iNodeC++] = 0;

      sequenceIdData[iNode++] = t+1;
    }

    if(dim[0]<1 || dim[1]<1)
      continue;

    auto matrix = this->GetInputArrayToProcess(0, correspondencesPC);
    if(!matrix)
      return !this->printErr("Unable to retrieve correspondence matrix.");

    switch(matrix->GetDataType()) {
      vtkTemplateMacro(
        (generateEdges<VTK_TT,int>(
          output,
          iEdge,
          outputCD,
          correspondencesPC->GetPointData(),
          ttkUtils::GetPointer<VTK_TT>(matrix),
          pointIndexOffset[t + 0],
          pointIndexOffset[t + 1],
          dim
        ))
      );
    }
  }

  this->printMsg(msg, 1, timer.getElapsedTime());

  return 1;
}

int ttkTrackingGraph::GenerateSpatialTrackingGraph(
  vtkPolyData* output,
  vtkMultiBlockDataSet* correspondences,
  vtkMultiBlockDataSet* features
) {
  ttk::Timer timer;
  std::string msg = "Initializing Output";
  this->printMsg(msg, 0, 0, ttk::debug::LineMode::REPLACE);

  const int nSteps = features->GetNumberOfBlocks();
  if(nSteps != (int)correspondences->GetNumberOfBlocks() + 1)
    return !this->printErr("Incompatible vtkMultiBlockDataSet structures for "
                          "features and correspondences.");

  auto outputPD = output->GetPointData();
  auto outputCD = output->GetCellData();

  typedef std::unordered_map<std::string, vtkAbstractArray *> ArrayMap;
  typedef std::vector<std::tuple<vtkFieldData *, ArrayMap *,int>> ArrayMapSet;

  int nNodes = 0;
  int nEdges = 0;
  std::vector<int> pointIndexOffset(1, 0);
  {
    ArrayMap iPointDataMap;
    ArrayMap iCellDataMap;
    ArrayMap iFieldDataMap;

    for(int t = 0; t < nSteps; t++) {
      auto featuresC = vtkPointSet::SafeDownCast(features->GetBlock(t));
      if(!featuresC)
        return !this->printErr("Features must be a list of vtkPointSets.");

      const int nFeaturesC = featuresC->GetNumberOfPoints();
      pointIndexOffset.push_back(
        t < 1 ? nFeaturesC
              : pointIndexOffset[pointIndexOffset.size() - 1] + nFeaturesC);
      nNodes += nFeaturesC;

      ArrayMapSet iSet;
      for(auto it : ArrayMapSet({
        {featuresC->GetPointData(), &iPointDataMap,0},
        {featuresC->GetFieldData(), &iFieldDataMap,0}
      })){
        for(int a = 0; a < std::get<0>(it)->GetNumberOfArrays(); a++) {
          auto array = std::get<0>(it)->GetArray(a);
          if(array)
            std::get<1>(it)->insert(std::make_pair(array->GetName(), array));
        }
      }

      if(t==0)
        continue;

      auto correspondencesPC = vtkImageData::SafeDownCast(correspondences->GetBlock(t-1));
      if(!correspondencesPC)
        return !this->printErr("Correspondences must be a list of vtkImageData objects.");

      int dim[3];
      correspondencesPC->GetDimensions(dim);
      if(dim[0] < 1 || dim[1] < 1)
        continue;

      auto matrix = this->GetInputArrayToProcess(0, correspondencesPC);
      if(!matrix)
        return !this->printErr("Unable to retrieve correspondence matrix.");

      auto featuresP = vtkPointSet::SafeDownCast(features->GetBlock(t-1));
      const int nFeaturesP = featuresP->GetNumberOfPoints();
      if(nFeaturesP<1)
        continue;

      auto labelsP = this->GetInputArrayToProcess(1, featuresP);
      auto labelsC = this->GetInputArrayToProcess(1, featuresC);
      if(!labelsP || !labelsC)
        return !this->printErr("Label lookup requires label arrays.");

      if(labelsP->GetDataType()!=labelsC->GetDataType())
        return !this->printErr("Labels must be of same data type.");

      auto cFD = correspondencesPC->GetFieldData();
      auto indexLabelMap0 = cFD->GetArray("IndexLabelMap0");
      auto indexLabelMap1 = cFD->GetArray("IndexLabelMap1");
      if(!indexLabelMap0 || !indexLabelMap1)
        return !this->printErr("Label lookup requires IndexLabelMaps for Correspondence Matrices.");

      switch(vtkTemplate2PackMacro(labelsC->GetDataType(), matrix->GetDataType())) {
        ttkTemplate2IdMacro(
          countEdges(
            nEdges,
            ttkUtils::GetPointer<VTK_T2>(matrix),
            dim,
            nFeaturesP,
            nFeaturesC,
            ttkUtils::GetPointer<VTK_T1>(labelsP),
            ttkUtils::GetPointer<VTK_T1>(labelsC),
            ttkUtils::GetPointer<int>(indexLabelMap0),
            ttkUtils::GetPointer<int>(indexLabelMap1)
          )
        );
      }

      auto data = correspondencesPC->GetPointData();
      for(int a = 0; a < data->GetNumberOfArrays(); a++) {
        auto iArray = data->GetArray(a);
        if(iArray)
          iCellDataMap.insert(std::make_pair(iArray->GetName(), iArray));
      }
    }

    ArrayMapSet oSet{
      {outputPD, &iPointDataMap,nNodes},
      {outputPD, &iFieldDataMap,nNodes},
      {outputCD, &iCellDataMap,nEdges}
    };
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

  auto points = vtkSmartPointer<vtkPoints>::New();
  points->SetDataTypeToFloat();
  points->SetNumberOfPoints(nNodes);
  auto pointCoords = points->GetData();

  output->SetPoints(points);
  output->AllocateExact(0, 0, nEdges, 2, 0, 0, 0, 0);

  this->printMsg(msg, 1, timer.getElapsedTime());
  timer.reStart();
  msg = "Generating Spatial Tracking Graph";
  this->printMsg(msg, 0, 0, ttk::debug::LineMode::REPLACE);

  const int nPointArrays = outputPD->GetNumberOfArrays();
  int iEdge = 0;
  for(int t = 0, q = 0; t < nSteps; t++) {
    auto featuresC = vtkPointSet::SafeDownCast(features->GetBlock(t));
    const int nFeaturesC = featuresC->GetNumberOfPoints();
    if(nFeaturesC < 1)
      continue;

    auto pointData = featuresC->GetPointData();
    auto fieldData = featuresC->GetFieldData();

    pointCoords->InsertTuples(q, nFeaturesC, 0, featuresC->GetPoints()->GetData());

    for(int a = 0; a < nPointArrays; a++) {
      auto oArray = outputPD->GetAbstractArray(a);
      auto iArray = pointData->GetAbstractArray(oArray->GetName());
      if(iArray) {
        oArray->InsertTuples(q, nFeaturesC, 0, iArray);
      } else {
        iArray = fieldData->GetAbstractArray(oArray->GetName());
        for(int i = 0; i < nFeaturesC; i++)
          oArray->InsertTuples(q + i, 1, 0, iArray);
      }
    }

    q += nFeaturesC;

    if(t==0)
      continue;

    auto correspondencesPC = vtkImageData::SafeDownCast(correspondences->GetBlock(t - 1));
    int dim[3];
    correspondencesPC->GetDimensions(dim);

    auto featuresP = vtkPointSet::SafeDownCast(features->GetBlock(t-1));
    const int nFeaturesP = featuresP->GetNumberOfPoints();
    if(nFeaturesP<1 || nFeaturesC<1)
      continue;

    auto matrix = this->GetInputArrayToProcess(0, correspondencesPC);
    auto labelsP = this->GetInputArrayToProcess(1, featuresP);
    auto labelsC = this->GetInputArrayToProcess(1, featuresC);
    auto cFD = correspondencesPC->GetFieldData();
    auto indexLabelMap0 = cFD->GetArray("IndexLabelMap0");
    auto indexLabelMap1 = cFD->GetArray("IndexLabelMap1");

    switch(vtkTemplate2PackMacro(labelsC->GetDataType(), matrix->GetDataType())) {
      ttkTemplate2IdMacro(
        (generateEdges<VTK_T2,VTK_T1>(
          output,
          iEdge,
          outputCD,
          correspondencesPC->GetPointData(),
          ttkUtils::GetPointer<VTK_T2>(matrix),
          pointIndexOffset[t - 1],
          pointIndexOffset[t + 0],
          dim,
          nFeaturesP,
          nFeaturesC,
          ttkUtils::GetPointer<VTK_T1>(labelsP),
          ttkUtils::GetPointer<VTK_T1>(labelsC),
          ttkUtils::GetPointer<int>(indexLabelMap0),
          ttkUtils::GetPointer<int>(indexLabelMap1)
        ))
      );
    }
  }

  this->printMsg(msg, 1, timer.getElapsedTime());
  return 1;
}

int ttkTrackingGraph::RequestData(vtkInformation *request,
                                  vtkInformationVector **inputVector,
                                  vtkInformationVector *outputVector) {

  auto correspondences = vtkMultiBlockDataSet::GetData(inputVector[0]);
  auto features = vtkMultiBlockDataSet::GetData(inputVector[1]);
  auto output = vtkPolyData::GetData(outputVector);

  if(correspondences && features){
    if(!this->GenerateSpatialTrackingGraph(output,correspondences,features))
      return 0;
  } else if (correspondences) {
    if(!this->GeneratePlanarTrackingGraph(output,correspondences))
      return 0;
  } else {
    return !this->printErr("Unable to retrieve input data objects.");
  }

  return 1;
}
