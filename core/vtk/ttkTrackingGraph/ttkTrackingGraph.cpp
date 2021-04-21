#include <ttkTrackingGraph.h>

#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPolyData.h>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkCellData.h>

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
    return 1;
  }
  return 0;
}

int ttkTrackingGraph::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
    return 1;
  }
  return 0;
}

template<typename DT>
int countEdges(
  int& nEdges,
  const DT* correlations,
  const int nCorrelations
){
  for(int i=0; i<nCorrelations; i++)
    if(correlations[i]>0)
      nEdges++;
  return 1;
}

template<typename DT>
int generateEdges(
  vtkPolyData* output,
  int& iEdge,
  vtkFieldData* outputCD,
  vtkFieldData* inputPD,
  const DT* correlations,
  const int& offset0,
  const int& offset1,
  const int& nLabels0,
  const int& nLabels1
){
  std::vector<std::pair<vtkAbstractArray*,vtkAbstractArray*>> arrayMap;
  for(int a=0; a<outputCD->GetNumberOfArrays(); a++){
    auto oArray = outputCD->GetArray(a);
    if(oArray)
      arrayMap.push_back({oArray,inputPD->GetArray(oArray->GetName())});
  }

  for(int i=0; i<nLabels0; i++){
    for(int j=0; j<nLabels1; j++){
      const int index = j*nLabels0 + i;
      const auto& c = correlations[index];
      if(c<=0)
        continue;

      vtkIdType points[2]{
        offset0+i,
        offset1+j
      };
      output->InsertNextCell(VTK_LINE, 2, points);

      for(auto& it: arrayMap)
        it.first->SetTuple(iEdge, index, it.second);

      iEdge++;
    }
  }

  return 1;
}

int ttkTrackingGraph::RequestData(vtkInformation *request,
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {

  auto nodes = vtkMultiBlockDataSet::GetData(inputVector[0]);
  auto edges = vtkMultiBlockDataSet::GetData(inputVector[1]);
  if(!nodes || !edges)
    return !this->printErr("Unable to retrieve input data objects.");

  const int nSteps = nodes->GetNumberOfBlocks();
  if(nSteps!=(int)edges->GetNumberOfBlocks()+1)
    return !this->printErr("Incompatible vtkMultiBlockDataSet structures for nodes and correlations.");

  auto output = vtkPolyData::GetData(outputVector);
  auto outputPD = output->GetPointData();
  auto outputCD = output->GetCellData();

  typedef std::unordered_map<std::string,vtkAbstractArray*> ArrayMap;
  typedef std::vector<std::pair<vtkFieldData*,ArrayMap*>> ArrayMapSet;

  int nNodes = 0;
  std::vector<int> pointIndexOffset(1,0);
  {
    ArrayMap iPointDataMap;
    ArrayMap iFieldDataMap;
    for(int t=0; t<nSteps; t++){
      auto nodesT = vtkPointSet::SafeDownCast(nodes->GetBlock(t));
      if(!nodesT)
        return !this->printErr("Nodes must be a list of vtkPointSets.");

      const int nNodesT = nodesT->GetNumberOfPoints();
      pointIndexOffset.push_back(
        t<1 ? nNodesT : pointIndexOffset[pointIndexOffset.size()-1] + nNodesT
      );
      nNodes+=nNodesT;

      ArrayMapSet iSet{
        {nodesT->GetPointData(),&iPointDataMap},
        {nodesT->GetFieldData(),&iFieldDataMap}
      };
      for(auto type: iSet){
        for(int a=0; a<type.first->GetNumberOfArrays(); a++){
          auto array = type.first->GetArray(a);
          if(array)
            type.second->insert(std::make_pair(array->GetName(), array));
        }
      }
    }

    ArrayMapSet oSet{
      {outputPD,&iPointDataMap},
      {outputPD,&iFieldDataMap}
    };
    for(auto type: oSet){
      for(auto it: *(type.second)){
        auto array = vtkSmartPointer<vtkAbstractArray>::Take(
          it.second->NewInstance()
        );
        array->SetName(it.second->GetName());
        array->SetNumberOfComponents(it.second->GetNumberOfComponents());
        array->SetNumberOfTuples(nNodes);
        type.first->AddArray(array);
      }
    }
  }

  int nEdges = 0;
  {
    ArrayMap iCellData;
    for(int t=0; t<nSteps-1; t++){
      auto edgesT = vtkImageData::SafeDownCast(edges->GetBlock(t));
      if(!edgesT)
        return !this->printErr("Edges must be a list of vtkImageData objects.");

      int dim[3];
      edgesT->GetDimensions(dim);
      if(dim[0]<1 || dim[1]<1)
        continue;

      auto array = this->GetInputArrayToProcess(0, edgesT);
      if(!array)
        return !this->printErr("Unable to retrieve input array");

      switch(array->GetDataType()){
        vtkTemplateMacro(countEdges(
          nEdges,
          ttkUtils::GetPointer<VTK_TT>(array),
          array->GetNumberOfTuples()
        ));
      }

      auto data = edgesT->GetPointData();
      for(int a=0; a<data->GetNumberOfArrays(); a++){
        auto iArray = data->GetArray(a);
        if(iArray)
          iCellData.insert(std::make_pair(iArray->GetName(), iArray));
      }
    }

    for(auto it: iCellData){
      auto array = vtkSmartPointer<vtkAbstractArray>::Take(
        it.second->NewInstance()
      );
      array->SetName(it.second->GetName());
      array->SetNumberOfComponents(it.second->GetNumberOfComponents());
      array->SetNumberOfTuples(nEdges);
      outputCD->AddArray(array);
    }
  }

  auto points = vtkSmartPointer<vtkPoints>::New();
  points->SetDataTypeToFloat();
  points->SetNumberOfPoints(nNodes);
  auto pointCoords = points->GetData();

  output->SetPoints(points);
  output->AllocateExact(0,0,nEdges,2,0,0,0,0);

  const int nPointArrays = outputPD->GetNumberOfArrays();

  int iEdge=0;
  for(int t=0,q=0; t<nSteps; t++){
    auto nodesT = vtkPointSet::SafeDownCast(nodes->GetBlock(t));
    const int nPoints = nodesT->GetNumberOfPoints();
    if(nPoints<1)
      continue;

    auto pointData = nodesT->GetPointData();
    auto fieldData = nodesT->GetFieldData();

    pointCoords->InsertTuples(q, nPoints, 0, nodesT->GetPoints()->GetData());

    for(int a=0; a<nPointArrays; a++){
      auto oArray = outputPD->GetAbstractArray(a);
      auto iArray = pointData->GetAbstractArray(oArray->GetName());
      if(iArray){
        oArray->InsertTuples(q, nPoints, 0, iArray);
      } else {
        iArray = fieldData->GetAbstractArray(oArray->GetName());
        for(int i=0; i<nPoints; i++)
          oArray->InsertTuples(q+i, 1, 0, iArray);
      }
    }

    q+=nPoints;

    if(t<1)
      continue;

    auto edgesT = vtkImageData::SafeDownCast(edges->GetBlock(t-1));
    int dim[3];
    edgesT->GetDimensions(dim);

    auto array = this->GetInputArrayToProcess(0, edgesT);
    switch(array->GetDataType()){
      vtkTemplateMacro(generateEdges(
        output,
        iEdge,
        outputCD,
        edgesT->GetPointData(),
        ttkUtils::GetPointer<VTK_TT>(array),
        pointIndexOffset[t-1],
        pointIndexOffset[t+0],
        dim[0], dim[1]
      ));
    }
  }

  return 1;
}
