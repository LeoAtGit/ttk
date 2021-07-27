#include <ttkBranchDecomposition.h>

#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include <ttkMacros.h>
#include <ttkUtils.h>
#include <TrackingGraph.h>

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkBranchDecomposition);

ttkBranchDecomposition::ttkBranchDecomposition() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkBranchDecomposition::~ttkBranchDecomposition() {
}

int ttkBranchDecomposition::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
    return 1;
  }
  return 0;
}

int ttkBranchDecomposition::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

template<typename IT, typename DT>
int computeBranchDecomposition(
  int* branchId,
  const IT* time,
  const DT* attribute,
  ttk::TrackingGraph& trackingGraph
){
  const int nNodes = trackingGraph.prevNodes.size();

  // sort all nodes by time in ascending order
  std::vector<int> nodesSortedByTime(nNodes);
  {
    for(int i=0; i<nNodes; i++)
      nodesSortedByTime[i]=i;
    std::sort(nodesSortedByTime.begin(),nodesSortedByTime.end(), [=](int a, int b){return time[a]<time[b];});
  }

  // sort all prev and next nodes by attribute in descending order
  {
    const auto compareAttribute = [=](int a, int b){return attribute[a]>attribute[b];};
    for(int i=0; i<nNodes; i++){
      const auto& nextNodes = trackingGraph.nextNodes[i];
      const auto& prevNodes = trackingGraph.prevNodes[i];
      if(nextNodes.size()>1)
        std::sort(nextNodes.begin(),nextNodes.end(), compareAttribute);
      if(prevNodes.size()>1)
        std::sort(prevNodes.begin(),prevNodes.end(), compareAttribute);
    }
  }

  std::vector<int> maxPrevNode(nNodes);
  std::vector<int> maxNextNode(nNodes);
  for(int i=0; i<nNodes; i++){
    maxPrevNode[i] = trackingGraph.prevNodes[i].size()>0 ? trackingGraph.prevNodes[i][0] : -1;
    maxNextNode[i] = trackingGraph.nextNodes[i].size()>0 ? trackingGraph.nextNodes[i][0] : -1;
  }

  // int branchIdCounter = 0;
  // initialize branch id at every birth node
  for(int i=0; i<nNodes; i++){
    branchId[i] = trackingGraph.prevNodes[i].size()<1 ? i : -1;
  }

  // propagate branch id along graph
  for(int i=0; i<nNodes; i++){
    const auto& v = nodesSortedByTime[i];

    // std::cout<<v<<std::endl;
    if(branchId[v]!=-1)
      continue;

    int maxPrevNodeV = -1;
    for(const auto& u: trackingGraph.prevNodes[v]){
      if(maxNextNode[u]==v){
        maxPrevNodeV = u;
        break;
      }
    }
    // for(const auto& u: trackingGraph.prevNodes[v]){
    //   std::cout<<" :"<<u<<" ("<<attribute[u]<<")"<<std::endl;
    // }
    // std::cout<<"  ->"<<maxPrevNodeV<<std::endl;

    branchId[v] = maxPrevNodeV<0 ? v : branchId[maxPrevNodeV];
  }

  return 1;
}

#define ttkTypeMacroAA(group0, group1, call) \
  switch(group1) {                                                          \
    ttkTypeMacroCase(VTK_FLOAT, float, 1, ttkTypeMacroA(group0, call));                           \
    ttkTypeMacroCase(VTK_DOUBLE, double, 1, ttkTypeMacroA(group0, call));                         \
    ttkTypeMacroCase(VTK_INT, int, 1, ttkTypeMacroA(group0, call));                               \
    ttkTypeMacroCase(VTK_UNSIGNED_INT, unsigned int, 1, ttkTypeMacroA(group0, call));             \
    ttkTypeMacroCase(VTK_CHAR, char, 1, ttkTypeMacroA(group0, call));                             \
    ttkTypeMacroCase(VTK_SIGNED_CHAR, signed char, 1, ttkTypeMacroA(group0, call));               \
    ttkTypeMacroCase(VTK_UNSIGNED_CHAR, unsigned char, 1, ttkTypeMacroA(group0, call));           \
    ttkTypeMacroCase(VTK_LONG, long, 1, ttkTypeMacroA(group0, call));                             \
    ttkTypeMacroCase(VTK_LONG_LONG, long long, 1, ttkTypeMacroA(group0, call));                   \
    ttkTypeMacroCase(VTK_UNSIGNED_LONG, unsigned long, 1, ttkTypeMacroA(group0, call));           \
    ttkTypeMacroCase(VTK_UNSIGNED_LONG_LONG, unsigned long long, 1, ttkTypeMacroA(group0, call)); \
    ttkTypeMacroCase(VTK_ID_TYPE, vtkIdType, 1, ttkTypeMacroA(group0, call));                     \
    ttkTypeMacroErrorCase(1, group1);                                       \
  }

int ttkBranchDecomposition::RequestData(vtkInformation *request,
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {

  auto input = vtkPointSet::GetData(inputVector[0]);
  if(!input)
    return 0;

  auto vtkTrackingGraph = vtkPointSet::GetData(outputVector);
  vtkTrackingGraph->ShallowCopy(input);

  auto timeArray = this->GetInputArrayToProcess(0, inputVector);
  if(!timeArray)
    return !this->printErr("Unable to retrieve time array.");

  auto attributeArray = this->GetInputArrayToProcess(1, inputVector);
  if(!attributeArray)
    return !this->printErr("Unable to retrieve attribute array.");

  if(this->GetInputArrayAssociation(0, inputVector) != 0 || this->GetInputArrayAssociation(1, inputVector) != 0)
    return !this->printErr("Input arrays need to be point data arrays.");
  if(timeArray->GetNumberOfComponents() != 1 || attributeArray->GetNumberOfComponents() != 1)
    return !this->printErr("Input arrays need to be scalar arrays.");

  const int nNodes = vtkTrackingGraph->GetNumberOfPoints();
  auto connectivityList =
    vtkTrackingGraph->IsA("vtkPolyData")
      ? static_cast<vtkPolyData*>(vtkTrackingGraph)->GetLines()->GetConnectivityArray()
      : vtkTrackingGraph->IsA("vtkUnstructuredGrid")
      ? static_cast<vtkUnstructuredGrid*>(vtkTrackingGraph)->GetCells()->GetConnectivityArray()
      : nullptr;

  if(!connectivityList)
    return !this->printErr("Unable to retrieve connectivity list.");

  const int nEdges = connectivityList->GetNumberOfValues()/2;

  ttk::TrackingGraph ttkTrackingGraph;
  ttkTrackingGraph.setDebugLevel(this->debugLevel_);
  ttkTypeMacroI(
    connectivityList->GetDataType(),
    ttkTrackingGraph.preconditionPrevAndNextNodes<T0>(
      nNodes,
      nEdges,
      ttkUtils::GetConstPointer<const T0>(connectivityList)
    )
  );

  auto branchId = vtkSmartPointer<vtkIntArray>::New();
  branchId->SetName("BranchId");
  branchId->SetNumberOfComponents(1);
  branchId->SetNumberOfTuples(nNodes);
  vtkTrackingGraph->GetPointData()->AddArray(branchId);

  ttkTypeMacroAA(
    timeArray->GetDataType(),
    attributeArray->GetDataType(),
    (computeBranchDecomposition<T0,T1>(
      ttkUtils::GetPointer<int>(branchId),
      ttkUtils::GetConstPointer<const T0>(timeArray),
      ttkUtils::GetConstPointer<const T1>(attributeArray),
      ttkTrackingGraph
    ))
  );

  // std::vector<int> nodesSortedByTime(nNodes);
  // for(int i=0; i<nNodes; i++)
  //   nodesSortedByTime[i]=i;

  // std::sort(nodesSortedByTime.begin(),nodesSortedByTime.end(), [=](int a, int b){return in});


  // return success
  return 1;
}
