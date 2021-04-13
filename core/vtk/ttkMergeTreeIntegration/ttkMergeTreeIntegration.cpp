#include <ttkMergeTreeIntegration.h>

#include <vtkInformation.h>
#include <vtkDataSet.h>

#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkCellData.h>

#include <vtkIntArray.h>
#include <vtkIdTypeArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkSignedCharArray.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkMergeTreeIntegration);

ttkMergeTreeIntegration::ttkMergeTreeIntegration(){
    this->SetNumberOfInputPorts(2);
    this->SetNumberOfOutputPorts(1);
}

ttkMergeTreeIntegration::~ttkMergeTreeIntegration(){}

int ttkMergeTreeIntegration::FillInputPortInformation(int port, vtkInformation* info) {
    switch(port){
        case 0:
            info->Set( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid" );
            return 1;
        case 1:
            info->Set( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet" );
            return 1;
        default:
            return 0;
    }
}

int ttkMergeTreeIntegration::FillOutputPortInformation(int port, vtkInformation* info) {
    switch(port){
        case 0:
            info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
            return 1;
        default:
            return 0;
    }
}

int ttkMergeTreeIntegration::RequestData(
    vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector
){
    ttk::Timer timer;
    this->printMsg("Integrating Arrays",0,0,1,ttk::debug::LineMode::REPLACE);

    auto mergeTree = vtkDataSet::GetData(outputVector);
    mergeTree->ShallowCopy( vtkDataSet::GetData(inputVector[0]) );
    const size_t nMerteTreeVertices = mergeTree->GetNumberOfPoints();

    auto segmentation = vtkDataSet::GetData(inputVector[1]);

    auto array = this->GetInputArrayToProcess(0, segmentation);
    if(!array)
      return !this->printErr("Unable to retrieve input array.");
    const size_t nSegVertices = array->GetNumberOfTuples();
    const size_t nComponents = array->GetNumberOfComponents();
    auto arrayData = ttkUtils::GetPointer<float>(array);

    auto nodeIdArray = this->GetInputArrayToProcess(1, segmentation);
    if(!nodeIdArray)
      return !this->printErr("Unable to retrieve node id array.");
    auto nodeIdArrayData = ttkUtils::GetPointer<int>(nodeIdArray);

    auto i0Array = vtkSmartPointer<vtkFloatArray>::New();
    i0Array->SetName((std::string(array->GetName())+"_I0").data());
    i0Array->SetNumberOfComponents(nComponents);
    i0Array->SetNumberOfTuples(nMerteTreeVertices);
    mergeTree->GetPointData()->AddArray(i0Array);
    auto i0ArrayData = ttkUtils::GetPointer<float>(i0Array);

    // clear
    for(size_t i=0, j=nMerteTreeVertices*nComponents; i<j; i++)
      i0ArrayData[i] = 0;

    for(size_t v=0; v<nSegVertices; v++){
      const auto& nodeId = nodeIdArrayData[v];

      for(size_t c=0; c<nComponents; c++){
        i0ArrayData[nodeId*nComponents+c] += arrayData[v*nComponents+c];
      }
    }

    auto nextIdArray = this->GetInputArrayToProcess(2, mergeTree);
    if(!nextIdArray)
      return !this->printErr("Unable to retrieve next id array.");
    auto nextIdArrayData = ttkUtils::GetPointer<int>(nextIdArray);

    auto orderArray = this->GetOrderArray(mergeTree, 3);
    if(!orderArray)
      return 0;
    const auto orderArrayData = ttkUtils::GetPointer<const ttk::SimplexId>(orderArray);

    std::vector<int> sortedMergeTreeVertices(nMerteTreeVertices);
    for(size_t i=0;i<nMerteTreeVertices; i++)
      sortedMergeTreeVertices[i] = i;

    std::sort(
      sortedMergeTreeVertices.begin(),
      sortedMergeTreeVertices.end(),
      [=](const int& a, const int& b) -> bool {
        return orderArrayData[a] > orderArrayData[b];
      }
    );

    auto i1Array = vtkSmartPointer<vtkFloatArray>::New();
    i1Array->DeepCopy(i0Array);
    i1Array->SetName((std::string(array->GetName())+"_I1").data());
    mergeTree->GetPointData()->AddArray(i1Array);
    auto i1ArrayData = ttkUtils::GetPointer<float>(i1Array);

    for(size_t v=0; v<nMerteTreeVertices; v++){
      auto curr = sortedMergeTreeVertices[v];
      auto next = nextIdArrayData[curr];
      if(next<0)
        continue;

      for(size_t c=0; c<nComponents; c++){
        i1ArrayData[next*nComponents+c] += i1ArrayData[curr*nComponents+c];
      }
    }

    this->printMsg("Integrating Arrays",1,timer.getElapsedTime(),this->threadNumber_);

    return 1;
}