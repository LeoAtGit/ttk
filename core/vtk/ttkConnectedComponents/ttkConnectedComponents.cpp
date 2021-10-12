#include <ttkConnectedComponents.h>

#include <vtkDataArray.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include <vtkFloatArray.h>
#include <vtkIntArray.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkConnectedComponents);

ttkConnectedComponents::ttkConnectedComponents() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(2);
}

ttkConnectedComponents::~ttkConnectedComponents() {
}

int ttkConnectedComponents::FillInputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  } else if(port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), 1);
    return 1;
  }
  return 0;
}

int ttkConnectedComponents::FillOutputPortInformation(int port,
                                                      vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  } else if(port == 1) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
    return 1;
  }
  return 0;
}

int ttkConnectedComponents::RequestData(vtkInformation *,
                                        vtkInformationVector **inputVector,
                                        vtkInformationVector *outputVector) {

  // Fetch Input Data
  auto inputDataSet = vtkDataSet::GetData(inputVector[0]);
  if(!inputDataSet)
    return 0;

  auto inputArray = this->GetInputArrayToProcess(0, inputVector);
  if(!inputArray) {
    this->printErr("Unable to retrieve input array.");
    return 0;
  }

  if(this->GetInputArrayAssociation(0, inputVector) != 0) {
    this->printErr("Input array needs to be a point data array.");
    return 0;
  }
  if(inputArray->GetNumberOfComponents() != 1) {
    this->printErr("Input array needs to be a scalar array.");
    return 0;
  }

  auto triangulation = ttkAlgorithm::GetTriangulation(inputDataSet);
  if(!triangulation)
    return 0;

  // Allocate Output Label Data
  auto outputArray = vtkSmartPointer<vtkIntArray>::New();
  outputArray->SetName("ComponentId");
  outputArray->SetNumberOfTuples(inputArray->GetNumberOfTuples());

  // Compute Connected Components
  std::vector<ttk::ConnectedComponents::Component> components;
  {

    this->preconditionTriangulation(triangulation);

    int status = 0;
    ttkVtkTemplateMacro(
      inputArray->GetDataType(), triangulation->getType(),
      (status = this->computeConnectedComponents<VTK_TT, TTK_TT>(
         components, ttkUtils::GetPointer<int>(outputArray),
         ttkUtils::GetPointer<const VTK_TT>(inputArray),
         static_cast<const TTK_TT *>(triangulation->getData()),
         this->UseSeedIdAsComponentId)));

    if(status != 1)
      return 0;
  }

  // Segmentation Output
  {
    auto outputDataSet = vtkDataSet::GetData(outputVector, 0);
    outputDataSet->ShallowCopy(inputDataSet);
    outputDataSet->GetPointData()->AddArray(outputArray);
  }

  // Components Output
  {
    const int nComponents = components.size();
    auto outputComponents = vtkPolyData::GetData(outputVector, 1);

    // points
    {
      auto sizeArray = vtkSmartPointer<vtkFloatArray>::New();
      sizeArray->SetName("Size");
      sizeArray->SetNumberOfTuples(nComponents);
      auto sizeArrayData = ttkUtils::GetPointer<float>(sizeArray);

      auto idArray = vtkSmartPointer<vtkIntArray>::New();
      idArray->SetName("ComponentId");
      idArray->SetNumberOfTuples(nComponents);
      auto idArrayData = ttkUtils::GetPointer<int>(idArray);

      auto points = vtkSmartPointer<vtkPoints>::New();
      points->SetDataTypeToFloat();
      points->SetNumberOfPoints(nComponents);
      auto pointsData = ttkUtils::GetPointer<float>(points->GetData());
      for(int i = 0, j = 0; i < nComponents; i++) {
        const auto &c = components[i];
        pointsData[j++] = c.center[0];
        pointsData[j++] = c.center[1];
        pointsData[j++] = c.center[2];

        sizeArrayData[i] = c.size;
        idArrayData[i] = this->UseSeedIdAsComponentId ? c.seed : i;
      }

      outputComponents->SetPoints(points);
      auto pd = outputComponents->GetPointData();
      pd->AddArray(sizeArray);
      pd->AddArray(idArray);
    }

    // cells
    {
      auto connectivityArray = vtkSmartPointer<vtkIntArray>::New();
      connectivityArray->SetNumberOfTuples(nComponents);
      auto connectivityArrayData = ttkUtils::GetPointer<int>(connectivityArray);
      for(int i = 0; i < nComponents; i++)
        connectivityArrayData[i] = i;

      auto offsetArray = vtkSmartPointer<vtkIntArray>::New();
      offsetArray->SetNumberOfTuples(nComponents + 1);
      auto offsetArrayData = ttkUtils::GetPointer<int>(offsetArray);
      for(int i = 0; i <= nComponents; i++)
        offsetArrayData[i] = i;

      auto cellArray = vtkSmartPointer<vtkCellArray>::New();
      cellArray->SetData(offsetArray, connectivityArray);

      outputComponents->SetVerts(cellArray);
    }

    // Copy Field Data
    outputComponents->GetFieldData()->ShallowCopy(inputDataSet->GetFieldData());
  }

  // return success
  return 1;
}
