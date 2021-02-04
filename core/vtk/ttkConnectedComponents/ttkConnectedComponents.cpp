#include <ttkConnectedComponents.h>

#include <vtkDataArray.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkConnectedComponents);

ttkConnectedComponents::ttkConnectedComponents() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(2);
}

ttkConnectedComponents::~ttkConnectedComponents() {
}

int ttkConnectedComponents::FillInputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

int ttkConnectedComponents::FillOutputPortInformation(int port,
                                                      vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
    return 1;
  } else if(port == 1) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkConnectedComponents::RequestData(vtkInformation *request,
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
  auto outputArray
    = vtkSmartPointer<vtkDataArray>::Take(inputArray->NewInstance());
  outputArray->DeepCopy(inputArray);

  // Compute Connected Components
  std::vector<ttk::ConnectedComponents::Component> components;
  {

    this->preconditionTriangulation(triangulation);

    int status = 0;
    ttkVtkTemplateMacro(
      inputArray->GetDataType(), triangulation->getType(),
      (status = this->computeConnectedComponents<VTK_TT, TTK_TT>(
         components,
         static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(outputArray)),
         static_cast<const VTK_TT *>(ttkUtils::GetVoidPointer(inputArray)),
         static_cast<const TTK_TT *>(triangulation->getData()))));

    if(status != 1)
      return 0;
  }

  // Components Output
  {
    const auto nComponents = components.size();
    auto outputComponents = vtkPolyData::GetData(outputVector, 0);

    // points
    {
      auto sizeArray = vtkSmartPointer<vtkIntArray>::New();
      sizeArray->SetName("Size");
      sizeArray->SetNumberOfTuples(nComponents);
      auto sizeArrayData
        = static_cast<int *>(ttkUtils::GetVoidPointer(sizeArray));

      auto points = vtkSmartPointer<vtkPoints>::New();
      points->SetNumberOfPoints(nComponents);
      auto pointsData = static_cast<float *>(ttkUtils::GetVoidPointer(points));
      for(size_t i = 0, j = 0; i < nComponents; i++) {
        const auto &c = components[i];
        pointsData[j++] = c.center[0];
        pointsData[j++] = c.center[1];
        pointsData[j++] = c.center[2];

        sizeArrayData[i] = c.size;
      }

      outputComponents->SetPoints(points);
      outputComponents->GetPointData()->AddArray(sizeArray);
    }

    // cells
    {
      auto connectivityArray = vtkSmartPointer<vtkIdTypeArray>::New();
      connectivityArray->SetNumberOfTuples(nComponents);
      auto connectivityArrayData
        = static_cast<vtkIdType *>(ttkUtils::GetVoidPointer(connectivityArray));
      for(size_t i = 0; i < nComponents; i++)
        connectivityArrayData[i] = i;

      auto offsetArray = vtkSmartPointer<vtkIdTypeArray>::New();
      offsetArray->SetNumberOfTuples(nComponents + 1);
      auto offsetArrayData
        = static_cast<vtkIdType *>(ttkUtils::GetVoidPointer(offsetArray));
      for(size_t i = 0; i <= nComponents; i++)
        offsetArrayData[i] = i;

      auto cellArray = vtkSmartPointer<vtkCellArray>::New();
      cellArray->SetData(offsetArray, connectivityArray);

      outputComponents->SetVerts(cellArray);
    }

    // Copy Field Data
    outputComponents->GetFieldData()->ShallowCopy(inputDataSet->GetFieldData());
  }

  // Segmentation Output
  {
    auto outputDataSet = vtkDataSet::GetData(outputVector, 1);
    outputDataSet->ShallowCopy(inputDataSet);
    outputDataSet->GetPointData()->AddArray(outputArray);
  }

  // return success
  return 1;
}