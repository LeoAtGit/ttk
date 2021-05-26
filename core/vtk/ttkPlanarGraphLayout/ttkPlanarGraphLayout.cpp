#include <ttkPlanarGraphLayout.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

#include <vtkAbstractArray.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

vtkStandardNewMacro(ttkPlanarGraphLayout);

ttkPlanarGraphLayout::ttkPlanarGraphLayout() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}
ttkPlanarGraphLayout::~ttkPlanarGraphLayout() {
}

int ttkPlanarGraphLayout::FillInputPortInformation(int port,
                                                   vtkInformation *info) {
  if(port == 0)
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
  else
    return 0;
  return 1;
}

int ttkPlanarGraphLayout::FillOutputPortInformation(int port,
                                                    vtkInformation *info) {
  if(port == 0)
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
  else
    return 0;
  return 1;
}

int ttkPlanarGraphLayout::RequestData(vtkInformation *request,
                                      vtkInformationVector **inputVector,
                                      vtkInformationVector *outputVector) {
  // Get input and output objects
  auto input = vtkUnstructuredGrid::GetData(inputVector[0]);
  auto output = vtkUnstructuredGrid::GetData(outputVector);

  if(input == nullptr || output == nullptr) {
    return -1;
  }

  // Copy input to output
  output->ShallowCopy(input);

  size_t nPoints = output->GetNumberOfPoints();
  size_t nEdges = output->GetNumberOfCells();

  // Get input arrays
  auto sequenceArray = this->GetInputArrayToProcess(0, inputVector);
  if(this->UseSequences && !sequenceArray)
    return !this->printErr("Unable to retrieve sequence array.");

  auto sizeArray = this->GetInputArrayToProcess(1, inputVector);
  if(this->UseSizes && (!sizeArray || sizeArray->GetDataType()!=VTK_FLOAT))
    return !this->printErr("Unable to retrieve float size array.");

  auto branchArray = this->GetInputArrayToProcess(2, inputVector);
  if(this->UseBranches && !branchArray)
    return !this->printErr("Unable to retrieve branch array.");

  auto levelArray = this->GetInputArrayToProcess(3, inputVector);
  if(this->UseLevels && !levelArray)
    return !this->printErr("Unable to retrieve level array.");

  // Initialize output array
  auto outputArray = vtkSmartPointer<vtkFloatArray>::New();
  outputArray->SetName(this->GetOutputArrayName().data());
  outputArray->SetNumberOfComponents(2); // (x,y) position
  outputArray->SetNumberOfValues(nPoints * 2);

  auto DT = this->UseSequences ? sequenceArray->GetDataType() : VTK_INT;
  auto IT = this->UseBranches ? branchArray->GetDataType()
                : this->UseLevels ? levelArray->GetDataType() : VTK_INT;



  int status = 1;
  if(this->Algorithm==ALGORITHM::DOT){
    ttkTypeMacroAI(DT,IT,
      (
        status = this->computeGraphLayout<T0, T1>(
          // Output
          ttkUtils::GetPointer<float>(outputArray),
          // Input
          output->GetCells()->GetData()->GetPointer(0),
          nPoints, nEdges,
          !this->UseSequences ? nullptr : ttkUtils::GetPointer<T0>(sequenceArray),
          !this->UseSizes ? nullptr : ttkUtils::GetPointer<float>(sizeArray),
          !this->UseBranches ? nullptr : ttkUtils::GetPointer<T1>(branchArray),
          !this->UseLevels ? nullptr : ttkUtils::GetPointer<T1>(levelArray)
        )
      )
    );
  } else {
    if(!this->UseSequences || !this->UseBranches)
      return !this->printErr("Merge Tree Layout requires Sequence and Branch Arrays.");

    ttkTypeMacroAI(DT,IT,
      (
        status = this->computeMergeTreeLayout<T0, T1>(
          // Output
          ttkUtils::GetPointer<float>(outputArray),
          // Input
          output->GetCells()->GetData()->GetPointer(0),
          nPoints, nEdges,
          ttkUtils::GetPointer<T0>(sequenceArray),
          ttkUtils::GetPointer<T1>(branchArray)
        )
      )
    );
  }

  if(status != 1)
    return 0;

  // Add output field to output
  output->GetPointData()->AddArray(outputArray);

  return 1;
}
