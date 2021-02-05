#include <ttkPlanarGraphLayout.h>

#include <ttkMacros.h>

#include <vtkAbstractArray.h>
#include <vtkCellArray.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include <ttkUtils.h>

#define ttkTypeMacroErrorCase(idx, type)                          \
  default: {                                                      \
    this->printErr("Unsupported " #idx "-th Template Data Type: " \
                   + std::to_string(type));                       \
  } break;

#define ttkTypeMacroCase(enum, type, number, call) \
  case enum: {                                     \
    typedef type T##number;                        \
    call;                                          \
  } break;

#define ttkTypeMacroR(target, call)                                 \
  switch(target) {                                                  \
    ttkTypeMacroCase(VTK_FLOAT, float, 0, call) ttkTypeMacroCase(   \
      VTK_DOUBLE, double, 0, call) ttkTypeMacroErrorCase(0, target) \
  }

#define ttkTypeMacroA(target, call)                                        \
  switch(target) {                                                         \
    ttkTypeMacroCase(VTK_FLOAT, float, 0, call);                           \
    ttkTypeMacroCase(VTK_DOUBLE, double, 0, call);                         \
    ttkTypeMacroCase(VTK_INT, int, 0, call);                               \
    ttkTypeMacroCase(VTK_UNSIGNED_INT, unsigned int, 0, call);             \
    ttkTypeMacroCase(VTK_CHAR, char, 0, call);                             \
    ttkTypeMacroCase(VTK_SIGNED_CHAR, signed char, 0, call);               \
    ttkTypeMacroCase(VTK_UNSIGNED_CHAR, unsigned char, 0, call);           \
    ttkTypeMacroCase(VTK_LONG, long, 0, call);                             \
    ttkTypeMacroCase(VTK_LONG_LONG, long long, 0, call);                   \
    ttkTypeMacroCase(VTK_UNSIGNED_LONG, unsigned long, 0, call);           \
    ttkTypeMacroCase(VTK_UNSIGNED_LONG_LONG, unsigned long long, 0, call); \
    ttkTypeMacroCase(VTK_ID_TYPE, vtkIdType, 0, call);                     \
    ttkTypeMacroErrorCase(0, target);                                      \
  }

#define ttkTypeMacroRR(target0, target1, call)                             \
  switch(target1) {                                                        \
    ttkTypeMacroCase(VTK_FLOAT, float, 1, ttkTypeMacroR(target0, call));   \
    ttkTypeMacroCase(VTK_DOUBLE, double, 1, ttkTypeMacroR(target0, call)); \
    ttkTypeMacroErrorCase(1, target1);                                     \
  }

#define ttkTypeMacroRRR(target0, target1, target2, call)              \
  switch(target2) {                                                   \
    ttkTypeMacroCase(                                                 \
      VTK_FLOAT, float, 2, ttkTypeMacroRR(target0, target1, call));   \
    ttkTypeMacroCase(                                                 \
      VTK_DOUBLE, double, 2, ttkTypeMacroRR(target0, target1, call)); \
    ttkTypeMacroErrorCase(2, target2);                                \
  }

#define ttkTypeMacroRRI(target0, target1, target2, call)                       \
  switch(target2) {                                                            \
    ttkTypeMacroCase(VTK_INT, int, 2, ttkTypeMacroRR(target0, target1, call)); \
    ttkTypeMacroCase(                                                          \
      VTK_LONG_LONG, long long, 2, ttkTypeMacroRR(target0, target1, call));    \
    ttkTypeMacroCase(                                                          \
      VTK_ID_TYPE, vtkIdType, 2, ttkTypeMacroRR(target0, target1, call));      \
    ttkTypeMacroErrorCase(2, target2);                                         \
  }

#define ttkTypeMacroAI(target0, target1, call)                                 \
  switch(target1) {                                                            \
    ttkTypeMacroCase(VTK_INT, int, 1, ttkTypeMacroA(target0, call));           \
    ttkTypeMacroCase(                                                          \
      VTK_LONG_LONG, long long, 1, ttkTypeMacroA(target0, call));              \
    ttkTypeMacroCase(VTK_ID_TYPE, vtkIdType, 1, ttkTypeMacroA(target0, call)); \
    ttkTypeMacroErrorCase(1, target1);                                         \
  }

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

  if(!input || !output)
    return 0;

  // Copy input to output
  output->ShallowCopy(input);

  size_t nPoints = output->GetNumberOfPoints();
  size_t nEdges = output->GetNumberOfCells();

  // Get input arrays
  auto sequenceArray = this->GetInputArrayToProcess(0, inputVector);
  if(this->GetUseSequences() && !sequenceArray) {
    this->printErr("Unable to retrieve sequence array.");
    return 0;
  }

  auto sizeArray = this->GetInputArrayToProcess(1, inputVector);
  if(this->GetUseSizes() && (!sizeArray || sizeArray->GetDataType()!=VTK_FLOAT)) {
    this->printErr("Unable to retrieve size array of type float.");
    return 0;
  }

  auto branchArray = this->GetInputArrayToProcess(2, inputVector);
  if(this->GetUseBranches() && (!branchArray || branchArray->GetDataType()!=VTK_INT)) {
    this->printErr("Unable to retrieve branch array of type int.");
    return 0;
  }

  auto levelArray = this->GetInputArrayToProcess(3, inputVector);
  if(this->GetUseLevels() && (!levelArray || levelArray->GetDataType()!=VTK_INT)) {
    this->printErr("Unable to retrieve level array of type int.");
    return 0;
  }

  // Initialize output array
  auto outputArray = vtkSmartPointer<vtkFloatArray>::New();
  outputArray->SetName(this->GetOutputArrayName().data());
  outputArray->SetNumberOfComponents(2); // (x,y) position
  outputArray->SetNumberOfValues(nPoints * 2);

  int status = 1;

  auto inputConnectivityList = input->GetCells()->GetConnectivityArray();

  ttkTypeMacroAI(
    this->GetUseSequences() ? sequenceArray->GetDataType() : VTK_INT,
    inputConnectivityList->GetDataType(),
    (status = this->execute<T1,T0>(
    // Output
    static_cast<float *>(ttkUtils::GetVoidPointer(outputArray)),

    // Input
    static_cast<T1*>(ttkUtils::GetVoidPointer(inputConnectivityList)),
    nPoints, nEdges,
    static_cast<T0*>(this->GetUseSequences() ? ttkUtils::GetVoidPointer(sequenceArray) : nullptr),
    static_cast<float*>(this->GetUseSizes() ? ttkUtils::GetVoidPointer(sizeArray) : nullptr),
    static_cast<int*>(this->GetUseBranches() ? ttkUtils::GetVoidPointer(branchArray) : nullptr),
    static_cast<int*>(this->GetUseLevels() ? ttkUtils::GetVoidPointer(levelArray) : nullptr)
  )));
  if(status != 1)
    return 0;

  // Add output field to output
  output->GetPointData()->AddArray(outputArray);

  return 1;
}