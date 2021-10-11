#include <ttkScalarFieldFromPoints.h>


#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkFloatArray.h>
#include <vtkDataSet.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkInformationVector.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkPointSet.h>
#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkPolyData.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkScalarFieldFromPoints);

ttkScalarFieldFromPoints::ttkScalarFieldFromPoints() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkScalarFieldFromPoints::~ttkScalarFieldFromPoints() {
}

int ttkScalarFieldFromPoints::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0)
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
  else
    return 0;

  return 1;
}

int ttkScalarFieldFromPoints::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0)
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
  else
    return 0;

  return 1;
}

int ttkScalarFieldFromPoints::RequestInformation(vtkInformation *,
                                   vtkInformationVector **,
                                   vtkInformationVector * outputVector) {

  return 1;
}

int ttkScalarFieldFromPoints::RequestData(vtkInformation *request,
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {
  
  auto inputMB = vtkMultiBlockDataSet::GetData(inputVector[0]);
  if(!inputMB) {
    this->printErr("No input data available.");
    return 0;
  }

  // Get points for current timestep
  auto tBlock = inputMB->GetBlock(0);
  auto curPoints = vtkSmartPointer<vtkPolyData>::New();
  curPoints->ShallowCopy(tBlock);
  const size_t nPoints = curPoints->GetNumberOfPoints();

  auto image = vtkSmartPointer<vtkImageData>::New();
  image->SetDimensions(
    this->Resolution[0],
    this->Resolution[1],
    this->Resolution[2]
  );
  image->SetOrigin(
    this->ImageBounds[0],
    this->ImageBounds[2],
    this->ImageBounds[4]
  );
  image->SetSpacing(
    this->Resolution[0]>1 ? (this->ImageBounds[1]-this->ImageBounds[0])/(this->Resolution[0]-1) : 0,
    this->Resolution[1]>1 ? (this->ImageBounds[3]-this->ImageBounds[2])/(this->Resolution[1]-1) : 0,
    this->Resolution[2]>1 ? (this->ImageBounds[5]-this->ImageBounds[4])/(this->Resolution[2]-1) : 0
  );
  image->AllocateScalars(VTK_DOUBLE, 1);

  auto scalarArray = image->GetPointData()->GetArray(0);
  scalarArray->SetName("Scalars");
  auto scalarArrayData = ttkUtils::GetPointer<double>(scalarArray);

  // auto nPixels = scalarArray->GetNumberOfTuples();

  // Get a triangulation
  ttk::Triangulation *triangulation = ttkAlgorithm::GetTriangulation(image);
  if(!triangulation) {
    this->printErr("No triangulation could be created.");
    return 0;
  }

  int status = 0;
  switch(this->Kernel){
    case 0: {
      ttkVtkTemplateMacro(scalarArray->GetDataType(), triangulation->getType(),
      (status = this->computeKDE<TTK_TT, ScalarFieldFromPoints::Gaussian>(
        scalarArrayData,
        ttkUtils::GetPointer<double>(curPoints->GetPoints()->GetData()),
        nPoints,
        this->Bandwidth,
        (TTK_TT *)triangulation->getData()
      )));
      break;
    }
    case 1: {
      ttkVtkTemplateMacro(scalarArray->GetDataType(), triangulation->getType(),
      (status = this->computeKDE<TTK_TT, ScalarFieldFromPoints::Linear>(
        scalarArrayData,
        ttkUtils::GetPointer<double>(curPoints->GetPoints()->GetData()),
        nPoints,
        this->Bandwidth,
        (TTK_TT *)triangulation->getData()
      )));
      break;
    }
    case 2: {
      ttkVtkTemplateMacro(scalarArray->GetDataType(), triangulation->getType(),
      (status = this->computeKDE<TTK_TT, ScalarFieldFromPoints::Epanechnikov>(
        scalarArrayData,
        ttkUtils::GetPointer<double>(curPoints->GetPoints()->GetData()),
        nPoints,
        this->Bandwidth,
        (TTK_TT *)triangulation->getData()
      )));
      break;
    }
  }

  // On error cancel filter execution
  if(status == 0)
    return 0;
  
  auto outputMB = vtkMultiBlockDataSet::GetData(outputVector);
  // Set image to a block in the output dataset
  size_t nBlocks = outputMB->GetNumberOfBlocks();
  outputMB->SetBlock(nBlocks, image);


  return 1;
}
