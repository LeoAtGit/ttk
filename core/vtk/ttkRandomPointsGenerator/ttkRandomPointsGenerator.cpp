#include <ttkRandomPointsGenerator.h>

#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointSet.h>
#include <vtkPolyData.h>
#include <vtkIntArray.h>


#include <ttkMacros.h>
#include <ttkUtils.h>

// std includes
#include <random>

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkRandomPointsGenerator);

ttkRandomPointsGenerator::ttkRandomPointsGenerator() {
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
}

ttkRandomPointsGenerator::~ttkRandomPointsGenerator() {
}

int ttkRandomPointsGenerator::FillInputPortInformation(int port, vtkInformation *info) {
  return 0;
}

int ttkRandomPointsGenerator::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
    return 1;
  }
  return 0;
}


int ttkRandomPointsGenerator::RequestData(vtkInformation *request,
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {
  
  // Initialize points
  auto points = vtkSmartPointer<vtkPoints>::New();
  points->SetDataType(VTK_DOUBLE);
  points->SetNumberOfPoints(nPoints);

  // Get random number engine and seed it with user provided seed
  std::mt19937 posGen(RandomSeed);

  std::uniform_real_distribution<> disX(0.0, (double)PositionDomain[0]);
  std::uniform_real_distribution<> disY(0.0, (double)PositionDomain[1]);
  std::uniform_real_distribution<> disZ(0.0, (double)PositionDomain[2]);

  // Create array with point Ids for viz purposes
  vtkSmartPointer<vtkIntArray> idArray = vtkIntArray::New();
  idArray->SetNumberOfComponents(1);
  idArray->SetNumberOfTuples(nPoints);
  idArray->SetName("PointId");

  for (int i = 0; i < nPoints; i++) {
    double pos[3] = {disX(posGen), disY(posGen), disZ(posGen)};
    points->SetPoint(i, pos);
    //int data = i;
    idArray->SetTuple1(i, i);
  }

  // Get output image
  auto output = vtkPolyData::GetData(outputVector);
  output->SetPoints(points);
  output->GetPointData()->AddArray(idArray);


  // return success
  return 1;
}
