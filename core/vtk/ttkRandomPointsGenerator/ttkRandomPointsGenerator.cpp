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

  // Position distributions
  std::uniform_real_distribution<> disX(0.0, (double)this->PositionDomain[0]);
  std::uniform_real_distribution<> disY(0.0, (double)this->PositionDomain[1]);
  std::uniform_real_distribution<> disZ(0.0, (double)this->PositionDomain[2]);

  // Point data distributions
  std::uniform_real_distribution<> disAmp(this->Amplitude[0], this->Amplitude[1]);
  std::uniform_real_distribution<> disSpread(this->Spread[0], this->Spread[1]);
  std::uniform_real_distribution<> disRate(0, 1);


  // Function for formatting data arrays
  auto prepArray = [](vtkDataArray* array, std::string name, int nTuples, int nComponents){
    array->SetName(name.data());
    array->SetNumberOfComponents(nComponents);
    array->SetNumberOfTuples(nTuples);
  };  

  // Create an array to store dimensions in
  auto dimArray = vtkSmartPointer<vtkIntArray>::New();
  prepArray(dimArray, "DomainDimension", 1, 3);
  dimArray->SetTuple3(0, PositionDomain[0], PositionDomain[1], PositionDomain[2]);

  // Create array with point Ids for viz purposes
  auto idArray = vtkSmartPointer<vtkIntArray>::New();
  prepArray(idArray, "PointId", nPoints, 1);

  // Create amplitude, spread and rate array to be used when creating scalar fields
  auto ampArray = vtkSmartPointer<vtkDoubleArray>::New();
  prepArray(ampArray, "Amplitude", nPoints, 1);
  auto spreadArray = vtkSmartPointer<vtkDoubleArray>::New();
  prepArray(spreadArray, "Spread", nPoints, 1);
  auto rateArray = vtkSmartPointer<vtkDoubleArray>::New();
  prepArray(rateArray, "Rate", nPoints, 1);

  for (int i = 0; i < nPoints; i++) {
    double pos[3] = {disX(posGen), disY(posGen), disZ(posGen)};
    points->SetPoint(i, pos);

    // Set point data
    idArray->SetTuple1(i, i);
    ampArray->SetTuple1(i, disAmp(posGen));
    spreadArray->SetTuple1(i, disSpread(posGen));
    rateArray->SetTuple1(i, disRate(posGen));
  }

  // Get output data and add arrays
  auto output = vtkPolyData::GetData(outputVector);
  output->SetPoints(points);
  auto outputPD = output->GetPointData();
  outputPD->AddArray(idArray);
  outputPD->AddArray(ampArray);
  outputPD->AddArray(spreadArray);
  outputPD->AddArray(rateArray);
  output->GetFieldData()->AddArray(dimArray);

  // return success
  return 1;
}
