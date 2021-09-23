#include <ttkPerlinPathIntegrator.h>

#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkPerlinPathIntegrator);

ttkPerlinPathIntegrator::ttkPerlinPathIntegrator() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkPerlinPathIntegrator::~ttkPerlinPathIntegrator() {
}

int ttkPerlinPathIntegrator::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    return 1;
  }
  return 0;
}

int ttkPerlinPathIntegrator::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkPerlinPathIntegrator::RequestData(vtkInformation *request,
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {

  // Get input object from input vector
  auto input = vtkPolyData::GetData(inputVector[0]);
  if(!input) {
    this->printErr("No valid vtkPolyData input.");
    return 0;
  }

  // Clear previous data
  pointsPerTimestep.clear();

  // Check if input has a domain dimension vector we can use for
  // the integration
  auto dimArray = GetInputArrayToProcess(1, input);
  if(!dimArray) {
    this->printErr("No field array with dimension data was provided.");
    return 0;
  }
  int dims[3];
  dims[0] = dimArray->GetTuple(0)[0];
  dims[1] = dimArray->GetTuple(0)[1];
  dims[2] = dimArray->GetTuple(0)[2];

  // Check if input has a pointId array. If there exists none, we'll create our
  // own
  auto idArray = GetInputArrayToProcess(0, input);
  if(!idArray) {
    this->printErr("No id array with point ids was provided.");
    return 0;
  }

  // Initialize size of vector of points per timestep
  pointsPerTimestep.resize(this->nTimesteps);

  // Create initial vector of points for timestep 0 from input seed points
  int nInitPoints = input->GetNumberOfPoints();
  auto initPoints = input->GetPoints();
  pointsPerTimestep[0].resize(nInitPoints);
  for (int i=0; i < nInitPoints; i++) {
    auto& p =  pointsPerTimestep[0][i];
    p.pointId = idArray->GetTuple(i)[0];
    p.timestep = 0;
    double pos[3];
    initPoints->GetPoint(i, pos);
    p.x = pos[0];
    p.y = pos[1];
    p.z = pos[2];
  }

  // Call base-layer with the input seed points, the output data structure,
  // and the dimension of the input domain, and time info
  int status = 0;
  switch(idArray->GetDataType()) {
    vtkTemplateMacro(
      status = this->integrate<VTK_TT>(
        pointsPerTimestep[0],
        pointsPerTimestep,
        this->nTimesteps,
        this->TimeInterval,
        dims,
        this->StepLength,
        this->PerlinScaleFactor
      )
    );
  }
  if (!status) {
    this->printErr("Integration could not be executed");
    return 0;
  }

  // Create all new points and arrays
  int nIntegratedPoints = 0;
  for (int i = 1; i < this->nTimesteps; i++) {
    nIntegratedPoints += pointsPerTimestep[i].size();
  }
  auto newPoints = vtkSmartPointer<vtkPoints>::New();
  newPoints->SetDataType(VTK_DOUBLE);
  newPoints->SetNumberOfPoints(nInitPoints+nIntegratedPoints);

  auto prepArray = [](vtkDataArray* array, std::string name, int nTuples, int nComponents){
    array->SetName(name.data());
    array->SetNumberOfComponents(nComponents);
    array->SetNumberOfTuples(nTuples);
    return ttkUtils::GetVoidPointer(array);
  };

  auto timeArray = vtkSmartPointer<vtkIntArray>::New();
  auto timeArrayData = static_cast<int*>(prepArray(timeArray,"Timestep",(nInitPoints+nIntegratedPoints),1));
  auto id2Array = vtkSmartPointer<vtkIntArray>::New();
  auto id2ArrayData = static_cast<int*>(prepArray(id2Array,"PointId",(nInitPoints+nIntegratedPoints),1));

  // Add data to points and arrays
  int idx = 0;
  for (int i = 0; i < this->nTimesteps; i++) {
    for (size_t j = 0; j < pointsPerTimestep[i].size(); j++) {
      auto p = pointsPerTimestep[i][j];
      double pos[3] = {p.x, p.y, p.z};
      newPoints->SetPoint(idx, pos);
      timeArrayData[idx] = p.timestep;
      id2ArrayData[idx] = p.pointId;
      idx++;
    }
  } 

  // Format the output data structure into a dataset of vtkPolyData
  auto output = vtkPolyData::GetData(outputVector);
  //auto points = output->GetPoints();
  output->SetPoints(newPoints);
  //idArray->InsertTuples(nInitPoints, nIntegratedPoints,nInitPoints, id2Array);
  output->GetPointData()->AddArray(id2Array);
  output->GetPointData()->AddArray(timeArray);


  // return success
  return 1;
}
