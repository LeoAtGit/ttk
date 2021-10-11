#include <ttkPerlinPathIntegrator.h>

#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkDoubleArray.h>
#include <vtkMultiBlockDataSet.h>

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
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
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

  // Create output
  auto outputMB = vtkMultiBlockDataSet::GetData(outputVector);

  // Function for formatting data arrays
  auto prepArray = [](vtkDataArray* array, std::string name, int nTuples, int nComponents){
    array->SetName(name.data());
    array->SetNumberOfComponents(nComponents);
    array->SetNumberOfTuples(nTuples);
    return ttkUtils::GetVoidPointer(array);
  };  

  for (int i = 0; i < this->nTimesteps; i++) {
    // Create vtkPolyData for current timestep
    auto polyData = vtkSmartPointer<vtkPolyData>::New();
    auto dataPoints = vtkSmartPointer<vtkPoints>::New();
    int nPoints = pointsPerTimestep[i].size();

    dataPoints->SetDataType(VTK_DOUBLE);
    dataPoints->SetNumberOfPoints(nPoints);

    // Create cell arrays off offset array and connectivity array.
    // Our cells are just vertices
    auto offsetArray = vtkSmartPointer<vtkIntArray>::New();
    offsetArray->SetNumberOfTuples(nPoints+1);
    auto offsetArrayData = static_cast<int*>(ttkUtils::GetVoidPointer(offsetArray));

    auto connectivityArray = vtkSmartPointer<vtkIntArray>::New();
    connectivityArray->SetNumberOfTuples(nPoints);
    auto connectivityArrayData = static_cast<int*>(ttkUtils::GetVoidPointer(connectivityArray));

    auto cellArray = vtkSmartPointer<vtkCellArray>::New();
    cellArray->SetData(offsetArray, connectivityArray);

    auto timeArray = vtkSmartPointer<vtkIntArray>::New();
    auto timeArrayData = static_cast<int*>(prepArray(timeArray, "Timestep", nPoints, 1));
    auto id2Array = vtkSmartPointer<vtkIntArray>::New();
    auto id2ArrayData = static_cast<int*>(prepArray(id2Array, "PointId", nPoints, 1));
    auto velocityArray = vtkSmartPointer<vtkDoubleArray>::New();
    auto velocityArrayData = static_cast<double*>(prepArray(velocityArray,"Velocity", nPoints, 3));

    // Add data to points and arrays
    int idx = 0;
    for (int j = 0; j < nPoints; j++) {
      auto p = pointsPerTimestep[i][j];

      // Set position
      double pos[3] = {p.x, p.y, p.z};
      dataPoints->SetPoint(idx, pos);

      // Set data arrays
      timeArrayData[idx] = p.timestep;
      id2ArrayData[idx] = p.pointId;

      velocityArrayData[3 * idx] = p.v[0];
      velocityArrayData[3 * idx + 1] = p.v[1];
      velocityArrayData[3 * idx + 2] = p.v[2];
      
      // Set connectivity and offset array
      connectivityArrayData[idx] = idx;
      offsetArrayData[idx] = idx;
      idx++;
    }

    // Set offset array last index to the number of elements in the connectivity array
    offsetArrayData[nPoints] = nPoints;


    // Format the output data structure into a dataset of vtkPolyData
    polyData->SetPoints(dataPoints);
    polyData->SetVerts(cellArray);

    auto pointData = polyData->GetPointData();
    pointData->AddArray(id2Array);
    pointData->AddArray(timeArray);
    pointData->AddArray(velocityArray);

    // Set data to a block in the output dataset
    size_t nBlocks = outputMB->GetNumberOfBlocks();
    outputMB->SetBlock(nBlocks, polyData);
  }

  // return success
  return 1;
}
