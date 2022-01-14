#include <ttkPathIntegrator.h>

#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

// std includes
#include <limits>
#include <random>

vtkStandardNewMacro(ttkPathIntegrator);

ttkPathIntegrator::ttkPathIntegrator() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkPathIntegrator::~ttkPathIntegrator() {
}

int ttkPathIntegrator::FillInputPortInformation(int port,
                                                vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    return 1;
  }
  return 0;
}

int ttkPathIntegrator::FillOutputPortInformation(int port,
                                                 vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

int ttkPathIntegrator::RequestData(vtkInformation *request,
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

  // Check if all required input arrays exist
  auto idArray = GetInputArrayToProcess(0, input);
  if(!idArray) {
    this->printErr("No id array with point ids was provided.");
    return 0;
  }

  auto ampArray = GetInputArrayToProcess(2, input);
  if(!ampArray) {
    this->printErr("No amplitude array was provided.");
    return 0;
  }

  auto spreadArray = GetInputArrayToProcess(3, input);
  if(!spreadArray) {
    this->printErr("No spread array was provided.");
    return 0;
  }

  auto rateArray = GetInputArrayToProcess(4, input);
  if(!rateArray) {
    this->printErr("No rate array was provided.");
    return 0;
  }

  // Initialize size of vector of points per timestep
  pointsPerTimestep.resize(this->nTimesteps);

  // Create initial vector of points for timestep 0 from input seed points,
  // and find largest spread
  int nInitPoints = input->GetNumberOfPoints();
  auto initPoints = input->GetPoints();
  pointsPerTimestep[0].resize(nInitPoints);
  double maxSpread = std::numeric_limits<double>::min();
  for(int i = 0; i < nInitPoints; i++) {
    auto &p = pointsPerTimestep[0][i];
    p.pointId = idArray->GetTuple(i)[0];
    p.amplitude = ampArray->GetTuple(i)[0];
    p.rate = rateArray->GetTuple(i)[0];
    p.timestep = 0;
    p.spread = spreadArray->GetTuple(i)[0];

    // Check if spread is over max, if so save
    if(p.spread > maxSpread)
      maxSpread = p.spread;

    double pos[3];
    initPoints->GetPoint(i, pos);
    p.x = pos[0];
    p.y = pos[1];
    p.z = pos[2];
  }

  // Determine what vector field to use
  PathIntegrator::VectorField vf = PathIntegrator::VectorField::PerlinPerturbed;
  ;
  switch(this->VecField) {
    case 0: {
      vf = PathIntegrator::VectorField::PerlinPerturbed;
      break;
    }
    case 1: {
      vf = PathIntegrator::VectorField::PerlinGradient;
      break;
    }
    case 2: {
      vf = PathIntegrator::VectorField::PosDiagonal;
      break;
    }
    case 3: {
      vf = PathIntegrator::VectorField::PosX;
    }
  }

  // Call base-layer with the input seed points, the output data structure,
  // and the dimension of the input domain, and time info
  int status = 0;
  switch(idArray->GetDataType()) {
    vtkTemplateMacro(status = this->integrate<VTK_TT>(
                       pointsPerTimestep, this->nTimesteps, this->TimeInterval,
                       dims, this->StepLength, this->PerlinScaleFactor,
                       maxSpread, vf));
  }
  if(!status) {
    this->printErr("Integration could not be executed");
    return 0;
  }

  // Create output
  auto outputMB = vtkMultiBlockDataSet::GetData(outputVector);

  // Function for formatting data arrays
  auto prepArray
    = [](vtkDataArray *array, std::string name, int nTuples, int nComponents) {
        array->SetName(name.data());
        array->SetNumberOfComponents(nComponents);
        array->SetNumberOfTuples(nTuples);
        return ttkUtils::GetVoidPointer(array);
      };

  // Go through all timesteps and create output vtkPolyData
  for(int i = 0; i < this->nTimesteps; i++) {
    // Create vtkPolyData for current timestep
    auto polyData = vtkSmartPointer<vtkPolyData>::New();
    auto dataPoints = vtkSmartPointer<vtkPoints>::New();
    int nPoints = pointsPerTimestep[i].size();

    dataPoints->SetDataType(VTK_DOUBLE);
    dataPoints->SetNumberOfPoints(nPoints);

    // Create cell arrays off offset array and connectivity array.
    // Our cells are just vertices
    auto offsetArray = vtkSmartPointer<vtkIntArray>::New();
    offsetArray->SetNumberOfTuples(nPoints + 1);
    auto offsetArrayData
      = static_cast<int *>(ttkUtils::GetVoidPointer(offsetArray));

    auto connectivityArray = vtkSmartPointer<vtkIntArray>::New();
    connectivityArray->SetNumberOfTuples(nPoints);
    auto connectivityArrayData
      = static_cast<int *>(ttkUtils::GetVoidPointer(connectivityArray));

    auto cellArray = vtkSmartPointer<vtkCellArray>::New();
    cellArray->SetData(offsetArray, connectivityArray);

    // Create arrays for output data
    auto timeArray = vtkSmartPointer<vtkIntArray>::New();
    auto timeArrayData
      = static_cast<int *>(prepArray(timeArray, "Timestep", nPoints, 1));
    auto id2Array = vtkSmartPointer<vtkIntArray>::New();
    auto id2ArrayData
      = static_cast<int *>(prepArray(id2Array, "PointId", nPoints, 1));
    auto amp2Array = vtkSmartPointer<vtkDoubleArray>::New();
    auto amp2ArrayData
      = static_cast<double *>(prepArray(amp2Array, "Amplitude", nPoints, 1));
    auto spread2Array = vtkSmartPointer<vtkDoubleArray>::New();
    auto spread2ArrayData
      = static_cast<double *>(prepArray(spread2Array, "Spread", nPoints, 1));
    auto rate2Array = vtkSmartPointer<vtkDoubleArray>::New();
    auto rate2ArrayData
      = static_cast<double *>(prepArray(rate2Array, "Rate", nPoints, 1));
    auto velocityArray = vtkSmartPointer<vtkDoubleArray>::New();
    auto velocityArrayData
      = static_cast<double *>(prepArray(velocityArray, "Velocity", nPoints, 3));
    auto birthTimeArray = vtkSmartPointer<vtkIntArray>::New();
    auto birthTimeArrayData
      = static_cast<int *>(prepArray(birthTimeArray, "BirthTime", nPoints, 1));
    auto deathTimeArray = vtkSmartPointer<vtkIntArray>::New();
    auto deathTimeArrayData
      = static_cast<int *>(prepArray(deathTimeArray, "DeathTime", nPoints, 1));

    // Create a random engine which is used to randomly sample different birth
    // and death times for the points
    std::mt19937 randGen(1);
    std::uniform_int_distribution<> disBirth(
      -floor(this->nTimesteps / 2), this->nTimesteps);
    std::uniform_int_distribution<> disDeath(0, this->nTimesteps);

    // Add data to points and arrays
    int idx = 0;
    for(int j = 0; j < nPoints; j++) {
      auto p = pointsPerTimestep[i][j];

      // Set position
      double pos[3] = {p.x, p.y, p.z};
      dataPoints->SetPoint(idx, pos);

      // Set data arrays
      timeArrayData[idx] = p.timestep;
      id2ArrayData[idx] = p.pointId;
      amp2ArrayData[idx] = p.amplitude;
      spread2ArrayData[idx] = p.spread;
      rate2ArrayData[idx] = p.rate;

      // Generate birth and death data randomly
      int birth = disBirth(randGen);
      int death = disDeath(randGen);
      while(death < birth) {
        death = disDeath(randGen);
      }
      birthTimeArrayData[idx] = birth;
      deathTimeArrayData[idx] = death;

      // Store velocity
      velocityArrayData[3 * idx] = p.v[0];
      velocityArrayData[3 * idx + 1] = p.v[1];
      velocityArrayData[3 * idx + 2] = p.v[2];

      // Set connectivity and offset array
      connectivityArrayData[idx] = idx;
      offsetArrayData[idx] = idx;
      idx++;
    }

    // Set offset array last index to the number of elements in the connectivity
    // array
    offsetArrayData[nPoints] = nPoints;

    // Format the output data structure into a dataset of vtkPolyData
    polyData->SetPoints(dataPoints);
    polyData->SetVerts(cellArray);

    auto pointData = polyData->GetPointData();
    pointData->AddArray(id2Array);
    pointData->AddArray(amp2Array);
    pointData->AddArray(spread2Array);
    pointData->AddArray(rate2Array);
    pointData->AddArray(timeArray);
    pointData->AddArray(velocityArray);
    pointData->AddArray(birthTimeArray);
    pointData->AddArray(deathTimeArray);

    // Set data to a block in the output dataset
    size_t nBlocks = outputMB->GetNumberOfBlocks();
    outputMB->SetBlock(nBlocks, polyData);
  }

  // return success
  return 1;
}
