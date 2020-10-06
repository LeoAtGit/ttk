#include <ttkCinemaImaging.h>

#include <vtkInformation.h>

#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>

#include <vtkPointData.h>
#include <vtkCellData.h>

#include <vtkDoubleArray.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkCinemaImaging);

ttkCinemaImaging::ttkCinemaImaging() {
  this->setDebugMsgPrefix("CinemaImaging");
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
};

ttkCinemaImaging::~ttkCinemaImaging() {
};

int ttkCinemaImaging::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0){
    info->Remove(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE());
    info->Append(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
    info->Append(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
  }
  else if(port == 1)
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  else
    return 0;
  return 1;
};

int ttkCinemaImaging::FillOutputPortInformation(int port,
                                                vtkInformation *info) {
  if(port == 0){
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
};

int addGobalArray(
    std::vector<vtkSmartPointer<vtkAbstractArray>> &globalArrays,
    std::string name,
    size_t nValues,
    const double *values
){
    auto array = vtkSmartPointer<vtkDoubleArray>::New();
    array->SetName(name.data());
    array->SetNumberOfComponents(nValues);
    array->SetNumberOfTuples(1);
    for(size_t i = 0; i < nValues; i++)
        array->SetValue(i, values[i]);
    globalArrays.push_back(array);

    return 1;
};

int ttkCinemaImaging::RequestData(vtkInformation *request,
                                  vtkInformationVector **inputVector,
                                  vtkInformationVector *outputVector) {
  ttk::Timer globalTimer;

  // ---------------------------------------------------------------------------
  // Get Input / Output
  // ---------------------------------------------------------------------------
  auto inputObject = vtkPointSet::GetData(inputVector[0]);
  auto inputGrid = vtkPointSet::GetData(inputVector[1]);
  auto outputImages = vtkMultiBlockDataSet::GetData(outputVector);

  // ---------------------------------------------------------------------------
  // Parameters
  // ---------------------------------------------------------------------------

  // get auto data
  double defaultAngle = this->CamAngle;
  double defaultHeight = this->CamHeight;
  std::vector<double> defaultFocus = {this->CamFocus[0],this->CamFocus[1],this->CamFocus[2]};
  std::vector<double> defaultNearFar = {this->CamNearFar[0],this->CamNearFar[1]};

  if(this->GetCamFocusAuto() || this->GetCamNearFarAuto()
    || this->GetCamHeightAuto()) {
    double dxO, dyO, dzO, dO; // object
    double dxG, dyG, dzG, dG; // object

    auto getBoxDiameter = [](
      const double *bounds, double &dx, double &dy, double &dz, double &d) {
      dx = bounds[1] - bounds[0];
      dy = bounds[3] - bounds[2];
      dz = bounds[5] - bounds[4];
      d = std::sqrt(dx * dx + dy * dy + dz * dz);
      return 1;
    };

    double gridBounds[6];
    inputGrid->GetBounds(gridBounds);
    getBoxDiameter(gridBounds, dxG, dyG, dzG, dG);

    double objectBounds[6];
    inputObject->GetBounds(objectBounds);
    getBoxDiameter(objectBounds, dxO, dyO, dzO, dO);

    if(this->GetCamFocusAuto()) {
      defaultFocus[0] = objectBounds[0] + dxO * 0.5;
      defaultFocus[1] = objectBounds[2] + dyO * 0.5;
      defaultFocus[2] = objectBounds[4] + dzO * 0.5;
    }

    if(this->GetCamNearFarAuto()) {
      double gDiameter = std::max(dxG, std::max(dyG, dzG));
      defaultNearFar[0] = gDiameter / 2 - dO / 2;
      // Ensure that camNear is not 0
      defaultNearFar[0] = std::max(defaultNearFar[0], 0.00001);
      defaultNearFar[1] = defaultNearFar[0] + dO;
    }

    if(this->GetCamHeightAuto()) {
      defaultHeight = dO;
    }
  }

  // ensure grid contains all camera parameters as field data
  auto aInputGrid = vtkSmartPointer<vtkPointSet>::Take( inputGrid->NewInstance() );
  aInputGrid->ShallowCopy(inputGrid);
  int n = aInputGrid->GetNumberOfPoints();
  auto aInputGridPD = aInputGrid->GetPointData();
  ttkCinemaImaging::ensureGridData( aInputGridPD, "Resolution", n, {(double)this->Resolution[0],(double)this->Resolution[1]} );
  ttkCinemaImaging::ensureGridData( aInputGridPD, "ProjectionMode", n, {(double)this->CamProjectionMode} );
  if(aInputGridPD->HasArray("CamDir"))
    ttkCinemaImaging::ensureGridData( aInputGridPD, "CamDir", n, {0,0,-1} );
  else {
    ttkCinemaImaging::ensureGridData( aInputGridPD, "CamFocus", n, defaultFocus );
    ttkCinemaImaging::computeDirFromFocus(aInputGrid);
  }
  ttkCinemaImaging::ensureGridData( aInputGridPD, "CamNearFar", n, defaultNearFar );
  ttkCinemaImaging::ensureGridData( aInputGridPD, "CamHeight", n, {defaultHeight} );
  ttkCinemaImaging::ensureGridData( aInputGridPD, "CamAngle", n, {defaultAngle} );
  ttkCinemaImaging::ensureGridData( aInputGridPD, "CamUp", n, {0,1,0} );

  // print parameters
  {
    this->printMsg({
      {"Backend", this->Backend==0 ? std::string("VTK_OPENGL") : std::string("EMBREE")},
      {"Resolution", std::to_string(this->Resolution[0]) + " x "
                      + std::to_string(this->Resolution[1])},
      {"Projection", this->GetCamProjectionMode() ? "Perspective" : "Orthographic"},
      {"CamNearFar",
      std::to_string(defaultNearFar[0]) + " / " + std::to_string(defaultNearFar[1])},
      {"CamFocus", "(" + std::to_string(defaultFocus[0]) + ", "
                    + std::to_string(defaultFocus[1]) + ", "
                    + std::to_string(defaultFocus[2]) + ")"},
      {
          std::string(this->CamProjectionMode==0 ? "CamHeight" : "CamAngle"),
          std::to_string(this->CamProjectionMode==0 ? defaultHeight : defaultAngle)
      }
    });
    this->printMsg(ttk::debug::Separator::L1);
  }

  vtkCellArray* inputObjectCells = nullptr;
  if(inputObject->IsA("vtkPolyData")){
      inputObjectCells = static_cast<vtkPolyData*>(inputObject)->GetPolys();
  } else if(inputObject->IsA("vtkUnstructuredGrid")){
      inputObjectCells = static_cast<vtkUnstructuredGrid*>(inputObject)->GetCells();
  } else {
      this->printErr("Input has to be a vtkPolyData or vtkUnstructuredGrid.");
      return 0;
  }

  size_t nTriangles = inputObjectCells->GetNumberOfCells();
  // make sure that cells consists only of triangles
  {
      auto offsets = static_cast<vtkIdType*>(
          ttkUtils::GetVoidPointer(
              inputObjectCells->GetOffsetsArray()
          )
      );
      vtkIdType j=0;
      for(size_t i=0; i<nTriangles; i++, j+=3){
        if(j!=offsets[i]){
            this->printErr("Connectivity list has to contain only triangles.");
            return 0;
        }
      }
  }

  int status = 0;
  if(this->Backend==0){
//       ttk::ttkCinemaImagingVTK ciVTK;
//       ciVTK.setDebugLevel( this->debugLevel_ );
//       ciVTK.setThreadNumber( this->threadNumber_ );
//       status = ciVTK.RenderImages(
//         outputImages,

//         inputObject,
//         inputGrid,

//         resolution,
//         this->CamProjectionMode,
//         this->CamAngle,
//         camFocus,
//         camNearFar,
//         camHeight,
//         camDir,
//         camUp,
//         this->GenerateScalarImages
//       );
  } else {
    //   status = this->RenderImagesViaEmbree(
    //     outputImages,

    //     inputObjectCells,
    //     inputObject,
    //     aInputGrid
    //   );
  }
  if(!status)
    return 0;

//   // print stats
//   this->printMsg(ttk::debug::Separator::L2);
//   this->printMsg(
//      "Complete (#i: " + std::to_string(inputGrid->GetNumberOfPoints())+"|#t:"+"-"+")",
//      1, globalTimer.getElapsedTime());
//   this->printMsg(ttk::debug::Separator::L1);

  return 1;
}

int ttkCinemaImaging::addFieldDataArray(
  vtkFieldData* fd,
  vtkDataArray* array,
  int tupelIdx,
  std::string name
){
  if(!array)
    return 0;

  auto newArray = vtkSmartPointer<vtkDataArray>::Take( array->NewInstance() );
  newArray->SetName( name.empty() ? array->GetName() : name.data() );
  newArray->SetNumberOfComponents( array->GetNumberOfComponents() );
  newArray->SetNumberOfTuples( 1 );
  newArray->SetTuple( 0, tupelIdx, array );

  fd->AddArray( newArray );

  return 1;
};

int ttkCinemaImaging::addAllFieldDataArrays(
  vtkPointSet* inputGrid,
  vtkImageData* image,
  int tupelIdx
){
  auto imageFD = image->GetFieldData();

  auto inputGridPD = inputGrid->GetPointData();
  for(int i=0; i<inputGridPD->GetNumberOfArrays(); i++){
    addFieldDataArray(
      imageFD,
      inputGridPD->GetArray(i),
      tupelIdx
    );
  }

  addFieldDataArray(
    imageFD,
    inputGrid->GetPoints()->GetData(),
    tupelIdx,
    "CamPos"
  );

  return 1;
};

int ttkCinemaImaging::computeDirFromFocus(
  vtkPointSet* inputGrid
){
  auto pos = static_cast<float*>(ttkUtils::GetVoidPointer(inputGrid->GetPoints()));
  auto focus = static_cast<double*>(ttkUtils::GetVoidPointer(inputGrid->GetPointData()->GetArray("CamFocus")));

  int nTuples = inputGrid->GetNumberOfPoints();

  auto newArray = vtkSmartPointer<vtkDoubleArray>::New();
  newArray->SetName("CamDir");
  newArray->SetNumberOfComponents( 3 );
  newArray->SetNumberOfTuples( nTuples );

  auto dir = static_cast<double*>(ttkUtils::GetVoidPointer(newArray));

  for(int i=0,j=nTuples*3; i<j; i++)
    dir[i] = focus[i]-pos[i];

  inputGrid->GetPointData()->AddArray(newArray);

  return 1;
};

int ttkCinemaImaging::ensureGridData(
  vtkPointData* fd,
  std::string name,
  int nTuples,
  const std::vector<double>& defaultValues
){
  auto array = vtkDoubleArray::SafeDownCast(fd->GetArray(name.data()));

  if(!array){
    int nComponents = defaultValues.size();

    auto newArray = vtkSmartPointer<vtkDoubleArray>::New();
    newArray->SetName(name.data());
    newArray->SetNumberOfComponents( nComponents );
    newArray->SetNumberOfTuples( nTuples );

    auto data = static_cast<double*>(ttkUtils::GetVoidPointer(newArray));
    for(int i=0,q=0; i<nTuples; i++)
      for(int j=0; j<nComponents; j++)
        data[q++] = defaultValues[j];

    fd->AddArray(newArray);
  }

  return 1;
};