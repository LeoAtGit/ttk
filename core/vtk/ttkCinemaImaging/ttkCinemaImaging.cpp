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

#include <ttkCinemaImagingVTK.h>
#include <ttkCinemaImagingEmbree.h>
#include <ttkCinemaImagingNative.h>

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

  if(!inputObject || !inputGrid)
    return 0;

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
      defaultNearFar[0] = std::max(defaultNearFar[0], 0.01);
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
  ttkCinemaImaging::EnsureGridData( aInputGridPD, "Resolution", n, {(double)this->Resolution[0],(double)this->Resolution[1]} );
  ttkCinemaImaging::EnsureGridData( aInputGridPD, "ProjectionMode", n, {(double)this->CamProjectionMode} );
  if(aInputGridPD->HasArray("CamDirection"))
    ttkCinemaImaging::EnsureGridData( aInputGridPD, "CamDirection", n, {0,0,-1} );
  else {
    ttkCinemaImaging::EnsureGridData( aInputGridPD, "CamFocus", n, defaultFocus );
    ttkCinemaImaging::ComputeDirFromFocus(aInputGrid);
  }
  ttkCinemaImaging::EnsureGridData( aInputGridPD, "CamNearFar", n, defaultNearFar );
  ttkCinemaImaging::EnsureGridData( aInputGridPD, "CamHeight", n, {defaultHeight} );
  ttkCinemaImaging::EnsureGridData( aInputGridPD, "CamAngle", n, {defaultAngle} );
  ttkCinemaImaging::EnsureGridData( aInputGridPD, "CamUp", n, {0,1,0} );

  // print parameters
  {
    this->printMsg({
      {"Backend", this->Backend==0 ? std::string("VTK_OPENGL") : this->Backend==1 ? std::string("EMBREE") : std::string("NATIVE")},
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

  auto cells = ttkCinemaImaging::GetCells(inputObject);
  if(!cells)
    return 0;

  size_t nTriangles = cells->GetNumberOfCells();
  // make sure that cells consists only of triangles
  {
      auto offsets = static_cast<vtkIdType*>(
          ttkUtils::GetVoidPointer(
              cells->GetOffsetsArray()
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
    ttk::ttkCinemaImagingVTK renderer;
    renderer.setDebugLevel(this->debugLevel_);
    renderer.setThreadNumber(this->threadNumber_);
    status = renderer.RenderVTKObject(
        outputImages,
        inputObject,
        aInputGrid
      );
  } else if(this->Backend==1) {
    ttk::ttkCinemaImagingEmbree renderer;
    renderer.setDebugLevel(this->debugLevel_);
    renderer.setThreadNumber(this->threadNumber_);
    status = renderer.RenderVTKObject(
        outputImages,
        inputObject,
        aInputGrid
      );
  } else {
    ttk::ttkCinemaImagingNative renderer;
    renderer.setDebugLevel(this->debugLevel_);
    renderer.setThreadNumber(this->threadNumber_);
    status = renderer.RenderVTKObject(
        outputImages,
        inputObject,
        aInputGrid
      );
  }

  if(!status)
    return 0;

  return 1;
}

vtkCellArray* ttkCinemaImaging::GetCells(vtkPointSet* pointSet){
  switch(pointSet->GetDataObjectType()) {
    case VTK_UNSTRUCTURED_GRID:
      return static_cast<vtkUnstructuredGrid*>(pointSet)->GetCells();
    case VTK_POLY_DATA:
      return static_cast<vtkPolyData*>(pointSet)->GetPolys();
  }
  return nullptr;
};

int ttkCinemaImaging::AddFieldDataArray(
  vtkFieldData* fd,
  vtkDataArray* array,
  int tupelIdx,
  std::string name
){
  if(!array)
    return 0;

  size_t nComponents = array->GetNumberOfComponents();

  auto newArray = vtkSmartPointer<vtkDoubleArray>::New();
  newArray->SetName( name.empty() ? array->GetName() : name.data() );
  newArray->SetNumberOfComponents( nComponents );
  newArray->SetNumberOfTuples( 1 );

  if(newArray->GetDataType()==array->GetDataType()){
    newArray->SetTuple( 0, tupelIdx, array );
  } else {
    for(size_t i=0; i<nComponents; i++)
        newArray->SetValue(i, array->GetVariantValue(tupelIdx*nComponents+i).ToDouble());
  }

  fd->AddArray( newArray );

  return 1;
};

int ttkCinemaImaging::AddAllFieldDataArrays(
  vtkPointSet* inputGrid,
  vtkImageData* image,
  int tupelIdx
){
  auto imageFD = image->GetFieldData();

  auto inputGridPD = inputGrid->GetPointData();
  for(int i=0; i<inputGridPD->GetNumberOfArrays(); i++){
    ttkCinemaImaging::AddFieldDataArray(
      imageFD,
      inputGridPD->GetArray(i),
      tupelIdx
    );
  }

  ttkCinemaImaging::AddFieldDataArray(
    imageFD,
    inputGrid->GetPoints()->GetData(),
    tupelIdx,
    "CamPosition"
  );

  return 1;
};

int ttkCinemaImaging::ComputeDirFromFocus(
  vtkPointSet* inputGrid
){
  auto pos = static_cast<float*>(ttkUtils::GetVoidPointer(inputGrid->GetPoints()));
  auto focus = static_cast<double*>(ttkUtils::GetVoidPointer(inputGrid->GetPointData()->GetArray("CamFocus")));

  int nTuples = inputGrid->GetNumberOfPoints();

  auto newArray = vtkSmartPointer<vtkDoubleArray>::New();
  newArray->SetName("CamDirection");
  newArray->SetNumberOfComponents( 3 );
  newArray->SetNumberOfTuples( nTuples );

  auto dir = static_cast<double*>(ttkUtils::GetVoidPointer(newArray));

  for(int i=0,j=nTuples*3; i<j; i++)
    dir[i] = focus[i]-pos[i];

  inputGrid->GetPointData()->AddArray(newArray);

  return 1;
};

int ttkCinemaImaging::EnsureGridData(
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
