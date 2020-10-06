#include <ttkCinemaImagingVTK.h>

#include <vtkSmartPointer.h>
#include <vtkDataSetSurfaceFilter.h>

#include <vtkPointSet.h>
#include <vtkMultiBlockDataSet.h>

#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkRenderWindow.h>
#include <vtkWindowToImageFilter.h>
#include <vtkCameraPass.h>
#include <vtkSequencePass.h>
#include <vtkValuePass.h>
#include <vtkRenderPassCollection.h>
#include <vtkRenderer.h>
#include <vtkOpenGLRenderer.h>

#include <vtkSignedCharArray.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>

#include <vtkPointData.h>
#include <vtkCellData.h>

ttk::ttkCinemaImagingVTK::ttkCinemaImagingVTK(){
    this->setDebugMsgPrefix("CinemaImaging(VTK)");
};
ttk::ttkCinemaImagingVTK::~ttkCinemaImagingVTK(){
}

int ttk::ttkCinemaImagingVTK::setupRenderer(
    vtkRenderer *renderer,
    vtkPointSet *object,
    vtkCamera *camera
) const {
  auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputDataObject(object);

  auto actor = vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);

  renderer->SetBackground(0.0, 0.0, 0.0);
  renderer->GradientBackgroundOff();
  renderer->AddActor(actor);
  renderer->SetActiveCamera(camera);

  return 1;
}

int ttk::ttkCinemaImagingVTK::setupWindow(
    vtkRenderWindow *window,
    vtkRenderer *renderer,
    const double resolution[2]
) const {
  window->SetSize(resolution[0], resolution[1]);
  window->SetMultiSamples(0); // Disable AA
  window->OffScreenRenderingOn();
  window->AddRenderer(renderer);

  return 1;
};

int ttk::ttkCinemaImagingVTK::addValuePass(
    vtkPointSet* object,
    int fieldType,
    vtkRenderPassCollection *valuePassCollection,
    std::vector<std::string> &valuePassNames
) const {
  auto fd = static_cast<vtkFieldData*>(
      fieldType == 0
        ? (vtkFieldData *)object->GetPointData()
        : (vtkFieldData *)object->GetCellData()
  );

  for(size_t i=0, j=fd->GetNumberOfArrays(); i<j; i++){
    auto array = fd->GetArray(i);
    if(!array)
      continue;

    std::string name(array->GetName());

    double minmax[2];
    array->GetRange(minmax);

    size_t nComponents = array->GetNumberOfComponents();
    for(size_t c = 0; c < nComponents; c++){
      auto valuePass = vtkSmartPointer<vtkValuePass>::New();
      valuePass->SetRenderingMode(vtkValuePass::FLOATING_POINT);
      valuePass->SetInputArrayToProcess(fieldType == 0
                                          ? VTK_SCALAR_MODE_USE_POINT_FIELD_DATA
                                          : VTK_SCALAR_MODE_USE_CELL_FIELD_DATA,
                                        name.data());
      valuePass->SetInputComponentToProcess(c);

      valuePassCollection->AddItem(valuePass);
      valuePassNames.push_back(nComponents == 1 ? name : name + "_"
                                                           + std::to_string(c));
    }
  }

  return 1;
};

int ttk::ttkCinemaImagingVTK::addGobalArray(
    std::vector<vtkAbstractArray*>& globalArrays,
    std::string name,
    size_t nValues,
    const double *values
) const {
    auto array = vtkSmartPointer<vtkDoubleArray>::New();
    array->SetName(name.data());
    array->SetNumberOfComponents(nValues);
    array->SetNumberOfTuples(1);
    for(size_t i = 0; i < nValues; i++)
        array->SetValue(i, values[i]);
    globalArrays.push_back(array);

    return 1;
};

int ttk::ttkCinemaImagingVTK::RenderImages(
    vtkMultiBlockDataSet* outputImages,

    vtkPointSet* inputObject,
    vtkPointSet* inputGrid,

    const double defaultResolution[2],
    const int    projectionMode,
    const double defaultCamAngle,
    const double defaultCamFocus[3],
    const double defaultCamNearFar[2],
    const double defaultCamHeight,
    const double defaultCamDir[3],
    const double defaultCamUp[3],
    const bool   generateScalarImages
) const {
  ttk::Timer timer;

  double resolution[2] = {defaultResolution[0],defaultResolution[1]};
  double camFocus[3] = {defaultCamFocus[0],defaultCamFocus[1],defaultCamFocus[2]};
  double camNearFar[2] = {defaultCamNearFar[0],defaultCamNearFar[1]};
  double camHeight = defaultCamHeight;
  double camDir[3] = {defaultCamDir[0],defaultCamDir[1],defaultCamDir[2]};
  double camUp[3] = {defaultCamUp[0],defaultCamUp[1],defaultCamUp[2]};

  // Camera
  auto camera = vtkSmartPointer<vtkCamera>::New();
  camera->SetClippingRange(camNearFar);
  camera->SetFocalPoint(camFocus);
  if(projectionMode == 0) {
    camera->SetParallelProjection(true);
    camera->SetParallelScale(
      camHeight * 0.5); // *0.5 to convert CamHeight to weird VTK convention
  } else {
    camera->SetParallelProjection(false);
    camera->SetViewAngle( defaultCamAngle );
  }

  // Initialize Depth Renderer and Components
  this->printMsg(
    "Initializing Rendering Pipeline",
    0, ttk::debug::LineMode::REPLACE
  );

  auto inputObjectAsPD = vtkSmartPointer<vtkPolyData>::New();
  if(inputObject->IsA("vtkPolyData")){
    inputObjectAsPD->ShallowCopy(inputObject);
  } else {
    this->printWrn("VTK backend can only handle vtkPolyData input.");
    this->printWrn(" -> Extracting input surface");
    auto surfaceFilter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    surfaceFilter->SetInputDataObject( inputObject );
    surfaceFilter->Update();
    inputObjectAsPD->ShallowCopy(surfaceFilter->GetOutput());
  }

  // Depth Pass Elements
  auto rendererDepth = vtkSmartPointer<vtkRenderer>::New();
  this->setupRenderer(rendererDepth, inputObjectAsPD, camera);
  auto windowDepth = vtkSmartPointer<vtkRenderWindow>::New();
  this->setupWindow(windowDepth, rendererDepth, resolution);
  auto windowDepthToImageFilter
    = vtkSmartPointer<vtkWindowToImageFilter>::New();
  windowDepthToImageFilter->SetInput(windowDepth);
  windowDepthToImageFilter->SetInputBufferTypeToZBuffer();

  // Value passes Elements
  size_t nValuePasses = 0;

  auto rendererScalars = vtkSmartPointer<vtkRenderer>::New();
  this->setupRenderer(rendererScalars, inputObjectAsPD, camera);

  auto windowScalars = vtkSmartPointer<vtkRenderWindow>::New();
  this->setupWindow(windowScalars, rendererScalars, resolution);

  auto valuePassCollection = vtkSmartPointer<vtkRenderPassCollection>::New();
  std::vector<std::string> valuePassNames;
  size_t firstValuePassIndex = 0;
  if( generateScalarImages ) {
    auto preventVTKBug = [](vtkPointSet *object) {
      auto pd = object->GetPointData();

      if(pd->GetNumberOfArrays()<1) {
        size_t nP = object->GetNumberOfPoints();

        auto fakeArray = vtkSmartPointer<vtkSignedCharArray>::New();
        fakeArray->SetName("Fake");
        fakeArray->SetNumberOfComponents(1);
        fakeArray->SetNumberOfTuples(nP);
        auto fakeArrayData = (signed char *)fakeArray->GetVoidPointer(0);
        for(size_t i = 0; i<nP; i++)
          fakeArrayData[i] = 0;
        pd->AddArray(fakeArray);
        return 1;
      }
      return 0;
    };

    if(preventVTKBug(inputObjectAsPD)) {
      firstValuePassIndex = 1;
    };

    this->addValuePass(inputObjectAsPD, 0, valuePassCollection, valuePassNames);
    this->addValuePass(inputObjectAsPD, 1, valuePassCollection, valuePassNames);
    nValuePasses = valuePassNames.size();

    auto sequence = vtkSmartPointer<vtkSequencePass>::New();
    sequence->SetPasses(valuePassCollection);

    auto cameraPass = vtkSmartPointer<vtkCameraPass>::New();
    cameraPass->SetDelegatePass(sequence);

    auto glRenderer = vtkOpenGLRenderer::SafeDownCast(rendererScalars);
    glRenderer->SetPass(cameraPass);

    // First pass to setup everything
    windowScalars->Render();
  }

  this->printMsg(
    "Initializing Rendering Pipeline",
    1, timer.getElapsedTime()
  );
  timer.reStart();

  // ---------------------------------------------------------------------------
  // Prepare Field Data for Depth Values
  // ---------------------------------------------------------------------------
  std::vector<vtkAbstractArray*> globalArrays;
  {
    this->addGobalArray(globalArrays, "CamHeight", 1, &camHeight);
    this->addGobalArray(globalArrays, "CamNearFar", 2, camNearFar);
    this->addGobalArray(globalArrays, "Resolution", 2, resolution);

    auto inputObjectFD = inputObjectAsPD->GetFieldData();
    for(size_t i = 0, j = inputObjectFD->GetNumberOfArrays(); i < j; i++) {
      auto array = inputObjectFD->GetAbstractArray(i);
      auto copy = vtkSmartPointer<vtkAbstractArray>::Take(array->NewInstance());
      copy->DeepCopy(array);
      globalArrays.push_back(array);
    }
  }
  size_t nGlobalArrays = globalArrays.size();

  // ---------------------------------------------------------------------------
  // Iterate over Locations
  // ---------------------------------------------------------------------------
  auto inputGridPD = inputGrid->GetPointData();
  size_t nInputGridPD = inputGridPD->GetNumberOfArrays();

  double *camFocusData = nullptr;
  double *camDirData = nullptr;
  double *camUpData = nullptr;

  {
    auto checkGridArray
      = [](vtkFieldData *fd, std::string name, std::string &gridFieldNames,
          double *&data, const ttkCinemaImagingVTK *caller) {
          if(fd->HasArray(name.data())) {
            auto field
              = vtkDoubleArray::SafeDownCast(fd->GetAbstractArray(name.data()));
            if(!field) {
              caller->printErr(
                "Sampling grid field '" + name
                + "' is not a vtkDoubleArray with three components.");
              return 0;
            }
            data = (double *)field->GetVoidPointer(0);
            gridFieldNames += " '" + name + "'";
          }
          return 1;
        };

    std::string gridFieldNames = "";
    if(!checkGridArray(
        inputGridPD, "CamFocus", gridFieldNames, camFocusData, this))
      return 0;
    if(!checkGridArray(
        inputGridPD, "CamDirection", gridFieldNames, camDirData, this))
      return 0;
    if(!checkGridArray(inputGridPD, "CamUp", gridFieldNames, camUpData, this))
      return 0;

    if(gridFieldNames.compare("") != 0)
      this->printMsg("Sampling grid has field(s): " + gridFieldNames);
  }

  // -------------------------------------------------------------------------
  // Render Images for all Camera Locations
  // -------------------------------------------------------------------------
  {
    timer.reStart();
    size_t n = inputGrid->GetNumberOfPoints();
    this->printMsg(
      "Rendering " + std::to_string(n) + " images with "
        + std::to_string(nValuePasses + 1) + " fields",
      0, ttk::debug::LineMode::REPLACE
    );

    double camPosition[3] = {0, 0, 0};
    auto readCameraData = [](double target[3], double *src, int index) {
      target[0] = src[index];
      target[1] = src[index + 1];
      target[2] = src[index + 2];
    };

    auto addCamFieldData
      = [](vtkFieldData *fd, std::string name, double *data) {
          auto array = vtkSmartPointer<vtkDoubleArray>::New();
          array->SetName(name.data());
          array->SetNumberOfComponents(3);
          array->SetNumberOfTuples(1);
          array->SetValue(0, data[0]);
          array->SetValue(1, data[1]);
          array->SetValue(2, data[2]);
          fd->AddArray(array);
        };

    for(size_t i = 0; i < n; i++) {
      // Set Camera Position
      inputGrid->GetPoint(i, camPosition);

      // Cam Up Fix
      if(camPosition[0] == 0 && camPosition[2] == 0) {
        camPosition[0] = 0.00000000001;
        camPosition[2] = 0.00000000001;
      }
      camera->SetPosition(camPosition);

      if(camUpData != nullptr) {
        readCameraData(camUp, camUpData, i * 3);
        vtkMath::Normalize(camUp);
      }
      camera->SetViewUp(camUp);

      if(camFocusData != nullptr) {
        readCameraData(camFocus, camFocusData, i * 3);
        camera->SetFocalPoint(camFocus);
      }
      if(camDirData != nullptr) {
        readCameraData(camDir, camDirData, i * 3);
        camFocus[0] = camPosition[0] + camDir[0];
        camFocus[1] = camPosition[1] + camDir[1];
        camFocus[2] = camPosition[2] + camDir[2];
        camera->SetFocalPoint(camFocus);
      } else {
        camDir[0] = camFocus[0] - camPosition[0];
        camDir[1] = camFocus[1] - camPosition[1];
        camDir[2] = camFocus[2] - camPosition[2];
        vtkMath::Normalize(camDir);
      }

      // Initialize Output Image
      auto outputImage = vtkSmartPointer<vtkImageData>::New();
      auto outputImagePD = outputImage->GetPointData();

      // Initialize as depth image
      {
        windowDepthToImageFilter->Modified();
        windowDepthToImageFilter->Update();
        outputImage->DeepCopy(windowDepthToImageFilter->GetOutput());
        outputImagePD->GetAbstractArray(0)->SetName("Depth");
      }

      // Render Scalar Images
      if(nValuePasses > firstValuePassIndex) {
        windowScalars->Render();

        for(size_t j = firstValuePassIndex; j < nValuePasses; j++) {
          auto valuePass = vtkValuePass::SafeDownCast(
            valuePassCollection->GetItemAsObject(j));
          auto newValueArray = vtkSmartPointer<vtkFloatArray>::New();
          newValueArray->DeepCopy(
            valuePass->GetFloatImageDataArray(rendererScalars));
          newValueArray->SetName(valuePassNames[j].data());
          outputImagePD->AddArray(newValueArray);
        }
      }

      // Add Field Data
      {
        auto outputImageFD = outputImage->GetFieldData();

        // Global Arrays
        for(size_t j = 0; j < nGlobalArrays; j++) {
          outputImageFD->AddArray(globalArrays[j]);
        }

        // Specific Arrays
        addCamFieldData(outputImageFD, "CamPosition", camPosition);
        addCamFieldData(outputImageFD, "CamDirection", camDir);
        addCamFieldData(outputImageFD, "CamUp", camera->GetViewUp());

        for(size_t j = 0; j < nInputGridPD; j++) {
          auto array = inputGridPD->GetAbstractArray(j);
          auto newArray
            = vtkSmartPointer<vtkAbstractArray>::Take(array->NewInstance());
          std::string name(array->GetName());
          if(outputImageFD->HasArray(name.data()))
            name += "FromGrid";
          newArray->SetName(name.data());
          newArray->SetNumberOfComponents(array->GetNumberOfComponents());
          newArray->SetNumberOfTuples(1);
          newArray->SetTuple(0, i, array);

          outputImageFD->AddArray(newArray);
        }
      }

      // Add Image to MultiBlock
      outputImages->SetBlock(i, outputImage);
    }

    this->printMsg("Rendering " + std::to_string(n) + " images with "
        + std::to_string(nValuePasses + 1) + " fields",
      1, timer.getElapsedTime()
    );
  }

  return 1;
};