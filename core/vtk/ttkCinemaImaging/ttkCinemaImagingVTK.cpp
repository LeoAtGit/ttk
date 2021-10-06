#include <ttkCinemaImagingVTK.h>

#include <ttkCinemaImaging.h>

#include <vtkDataSetSurfaceFilter.h>
#include <vtkSmartPointer.h>

#include <vtkMultiBlockDataSet.h>
#include <vtkPointSet.h>

#include <ttkCinemaDarkroomCamera.h>

#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCameraPass.h>
#include <vtkOpenGLRenderer.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderPassCollection.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkSequencePass.h>
#include <vtkValuePass.h>
#include <vtkWindowToImageFilter.h>

#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkSignedCharArray.h>

#include <ttkUtils.h>
#include <vtkCellData.h>
#include <vtkPointData.h>

ttk::ttkCinemaImagingVTK::ttkCinemaImagingVTK() {
  this->setDebugMsgPrefix("CinemaImaging(VTK)");
};
ttk::ttkCinemaImagingVTK::~ttkCinemaImagingVTK() {
}

int ttk::ttkCinemaImagingVTK::setupRenderer(vtkRenderer *renderer,
                                            vtkPointSet *object,
                                            vtkCamera *camera) const {
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

int ttk::ttkCinemaImagingVTK::setupWindow(vtkRenderWindow *window,
                                          vtkRenderer *renderer) const {
  window->SetMultiSamples(0); // Disable AA
  window->OffScreenRenderingOn();
  window->AddRenderer(renderer);

  return 1;
};

int ttk::ttkCinemaImagingVTK::addValuePass(
  vtkPointSet *object,
  int fieldType,
  vtkRenderPassCollection *valuePassCollection,
  std::vector<std::string> &valuePassNames) const {
  auto fd = fieldType == 0 ? static_cast<vtkFieldData *>(object->GetPointData())
                           : static_cast<vtkFieldData *>(object->GetCellData());

  for(size_t i = 0, j = fd->GetNumberOfArrays(); i < j; i++) {
    auto array = fd->GetArray(i);
    if(!array)
      continue;

    std::string name(array->GetName());

    double minmax[2];
    array->GetRange(minmax);

    size_t nComponents = array->GetNumberOfComponents();
    for(size_t c = 0; c < nComponents; c++) {
      auto valuePass = vtkSmartPointer<vtkValuePass>::New();
      valuePass->SetInputArrayToProcess(fieldType == 0
                                          ? VTK_SCALAR_MODE_USE_POINT_FIELD_DATA
                                          : VTK_SCALAR_MODE_USE_CELL_FIELD_DATA,
                                        name.data());
      valuePass->SetInputComponentToProcess(c);

      valuePassCollection->AddItem(valuePass);
      valuePassNames.push_back(
        nComponents == 1 ? name : name + "_" + std::to_string(c));
    }
  }

  return 1;
};

int ttk::ttkCinemaImagingVTK::RenderVTKObject(vtkMultiBlockDataSet *images,

                                              vtkPointSet *object,
                                              vtkPointSet *cameras) const {

  auto objectAsPD = vtkSmartPointer<vtkPolyData>::New();
  if(object->IsA("vtkPolyData")) {
    objectAsPD->ShallowCopy(object);
  } else {
    auto surfaceFilter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    surfaceFilter->SetInputDataObject(object);
    surfaceFilter->Update();
    objectAsPD->ShallowCopy(surfaceFilter->GetOutput());
  }

  // ---------------------------------------------------------------------------
  // Prepare Field Data for Depth Values
  // ---------------------------------------------------------------------------

  // iterate over sampling locations
  this->printMsg(ttk::debug::Separator::L2);

  // Camera
  auto camera = vtkSmartPointer<vtkCamera>::New();

  // Depth Pass Elements
  auto rendererDepth = vtkSmartPointer<vtkRenderer>::New();
  this->setupRenderer(rendererDepth, objectAsPD, camera);

  auto windowDepth = vtkSmartPointer<vtkRenderWindow>::New();
  this->setupWindow(windowDepth, rendererDepth);

  auto windowDepthToImageFilter
    = vtkSmartPointer<vtkWindowToImageFilter>::New();
  windowDepthToImageFilter->SetInput(windowDepth);
  windowDepthToImageFilter->SetInputBufferTypeToZBuffer();

  // Value passes Elements
  size_t nValuePasses = 0;

  auto rendererScalars = vtkSmartPointer<vtkRenderer>::New();
  this->setupRenderer(rendererScalars, objectAsPD, camera);

  auto windowScalars = vtkSmartPointer<vtkRenderWindow>::New();
  this->setupWindow(windowScalars, rendererScalars);

  auto valuePassCollection = vtkSmartPointer<vtkRenderPassCollection>::New();
  std::vector<std::string> valuePassNames;
  size_t firstValuePassIndex = 0;

  auto preventVTKBug = [](vtkPointSet *o) {
    auto pd = o->GetPointData();
    auto cd = o->GetCellData();

    if(pd->GetNumberOfArrays() < 1 && cd->GetNumberOfArrays() > 0) {
      size_t nP = o->GetNumberOfPoints();

      auto fakeArray = vtkSmartPointer<vtkSignedCharArray>::New();
      fakeArray->SetName("Fake");
      fakeArray->SetNumberOfComponents(1);
      fakeArray->SetNumberOfTuples(nP);
      auto fakeArrayData = (signed char *)fakeArray->GetVoidPointer(0);
      for(size_t i = 0; i < nP; i++)
        fakeArrayData[i] = 0;
      pd->AddArray(fakeArray);
      return 1;
    }
    return 0;
  };

  if(preventVTKBug(objectAsPD)) {
    firstValuePassIndex = 1;
  };

  this->addValuePass(objectAsPD, 0, valuePassCollection, valuePassNames);
  this->addValuePass(objectAsPD, 1, valuePassCollection, valuePassNames);
  nValuePasses = valuePassNames.size();

  auto sequence = vtkSmartPointer<vtkSequencePass>::New();
  sequence->SetPasses(valuePassCollection);

  auto cameraPass = vtkSmartPointer<vtkCameraPass>::New();
  cameraPass->SetDelegatePass(sequence);

  auto glRenderer = vtkOpenGLRenderer::SafeDownCast(rendererScalars);
  glRenderer->SetPass(cameraPass);

  // First pass to setup everything
  windowScalars->Render();

  // fetch camera calibrations
  ttkCinemaDarkroomCamera::Calibration calibration(cameras);

  for(int c = 0; c < calibration.nCameras; c++) {
    ttk::Timer timer;
    const bool orthographic = calibration.Scale.Array != nullptr;
    const auto &resolution = calibration.Resolution.Get(c);

    int resX = resolution[0];
    int resY = resolution[1];

    const std::string msg{
      "Rendering Image (" + std::string(orthographic ? "O" : "P") + "|"
      + std::to_string(resX) + "x" + std::to_string(resY) + ")"};
    this->printMsg(msg, 0, 0, ttk::debug::LineMode::REPLACE);

    if(orthographic) {
      camera->SetParallelProjection(true);
      camera->SetParallelScale(
        calibration.Scale.Get(c)[0]
        * 0.5); // *0.5 to convert CamScale to weird VTK convention
    } else {
      camera->SetParallelProjection(false);
      camera->SetViewAngle(calibration.Angle.Get(c)[0]);
    }

    auto pos = calibration.Position.Get(c);
    auto dir = calibration.Direction.Get(c);

    camera->SetPosition(pos.data());
    camera->SetViewUp(calibration.Up.Get(c).data());
    camera->SetFocalPoint(pos[0] + dir[0], pos[1] + dir[1], pos[2] + dir[2]);
    camera->SetClippingRange(calibration.NearFar.Get(c).data());

    // perform depth pass
    windowDepth->SetSize(resX, resY);
    windowDepthToImageFilter->Modified();
    windowDepthToImageFilter->Update();

    // Initialize output as depth image
    auto outputImage = vtkSmartPointer<vtkImageData>::New();
    outputImage->DeepCopy(windowDepthToImageFilter->GetOutput());

    auto outputImagePD = outputImage->GetPointData();
    outputImagePD->GetAbstractArray(0)->SetName("Depth");

    ttkCinemaImaging::AddAllFieldDataArrays(cameras, outputImage, c);

    // Render Scalar Images
    if(nValuePasses > firstValuePassIndex) {

      if(windowScalars->GetSize()[0] != resX
         || windowScalars->GetSize()[1] != resY) {
        windowScalars->SetSize(resX, resY);
        windowScalars->Render(); // render twice to prevent VTK Bug
      }
      windowScalars->Render();

      for(size_t p = firstValuePassIndex; p < nValuePasses; p++) {
        auto valuePass
          = vtkValuePass::SafeDownCast(valuePassCollection->GetItemAsObject(p));
        auto newValueArray = vtkSmartPointer<vtkFloatArray>::New();
        newValueArray->DeepCopy(
          valuePass->GetFloatImageDataArray(rendererScalars));
        newValueArray->SetName(valuePassNames[p].data());
        outputImagePD->AddArray(newValueArray);
      }
    }

    // Add to output collection
    images->SetBlock(c, outputImage);

    this->printMsg(msg, 1, timer.getElapsedTime());
  }
  this->printMsg(ttk::debug::Separator::L2);

  return 1;
};
