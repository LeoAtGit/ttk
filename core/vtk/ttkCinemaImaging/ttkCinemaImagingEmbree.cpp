#include <ttkCinemaImagingEmbree.h>

#include <ttkCinemaImaging.h>
#include <ttkUtils.h>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>

#include <vtkFloatArray.h>
#include <vtkUnsignedIntArray.h>

#include <ttkCinemaDarkroomCamera.h>

ttk::ttkCinemaImagingEmbree::ttkCinemaImagingEmbree() {
  this->setDebugMsgPrefix("CinemaImaging(Embree)");
};
ttk::ttkCinemaImagingEmbree::~ttkCinemaImagingEmbree() {
}

int ttk::ttkCinemaImagingEmbree::RenderVTKObject(vtkMultiBlockDataSet *images,

                                                 vtkPointSet *object,
                                                 vtkPointSet *cameras) const {
  int status = 0;

#if TTK_ENABLE_EMBREE

  // initialize device
  RTCDevice device;
  status = this->initializeDevice(device);
  if(!status)
    return 0;

  // get cells
  auto cells = ttkCinemaImaging::GetCells(object);
  if(!cells)
    return !this->printErr("Unable to retrieve cell array.");

  auto connectivityList
    = ttkUtils::GetPointer<vtkIdType>(cells->GetConnectivityArray());

  // initilize scene
  RTCScene scene;
  status = this->initializeScene<vtkIdType>(
    scene,

    device, object->GetNumberOfPoints(),
    ttkUtils::GetPointer<float>(object->GetPoints()->GetData()),
    cells->GetNumberOfCells(), connectivityList);
  if(!status)
    return 0;

  // ---------------------------------------------------------------------------
  // Prepare Field Data for Depth Values
  // ---------------------------------------------------------------------------

  // iterate over sampling locations
  this->printMsg(ttk::debug::Separator::L2);

  // fetch camera calibrations
  ttkCinemaDarkroomCamera::Calibration calibration(cameras);

  for(int c = 0; c < calibration.nCameras; c++) {

    auto resolution = calibration.Resolution.Get(c);

    // Initialize Output
    auto image = vtkSmartPointer<vtkImageData>::New();
    image->SetDimensions(resolution[0], resolution[1], 1);
    image->SetSpacing(1, 1, 1);
    image->SetOrigin(0, 0, 0);
    image->AllocateScalars(VTK_FLOAT, 1);

    size_t nPixels = resolution[0] * resolution[1];
    auto imagePD = image->GetPointData();

    auto depthBuffer = imagePD->GetArray(0);
    depthBuffer->SetName("Depth");

    auto primitiveIdArray = vtkSmartPointer<vtkUnsignedIntArray>::New();
    primitiveIdArray->SetName("PrimitiveId");
    primitiveIdArray->SetNumberOfComponents(1);
    primitiveIdArray->SetNumberOfTuples(nPixels);
    auto primitiveIdArrayData
      = ttkUtils::GetPointer<unsigned int>(primitiveIdArray);
    imagePD->AddArray(primitiveIdArray);

    auto barycentricCoordinates = vtkSmartPointer<vtkFloatArray>::New();
    barycentricCoordinates->SetName("BarycentricCoordinates");
    barycentricCoordinates->SetNumberOfComponents(2);
    barycentricCoordinates->SetNumberOfTuples(nPixels);
    auto barycentricCoordinatesData
      = ttkUtils::GetPointer<float>(barycentricCoordinates);
    imagePD->AddArray(barycentricCoordinates);

    // Render Object
    status = this->renderImage(
      ttkUtils::GetPointer<float>(depthBuffer), primitiveIdArrayData,
      barycentricCoordinatesData,

      scene, resolution.data(), calibration.Position.Get(c).data(),
      calibration.Direction.Get(c).data(), calibration.Up.Get(c).data(),
      calibration.Scale.Array ? calibration.Scale.Get(c)[0]
                              : calibration.Angle.Get(c)[0],
      calibration.Scale.Array);
    if(!status)
      return 0;

    ttkCinemaImaging::Normalize(depthBuffer, calibration.NearFar.Get(c).data());

    status = ttkCinemaImaging::MapPointAndCellData(
      image,

      object, this, primitiveIdArrayData, barycentricCoordinatesData,
      connectivityList);
    if(!status)
      return !this->printErr("Unable to Map Point and Cell Data.");

    ttkCinemaImaging::AddAllFieldDataArrays(cameras, image, c);

    images->SetBlock(c, image);
  }
  this->printMsg(ttk::debug::Separator::L2);

  this->deallocateScene(device, scene);

#else
  this->printErr("TTK was build without EMBREE backend.");
#endif

  return status;
};
