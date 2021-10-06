#include <ttkCinemaImaging.h>

#include <vtkInformation.h>

#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>

#include <vtkCellData.h>
#include <vtkPointData.h>

#include <ttkUtils.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>

#include <ttkCinemaImagingEmbree.h>
#include <ttkCinemaImagingNative.h>
#include <ttkCinemaImagingVTK.h>

#include <ttkCinemaDarkroomCamera.h>

#include <CinemaImaging.h>

vtkStandardNewMacro(ttkCinemaImaging);

ttkCinemaImaging::ttkCinemaImaging() {
  this->setDebugMsgPrefix("CinemaImaging");
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
};

ttkCinemaImaging::~ttkCinemaImaging(){};

int ttkCinemaImaging::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Remove(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE());
    info->Append(
      vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    info->Append(
      vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
    info->Append(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
  } else if(port == 1)
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  else
    return 0;
  return 1;
};

int ttkCinemaImaging::FillOutputPortInformation(int port,
                                                vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
};

int ttkCinemaImaging::RequestData(vtkInformation *request,
                                  vtkInformationVector **inputVector,
                                  vtkInformationVector *outputVector) {

  auto object = vtkDataObject::GetData(inputVector[0]);
  auto cameras = vtkPointSet::GetData(inputVector[1]);
  if(!object || !cameras)
    return !this->printErr("Unsupported input object types.");

  ttkCinemaDarkroomCamera::Calibration calibration(cameras);
  if(calibration.nCameras < 0)
    return !this->printErr("Invalid Camera Calibrations.");

  if(calibration.nCameras == 0)
    return this->printWrn("Empty Camera Set.");

  auto objectAsMB = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  if(object->IsA("vtkMultiBlockDataSet"))
    objectAsMB->ShallowCopy(object);
  else
    objectAsMB->SetBlock(0, object);
  const size_t nObjects = objectAsMB->GetNumberOfBlocks();

  auto outputAsMB = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  for(size_t b = 0; b < nObjects; b++) {
    auto images = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    auto objectAsPS = vtkPointSet::SafeDownCast(objectAsMB->GetBlock(b));
    int status = this->RequestDataSingle(images, objectAsPS, cameras);
    if(!status)
      return 0;
    outputAsMB->SetBlock(b, images);
  }

  auto output = vtkMultiBlockDataSet::GetData(outputVector);
  if(object->IsA("vtkMultiBlockDataSet"))
    output->ShallowCopy(outputAsMB);
  else
    output->ShallowCopy(outputAsMB->GetBlock(0));

  return 1;
}

template <typename RT>
int renderDispatch(int debugLevel,
                   int threadNumber,
                   vtkMultiBlockDataSet *images,
                   vtkPointSet *object,
                   vtkPointSet *cameras) {
  RT renderer;
  renderer.setDebugLevel(debugLevel);
  renderer.setThreadNumber(threadNumber);
  return renderer.RenderVTKObject(images, object, cameras);
}

int ttkCinemaImaging::RequestDataSingle(vtkMultiBlockDataSet *images,
                                        vtkPointSet *object,
                                        vtkPointSet *cameras) {
  ttk::Timer globalTimer;

  auto cells = ttkCinemaImaging::GetCells(object);
  if(!cells)
    return !this->printErr("Unable to retrieve cell array.");

  size_t nTriangles = cells->GetNumberOfCells();
  // make sure that cells consists only of triangles
  {
    auto offsets = static_cast<vtkIdType *>(
      ttkUtils::GetVoidPointer(cells->GetOffsetsArray()));
    vtkIdType j = 0;
    for(size_t i = 0; i < nTriangles; i++, j += 3) {
      if(j != offsets[i]) {
        return !this->printErr(
          "Connectivity list has to contain only triangles.");
      }
    }
  }

  // perform rendering
  if(this->Backend == 0) {
    return renderDispatch<ttk::ttkCinemaImagingVTK>(
      this->debugLevel_, this->threadNumber_, images, object, cameras);
  } else if(this->Backend == 1) {
    return renderDispatch<ttk::ttkCinemaImagingEmbree>(
      this->debugLevel_, this->threadNumber_, images, object, cameras);
  } else {
    return renderDispatch<ttk::ttkCinemaImagingNative>(
      this->debugLevel_, this->threadNumber_, images, object, cameras);
  }
}

vtkCellArray *ttkCinemaImaging::GetCells(vtkPointSet *pointSet) {
  switch(pointSet->GetDataObjectType()) {
    case VTK_UNSTRUCTURED_GRID:
      return static_cast<vtkUnstructuredGrid *>(pointSet)->GetCells();
    case VTK_POLY_DATA:
      return static_cast<vtkPolyData *>(pointSet)->GetPolys();
  }

  return nullptr;
};

int ttkCinemaImaging::AddFieldDataArray(vtkFieldData *fd,
                                        vtkDataArray *array,
                                        int tupelIdx,
                                        const std::string &name) {
  if(!array)
    return 0;

  size_t nComponents = array->GetNumberOfComponents();

  auto newArray = vtkSmartPointer<vtkDoubleArray>::New();
  newArray->SetName(name.empty() ? array->GetName() : name.data());
  newArray->SetNumberOfComponents(nComponents);
  newArray->SetNumberOfTuples(1);

  if(newArray->GetDataType() == array->GetDataType()) {
    newArray->SetTuple(0, tupelIdx, array);
  } else {
    for(size_t i = 0; i < nComponents; i++)
      newArray->SetValue(
        i, array->GetVariantValue(tupelIdx * nComponents + i).ToDouble());
  }

  fd->AddArray(newArray);

  return 1;
};

int ttkCinemaImaging::AddAllFieldDataArrays(vtkPointSet *inputGrid,
                                            vtkImageData *image,
                                            int tupelIdx) {
  auto imageFD = image->GetFieldData();

  auto inputGridPD = inputGrid->GetPointData();
  for(int i = 0; i < inputGridPD->GetNumberOfArrays(); i++) {
    ttkCinemaImaging::AddFieldDataArray(
      imageFD, inputGridPD->GetArray(i), tupelIdx);
  }

  ttkCinemaImaging::AddFieldDataArray(
    imageFD, inputGrid->GetPoints()->GetData(), tupelIdx, "CamPosition");

  return 1;
};

int ttkCinemaImaging::ComputeDirFromFocalPoint(vtkPointSet *inputGrid) {
  auto pos
    = static_cast<float *>(ttkUtils::GetVoidPointer(inputGrid->GetPoints()));
  auto focal = static_cast<double *>(ttkUtils::GetVoidPointer(
    inputGrid->GetPointData()->GetArray("CamFocalPoint")));

  int nTuples = inputGrid->GetNumberOfPoints();

  auto newArray = vtkSmartPointer<vtkDoubleArray>::New();
  newArray->SetName("CamDirection");
  newArray->SetNumberOfComponents(3);
  newArray->SetNumberOfTuples(nTuples);

  auto dir = static_cast<double *>(ttkUtils::GetVoidPointer(newArray));

  for(int i = 0, j = nTuples * 3; i < j; i++)
    dir[i] = focal[i] - pos[i];

  // normalize
  for(int i = 0; i < nTuples; i++) {
    double *d = &dir[i * 3];
    const double norm = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    d[0] /= norm;
    d[1] /= norm;
    d[2] /= norm;
  }

  inputGrid->GetPointData()->AddArray(newArray);

  return 1;
};

int ttkCinemaImaging::Normalize(vtkDataArray *depthArray,
                                const double nearFar[2]) {
  if(!depthArray->IsA("vtkFloatArray")
     || depthArray->GetNumberOfComponents() != 1)
    return 0;

  if(nearFar[0] == 0.0 && nearFar[1] == 0.0)
    return 1;

  const size_t nValues = depthArray->GetNumberOfTuples();
  auto data = static_cast<float *>(ttkUtils::GetVoidPointer(depthArray));

  const float near = (float)nearFar[0];
  const float far = (float)nearFar[1];
  const float delta = far - near;

  for(size_t i = 0; i < nValues; i++) {
    if(std::isnan(data[i])) {
      data[i] = 1.0;
    } else {
      data[i] = (data[i] - near) / delta;
      if(data[i] > 1.0)
        data[i] = 1.0;
      else if(data[i] < 0.0)
        data[i] = 0.0;
    }
  }

  return 1;
};

int ttkCinemaImaging::MapPointAndCellData(
  vtkImageData *outputImage,

  vtkPointSet *inputObject,
  const ttk::CinemaImaging *renderer,
  const unsigned int *primitiveIdArray,
  const float *barycentricCoordinates,
  const vtkIdType *inputObjectConnectivityList) {

  auto inputObjectPD = inputObject->GetPointData();
  auto inputObjectCD = inputObject->GetCellData();
  auto outputImagePD = outputImage->GetPointData();
  int dim[3];
  outputImage->GetDimensions(dim);
  size_t nPixels = dim[0] * dim[1];

  const size_t nInputObjectPDArrays = inputObjectPD->GetNumberOfArrays();
  const size_t nInputObjectCDArrays = inputObjectCD->GetNumberOfArrays();

  int status = 0;

  // Map Point Data
  for(size_t j = 0; j < nInputObjectPDArrays; j++) {
    auto inputArray = inputObjectPD->GetArray(j);
    auto outputArray = vtkSmartPointer<vtkFloatArray>::New();
    outputArray->SetName(inputArray->GetName());
    outputArray->SetNumberOfComponents(inputArray->GetNumberOfComponents());
    outputArray->SetNumberOfTuples(nPixels);

    outputImagePD->AddArray(outputArray);

    switch(inputArray->GetDataType()) {
      vtkTemplateMacro(status = renderer->interpolateArray(
                         ttkUtils::GetPointer<float>(outputArray),

                         primitiveIdArray, barycentricCoordinates,
                         inputObjectConnectivityList,

                         ttkUtils::GetPointer<const VTK_TT>(inputArray),
                         nPixels, inputArray->GetNumberOfComponents()));
    }

    if(!status)
      return 0;
  }

  // Map Cell Data
  for(size_t j = 0; j < nInputObjectCDArrays; j++) {
    auto inputArray = inputObjectCD->GetArray(j);
    auto outputArray = vtkSmartPointer<vtkFloatArray>::New();
    outputArray->SetName(inputArray->GetName());
    outputArray->SetNumberOfComponents(inputArray->GetNumberOfComponents());
    outputArray->SetNumberOfTuples(nPixels);

    outputImagePD->AddArray(outputArray);

    switch(inputArray->GetDataType()) {
      vtkTemplateMacro(status = renderer->lookupArray(
                         ttkUtils::GetPointer<float>(outputArray),

                         primitiveIdArray,
                         ttkUtils::GetPointer<const VTK_TT>(inputArray),
                         nPixels, inputArray->GetNumberOfComponents()));
    }

    if(!status)
      return 0;
  }

  return 1;
};
