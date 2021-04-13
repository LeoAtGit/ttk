#include <ttkKernelDensityEstimation.h>


#include <vtkInformation.h>
#include <vtkFloatArray.h>
#include <vtkDataSet.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkInformationVector.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkPointSet.h>
#include <vtkImageData.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkKernelDensityEstimation);

ttkKernelDensityEstimation::ttkKernelDensityEstimation() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkKernelDensityEstimation::~ttkKernelDensityEstimation() {
}

int ttkKernelDensityEstimation::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0)
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPointSet");
  else
    return 0;

  return 1;
}

int ttkKernelDensityEstimation::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
  else
    return 0;

  return 1;
}

int ttkKernelDensityEstimation::RequestInformation(vtkInformation *,
                                   vtkInformationVector **,
                                   vtkInformationVector * outputVector) {

  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // Extent
  int wholeExtent[6]
    = {0, (int)this->Resolution[0]-1,
       0, (int)this->Resolution[1]-1,
       0, (int)this->Resolution[2]-1};
  outInfo->Set(
      vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), wholeExtent, 6);

  return 1;
}

int ttkKernelDensityEstimation::RequestData(vtkInformation *request,
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {

  auto inputPointSet = vtkPointSet::GetData(inputVector[0]);
  const size_t nPoints = inputPointSet->GetNumberOfPoints();

  auto categoryIndex = inputPointSet->GetPointData()->GetArray("CategoryIndex");
  double range[2];
  categoryIndex->GetRange(range);
  size_t nTypes = range[1] + 1; // +1 inc range bounds

  auto output = vtkImageData::GetData(outputVector);
  output->SetDimensions(
    this->Resolution[0],
    this->Resolution[1],
    this->Resolution[2]
  );
  output->SetOrigin(
    this->ImageBounds[0],
    this->ImageBounds[2],
    this->ImageBounds[4]
  );
  output->SetSpacing(
    this->Resolution[0]>1 ? (this->ImageBounds[1]-this->ImageBounds[0])/(this->Resolution[0]-1) : 0,
    this->Resolution[1]>1 ? (this->ImageBounds[3]-this->ImageBounds[2])/(this->Resolution[1]-1) : 0,
    this->Resolution[2]>1 ? (this->ImageBounds[5]-this->ImageBounds[4])/(this->Resolution[2]-1) : 0
  );

  output->AllocateScalars( VTK_FLOAT, 1 );
  // output->AllocateScalars( VTK_FLOAT, nTypes );

  auto kdeArray = output->GetPointData()->GetArray(0);
  kdeArray->SetName("KDE");
  auto kdeArrayData = ttkUtils::GetPointer<float>(kdeArray);

  auto nPixels = kdeArray->GetNumberOfTuples();

  auto kdeCatArray = vtkSmartPointer<vtkFloatArray>::New();
  kdeCatArray->SetName("KDE_ByType");
  kdeCatArray->SetNumberOfComponents(nTypes);
  kdeCatArray->SetNumberOfTuples(nPixels);
  output->GetPointData()->AddArray(kdeCatArray);
  auto kdeCatArrayData = ttkUtils::GetPointer<float>(kdeCatArray);

  auto countArray = vtkSmartPointer<vtkIntArray>::New();
  countArray->SetName("Count");
  countArray->SetNumberOfComponents( nTypes );
  countArray->SetNumberOfTuples( nPixels );
  output->GetPointData()->AddArray(countArray);
  auto countArrayData = ttkUtils::GetPointer<int>(countArray);

  int status = 0;

  status = this->computeCounts(
    countArrayData,

    ttkUtils::GetPointer<float>(inputPointSet->GetPoints()->GetData()),
    ttkUtils::GetPointer<unsigned char>(categoryIndex),
    nTypes,
    this->ImageBounds,
    this->Resolution,
    nPoints,
    nPixels
  );
  if(!status)
    return 0;

  switch(this->Kernel){
    case 0: {
      status = this->computeKDE2<KernelDensityEstimation::Gaussian>(
        kdeCatArrayData,
        countArrayData,
        nTypes,
        this->Bandwidth,
        this->ImageBounds,
        this->Resolution
      );
      break;
    }
    case 1: {
      status = this->computeKDE2<KernelDensityEstimation::Linear>(
        kdeCatArrayData,
        countArrayData,
        nTypes,
        this->Bandwidth,
        this->ImageBounds,
        this->Resolution
      );
      break;
    }
    case 2: {
      status = this->computeKDE2<KernelDensityEstimation::Epanechnikov>(
        kdeCatArrayData,
        countArrayData,
        nTypes,
        this->Bandwidth,
        this->ImageBounds,
        this->Resolution
      );
      break;
    }
  }

  status = this->computeTotalKDE(
    kdeArrayData,
    kdeCatArrayData,
    nPixels,
    nTypes
  );

  // On error cancel filter execution
  if(status == 0)
    return 0;

  return 1;
}
