#include <ttkPerlinVectorField.h>

#include <vtkInformation.h>
#include <vtkInformationVector.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkDoubleArray.h>


#include <ttkMacros.h>
#include <ttkUtils.h>

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkPerlinVectorField);

ttkPerlinVectorField::ttkPerlinVectorField() {
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
}

ttkPerlinVectorField::~ttkPerlinVectorField() {
}

int ttkPerlinVectorField::FillInputPortInformation(int port, vtkInformation *info) {
  return 0;
}

int ttkPerlinVectorField::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
    return 1;
  }
  return 0;
}

int ttkPerlinVectorField::RequestInformation (
  vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector *outputVector
  ) {
    // For vtkImageData output we have to set extent, spacing and origin already in request information
    vtkInformation* outInfo = outputVector->GetInformationObject(0);
    //int dimX, dimY;
    //dimX = dimY = 5;

    int extent[6] = {0, ImageDimension[0] - 1, 0, ImageDimension[1] - 1, 0, 0};
    double spacing[3] = {1.0, 1.0, 1.0};
    double origin[3] = {0.0, 0.0, 0.0};
    outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), extent, 6);
    outInfo->Set(vtkDataObject::SPACING(), spacing, 3);
    outInfo->Set(vtkDataObject::ORIGIN(), origin, 3);
    vtkDataObject::SetPointDataActiveScalarInfo(outInfo, VTK_FLOAT, 1);

    return 1;
}


int ttkPerlinVectorField::RequestData(vtkInformation *request,
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {

  // Get output image
  vtkInformation* outInfo = outputVector->GetInformationObject(0);
  auto output = vtkImageData::GetData(outputVector);
  if (!output) {
    this->printErr("No output data available.");
    return 0;
  }

  // Set dimensions and extent
  int dimX = ImageDimension[0];
  int dimY = ImageDimension[1];
  int dimZ = ImageDimension[2];
  int extent[6] = {0, dimX - 1, 0, dimY - 1, 0, dimZ};

  // Number of tuples in image
  int nTuples = dimX * dimY;

  // Create timer for measuring execution
  ttk::Timer timer;

  // Set dimensions and extent for image data
  output->SetOrigin(0, 0, 0);
  output->SetDimensions(extent[1], extent[3], extent[5]);
  outInfo->Get
    (vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), extent);
  output->SetExtent(extent);
  
  this->printMsg(
    "Creating Perlin noise vector field of dim [" + std::to_string(dimX) + ", " + std::to_string(dimY) + ", " + std::to_string(dimZ)+"]",
    0,0,this->threadNumber_,ttk::debug::LineMode::REPLACE
  );

  // Create array to store noise in
  auto noiseArray = vtkSmartPointer<vtkDoubleArray>::New();
  noiseArray->SetName("vecField");
  noiseArray->SetNumberOfComponents(2);
  noiseArray->SetNumberOfTuples(nTuples);

  ttk::PerlinNoise pN;

  
  auto noiseData = static_cast<double*>(ttkUtils::GetVoidPointer(noiseArray));

  double dimD = (Scale > 0) ? ((double)dimX - 1.0) / Scale : ((double)dimX - 1.0) / 1;

  for (int y = 0; y < dimY; y++) {
    for (int x = 0; x < dimX; x++) {
      double xD = ((double)x)/dimD;
      double yD = ((double)y)/dimD;

      double xNoise;
      pN.perlin4D<double>(xD, yD, 0, 0, xNoise);

      // Calculate index and set noise to it in the array (nComp * x + nComp * dimX * y)
      int idx = 2 * x + 2 * dimX * y;
      noiseData[idx] = xNoise;
    }
  }

  for (int z = 0; z < dimZ; z++) {
    for (int y = 0; y < dimY; y++) {
      double yD = ((double)y)/dimD;
      double zD = ((double)z)/dimD;

      double yNoise;
      pN.perlin4D<double>(0, yD, zD, 0, yNoise);

      // Calculate index and set noise to it in the array (nComp * x + nComp * dimX * y)
      int idx = 2 * y + 2 * dimY * z;
      noiseData[idx + 1] = yNoise;
    }
  } 


  // Add arrays to image output
  output->GetPointData()->AddArray(noiseArray);


  // return success
  return 1;
}
