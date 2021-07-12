#include <ttkPerlinNoise.h>

#include <vtkInformation.h>
#include <vtkInformationVector.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkFloatArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkStreamingDemandDrivenPipeline.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkPerlinNoise);

ttkPerlinNoise::ttkPerlinNoise() {
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
}

ttkPerlinNoise::~ttkPerlinNoise() {
}


int ttkPerlinNoise::FillInputPortInformation(int port, vtkInformation *info) {
  return 0;
}

int ttkPerlinNoise::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    if (this->TimeProp < 2) {
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
    }
    else {
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    }
  }
  else
    return 0;
  return 1;
}

int ttkPerlinNoise::RequestInformation (
  vtkInformation* request, vtkInformationVector** inputVector, vtkInformationVector *outputVector
  ) {
    // For vtkImageData output we have to set extent, spacing and origin already in request information
    if (this->TimeProp < 2) {
      vtkInformation* outInfo = outputVector->GetInformationObject(0);
      int dimX, dimY, dimZ;
      dimX = dimY = dimZ = 0;
      switch (this->PerlinDim)
      {
      case 0:
        dimX = CubeDomain - 1;
        dimY = CubeDomain - 1;
        dimZ = 0;
        break;
      case 1:
        dimX = CubeDomain - 1;
        dimY = CubeDomain - 1;
        dimZ = CubeDomain - 1;
        break;
      default:
        this->printErr("Perlin noise dimension is off.");
        return 0;
      }
      int extent[6] = {0, dimX, 0, dimY, 0, dimZ};
      double spacing[3] = {1.0, 1.0, 1.0};
      double origin[3] = {0.0, 0.0, 0.0};
      outInfo->Set(vtkStreamingDemandDrivenPipeline::WHOLE_EXTENT(), extent, 6);
      outInfo->Set(vtkDataObject::SPACING(), spacing, 3);
      outInfo->Set(vtkDataObject::ORIGIN(), origin, 3);
      vtkDataObject::SetPointDataActiveScalarInfo(outInfo, VTK_FLOAT, 1);
    }

    return 1;
}

int ttkPerlinNoise::RequestData(vtkInformation *request,
                               vtkInformationVector **inputVector,
                               vtkInformationVector *outputVector) {
  
  vtkInformation* outInfo = outputVector->GetInformationObject(0);

  // Set dimensions, extent and number of tuples with values
  int cubeDim, dimX, dimY, dimZ, nTuples;
  cubeDim = dimX = dimY = this->CubeDomain;
  dimZ = nTuples = 0;
  int extent[6] = {0, dimX - 1, 0, dimY - 1, 0, -1};
  
  if (this->PerlinDim == 0) {
    dimZ = 0;
    nTuples = dimX * dimY;
    extent[5] = dimZ;
  }
  else if (this->PerlinDim == 1) {
    dimZ = cubeDim;
    nTuples = dimX * dimY * dimZ;
    extent[5] = dimZ - 1;
  }
  else {
    this->printErr("Invalid perlin dimension");
    return 0;
  }

  // Create timer for measuring execution
  ttk::Timer timer;

  // Switch on type of time-depdendency for the Perlin noise
  switch (this->TimeProp)
  {
  case 0: { // no time
    auto output = vtkImageData::GetData(outputVector);
    if (!output) {
      this->printErr("No output data available.");
      return 0;
    }

    this->printMsg(
      "Creating Perlin noise image [" + std::to_string(dimX) + ", " + std::to_string(dimY) + ", " + std::to_string(dimZ)+"]",
      0,0,this->threadNumber_,ttk::debug::LineMode::REPLACE
    );

    // Set dimensions and extent for image data
    output->SetOrigin(0, 0, 0);
    output->SetDimensions(extent[1], extent[3], extent[5]);
    outInfo->Get
      (vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), extent);
    output->SetExtent(extent);

    // Create array to store noise in
    auto noiseArray = vtkSmartPointer<vtkDoubleArray>::New();
    noiseArray->SetName("Field");
    noiseArray->SetNumberOfComponents(1);
    noiseArray->SetNumberOfTuples(nTuples);

    if (this->PerlinDim == 0) {
      // Calculate 2D noise for image
      int dims[2] = {dimX, dimY};
      switch(noiseArray->GetDataType()) {
        vtkTemplateMacro(this->perlin2Daux<VTK_TT>(
          dims,
          this->nOctaves,
          this->Scale,
          this->Frequency,
          this->Persistence,
          static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(noiseArray))
        ));
      }
    }
    else {
      // Calculate 3D noise for image
      int dims[3] = {dimX, dimY, dimZ};
      switch(noiseArray->GetDataType()) {
        vtkTemplateMacro(this->perlin3Daux<VTK_TT>(
          dims,
          this->nOctaves,
          this->Scale,
          this->Frequency,
          this->Persistence,
          static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(noiseArray))
        ));
      }
    }
    output->GetPointData()->AddArray(noiseArray);

    this->printMsg(
      "Creating Perlin noise image [" + std::to_string(dimX) + ", " + std::to_string(dimY) + ", " + std::to_string(dimZ)+"]",
      1,timer.getElapsedTime(),this->threadNumber_
    );
    timer.reStart();

    break;
  }
  case 1: { // singular time-step
    auto output = vtkImageData::GetData(outputVector);
    if (!output) {
      this->printErr("No output data available.");
      return 0;
    }

    this->printMsg(
      "Creating Perlin noise image [" + std::to_string(dimX) + ", " + std::to_string(dimY) + ", " + std::to_string(dimZ)+"], t = " + std::to_string(this->TimeStep),
      0,0,this->threadNumber_,ttk::debug::LineMode::REPLACE
    );

    // Set dimensions and extent for image data
    output->SetOrigin(0, 0, 0);
    output->SetDimensions(extent[1], extent[3], extent[5]);
    outInfo->Get
      (vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), extent);
    output->SetExtent(extent);

    // Create array to store noise in
    auto noiseArray = vtkSmartPointer<vtkDoubleArray>::New();
    noiseArray->SetName("Field");
    noiseArray->SetNumberOfComponents(1);
    noiseArray->SetNumberOfTuples(nTuples);

    // Create an array to store time-step in
    auto tsArray = vtkSmartPointer<vtkDoubleArray>::New();
    tsArray->SetName("Time");
    tsArray->SetNumberOfComponents(1);
    tsArray->InsertNextValue(this->TimeStep);

    // Check perlin dimension
    if (this->PerlinDim == 0) {
      // Get 2D+T noise for chosen time-step
      int dims[2] = {dimX, dimY};
      switch(noiseArray->GetDataType()) {
        vtkTemplateMacro(this->perlin2DTaux<VTK_TT>(
          dims,
          this->TimeStep,
          this->nOctaves,
          this->Scale,
          this->Frequency,
          this->Persistence,
          static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(noiseArray))
        ));
      }
    }
    else {
      // Get 3D+T noise for chosen time-step
      int dims[3] = {dimX, dimY, dimZ};
      switch(noiseArray->GetDataType()) {
        vtkTemplateMacro(this->perlin3DTaux<VTK_TT>(
          dims,
          this->TimeStep,
          this->nOctaves,
          this->Scale,
          this->Frequency,
          this->Persistence,
          static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(noiseArray))
        ));
      }
    }
    // Add arrays to image output
    output->GetPointData()->AddArray(noiseArray);
    output->GetFieldData()->AddArray(tsArray);

    this->printMsg(
      "Creating Perlin noise image [" + std::to_string(dimX) + ", " + std::to_string(dimY) + ", " + std::to_string(dimZ)+"], t = " + std::to_string(this->TimeStep),
      1,timer.getElapsedTime(),this->threadNumber_
    );
    timer.reStart();

    break;
  }
  case 2: { // Outputs a time-series
    auto output = vtkMultiBlockDataSet::GetData(outputVector);
    if (!output) {
      this->printErr("No output data available.");
      return 0;
    }

    // Loop through all time-steps
    this->printMsg(
      "Creating time series of " + std::to_string(this->TimeSeries) + " timesteps, interval = " + std::to_string(this->Interval) +", of Perlin noise images [" + std::to_string(dimX) + ", " + std::to_string(dimY) + ", " + std::to_string(dimZ)+"]",
      0,0,this->threadNumber_, ttk::debug::LineMode::REPLACE
    );
    for (int t = 0; t < this->TimeSeries; t++) {
      // Calculate actual time-step with Interval
      double timeStep = this->Interval * t;

      // Create a VTK image for current time-step using smart pointers
      auto image = vtkSmartPointer<vtkImageData>::New();

      // Set image origin, dims and extent
      image->SetOrigin(0, 0, 0);
      image->SetDimensions(extent[1], extent[3], extent[5]);
      image->SetExtent(extent);

      // Create an array to store the noise in and set attributes
      auto noiseArray = vtkSmartPointer<vtkDoubleArray>::New();
      noiseArray->SetName("Field");
      noiseArray->SetNumberOfComponents(1);
      noiseArray->SetNumberOfTuples(nTuples);

      // Create an array to store time-step in
      auto tsArray = vtkSmartPointer<vtkDoubleArray>::New();
      tsArray->SetName("Time");
      tsArray->SetNumberOfComponents(1);
      tsArray->InsertNextValue(timeStep);


      // Check perlin dimension
      if (this->PerlinDim == 0) {
        int dims[2] = {dimX, dimY};

        // Execute perlin for time-step
        switch(noiseArray->GetDataType()) {
          vtkTemplateMacro(this->perlin2DTaux<VTK_TT>(
            dims,
            timeStep,
            this->nOctaves,
            this->Scale,
            this->Frequency,
            this->Persistence,
            static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(noiseArray))
          ));
        }
      }
      else {
        int dims[3] = {dimX, dimY, dimZ};

        // Execute perlin for time-step
        switch(noiseArray->GetDataType()) {
          vtkTemplateMacro(this->perlin3DTaux<VTK_TT>(
            dims,
            timeStep,
            this->nOctaves,
            this->Scale,
            this->Frequency,
            this->Persistence,
            static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(noiseArray))
          ));
        }
      }

      // Add arrays to image
      image->GetPointData()->AddArray(noiseArray);
      image->GetFieldData()->AddArray(tsArray);

      // Set image to a block in the output dataset
      size_t nBlocks = output->GetNumberOfBlocks();
      output->SetBlock(nBlocks, image); 
    }
    this->printMsg(
      "Creating time series of " + std::to_string(this->TimeSeries) + " timesteps, interval = " + std::to_string(this->Interval) +", of Perlin noise images [" + std::to_string(dimX) + ", " + std::to_string(dimY) + ", " + std::to_string(dimZ)+"]",
      1,timer.getElapsedTime(),this->threadNumber_
    );
    timer.reStart();
    break;
  }
  default: {
    this->printErr("No time option selected.");
    break;
  }
  }

  // return success
  return 1;
}
