#include <ttkCinemaImagingNative.h>

#include <vtkMultiBlockDataSet.h>
#include <vtkImageData.h>
#include <vtkPointSet.h>

#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkUnsignedIntArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>

#include <ttkCinemaImaging.h>
#include <ttkUtils.h>

#include <BVH.h>

ttk::ttkCinemaImagingNative::ttkCinemaImagingNative(){
    this->setDebugMsgPrefix("CinemaImaging(Native)");
};

ttk::ttkCinemaImagingNative::~ttkCinemaImagingNative(){
}

int ttk::ttkCinemaImagingNative::RenderVTKObject(
  vtkMultiBlockDataSet* outputImages,

  vtkPointSet* inputObject,
  vtkPointSet* inputGrid
) const {
  int status = 0;

  auto inputObjectCells = ttkCinemaImaging::GetCells(inputObject);

  // get point and cell data of input object
  auto inputObjectPD = inputObject->GetPointData();
  size_t nInputObjectPDArrays = inputObjectPD->GetNumberOfArrays();
  auto inputObjectCD = inputObject->GetCellData();
  size_t nInputObjectCDArrays = inputObjectCD->GetNumberOfArrays();

  // ---------------------------------------------------------------------------
  // Prepare Field Data for Depth Values
  // ---------------------------------------------------------------------------

  // iterate over sampling locations
  this->printMsg(ttk::debug::Separator::L2);
  float* samplingPositions = static_cast<float*>(ttkUtils::GetVoidPointer(inputGrid->GetPoints()));
  int nSamplingPositions = inputGrid->GetNumberOfPoints();
  auto camParameters = inputGrid->GetPointData();
  auto camUp = static_cast<double*>(ttkUtils::GetVoidPointer(camParameters->GetArray("CamUp")));
  auto camDir = static_cast<double*>(ttkUtils::GetVoidPointer(camParameters->GetArray("CamDirection")));
  auto camHeight = static_cast<double*>(ttkUtils::GetVoidPointer(camParameters->GetArray("CamHeight")));
  auto camAngle = static_cast<double*>(ttkUtils::GetVoidPointer(camParameters->GetArray("CamAngle")));
  auto resolution = static_cast<double*>(ttkUtils::GetVoidPointer(camParameters->GetArray("Resolution")));
  auto projectionMode = static_cast<double*>(ttkUtils::GetVoidPointer(camParameters->GetArray("ProjectionMode")));

  auto inputObjectConnectivityList = static_cast<vtkIdType*>(
      ttkUtils::GetVoidPointer(
        inputObjectCells->GetConnectivityArray()
      )
  );

  ttk::Timer test;
  BVH<vtkIdType> bvh(
    static_cast<float*>(ttkUtils::GetVoidPointer(inputObject->GetPoints())),
    inputObjectConnectivityList,
    inputObjectCells->GetNumberOfCells()
  );

  this->printMsg("BVH",1,test.getElapsedTime(),1);

  for(int i=0; i<nSamplingPositions; i++) {

      double camPos[3]{samplingPositions[i*3],samplingPositions[i*3+1],samplingPositions[i*3+2]};

      // Initialize Output
      auto outputImage = vtkSmartPointer<vtkImageData>::New();
      outputImage->SetDimensions(resolution[0],resolution[1],1);
      outputImage->SetSpacing(1,1,1);
      outputImage->SetOrigin(0,0,0);
      outputImage->AllocateScalars(VTK_FLOAT,1);

      size_t nPixels = resolution[i*2]*resolution[i*2+1];
      auto outputImagePD = outputImage->GetPointData();

      auto depthBuffer = outputImagePD->GetArray(0);
      depthBuffer->SetName("Depth");

      auto primitiveIdArray = vtkSmartPointer<vtkUnsignedIntArray>::New();
      primitiveIdArray->SetName("PrimitiveId");
      primitiveIdArray->SetNumberOfComponents(1);
      primitiveIdArray->SetNumberOfTuples( nPixels );
      outputImagePD->AddArray(primitiveIdArray);

      auto barycentricCoordinates = vtkSmartPointer<vtkFloatArray>::New();
      barycentricCoordinates->SetName("BarycentricCoordinates");
      barycentricCoordinates->SetNumberOfComponents(2);
      barycentricCoordinates->SetNumberOfTuples( nPixels );
      outputImagePD->AddArray(barycentricCoordinates);

      // Render Object
      status = this->renderImage(
        static_cast<float*>(ttkUtils::GetVoidPointer(depthBuffer)),
        static_cast<unsigned int*>(ttkUtils::GetVoidPointer(primitiveIdArray)),
        static_cast<float*>(ttkUtils::GetVoidPointer(barycentricCoordinates)),

        inputObject->GetNumberOfPoints(),
        static_cast<float*>(ttkUtils::GetVoidPointer(inputObject->GetPoints())),
        inputObjectCells->GetNumberOfCells(),
        inputObjectConnectivityList,

        bvh,

        &resolution[i*2],
        camPos,
        &camDir[i*3],
        &camUp[i*3],
        camHeight[i],
        projectionMode[i]==0
      );
      if(!status)
        return 0;

      // Map Point Data
      for(size_t j=0; j<nInputObjectPDArrays; j++){
        auto inputArray = inputObjectPD->GetArray(j);
        auto outputArray = vtkSmartPointer<vtkDataArray>::Take( inputArray->NewInstance() );
        outputArray->SetName(inputArray->GetName());
        outputArray->SetNumberOfComponents(inputArray->GetNumberOfComponents());
        outputArray->SetNumberOfTuples( nPixels );

        outputImagePD->AddArray(outputArray);

        switch(outputArray->GetDataType()){
          vtkTemplateMacro(
            status = this->interpolateArray(
              (VTK_TT*) ttkUtils::GetVoidPointer(outputArray),

              (const unsigned int*) ttkUtils::GetVoidPointer(primitiveIdArray),
              (const float*) ttkUtils::GetVoidPointer(barycentricCoordinates),
              inputObjectConnectivityList,

              (const VTK_TT*) ttkUtils::GetVoidPointer(inputArray),
              nPixels,
              inputArray->GetNumberOfComponents()
            )
          );
        }

        if(!status)
          return 0;
      }

      // Map Cell Data
      for(size_t j=0; j<nInputObjectCDArrays; j++){
        auto inputArray = inputObjectCD->GetArray(j);
        auto outputArray = vtkSmartPointer<vtkDataArray>::Take( inputArray->NewInstance() );
        outputArray->SetName(inputArray->GetName());
        outputArray->SetNumberOfComponents(inputArray->GetNumberOfComponents());
        outputArray->SetNumberOfTuples( nPixels );

        outputImagePD->AddArray(outputArray);

        switch(outputArray->GetDataType()){
          vtkTemplateMacro(
            status = this->lookupArray(
              (VTK_TT*) ttkUtils::GetVoidPointer(outputArray),

              (const unsigned int*) ttkUtils::GetVoidPointer(primitiveIdArray),

              (const VTK_TT*) ttkUtils::GetVoidPointer(inputArray),
              nPixels,
              inputArray->GetNumberOfComponents()
            )
          );
        }

        if(!status)
          return 0;
      }

      ttkCinemaImaging::AddAllFieldDataArrays(
        inputGrid,
        outputImage,
        i
      );

      outputImages->SetBlock(i, outputImage);
  }
  this->printMsg(ttk::debug::Separator::L2);

  return status;
};