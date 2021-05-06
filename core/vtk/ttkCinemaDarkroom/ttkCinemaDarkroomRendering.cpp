#include <ttkCinemaDarkroomRendering.h>

#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>

#include <ttkCinemaDarkroomColorMapping.h>
#include <ttkCinemaDarkroomSSSAO.h>
#include <ttkCinemaDarkroomIBS.h>
#include <ttkCinemaDarkroomSSDoF.h>
#include <ttkCinemaDarkroomFXAA.h>

vtkStandardNewMacro(ttkCinemaDarkroomRendering);

ttkCinemaDarkroomRendering::ttkCinemaDarkroomRendering() : ttkAlgorithm() {
  this->setDebugMsgPrefix("CinemaDarkroomRendering");
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkCinemaDarkroomRendering::~ttkCinemaDarkroomRendering() {
}

int ttkCinemaDarkroomRendering::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
    return 1;
  }
  return 0;
}

int ttkCinemaDarkroomRendering::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkCinemaDarkroomRendering::RequestData(vtkInformation *request,
                                      vtkInformationVector **inputVector,
                                      vtkInformationVector *outputVector) {

  auto inputImage = vtkImageData::GetData(inputVector[0]);
  auto outputImage = vtkImageData::GetData(outputVector);

  vtkSmartPointer<ttkAlgorithm> lastUsed;

  if(this->UseSSSAO) {
    auto sssao = vtkSmartPointer<ttkCinemaDarkroomSSSAO>::New();
    sssao->SetInputData(inputImage);
    //input array, port number, input connection, point / cell / field data, name of the array
    sssao->SetInputArrayToProcess(0,0,0,0,"Depth");

    sssao->SetSamples(this->Samples);
    sssao->SetRadius(this->Radius);
    sssao->SetDiffArea(this->DiffArea);
    sssao->Update();

    lastUsed = sssao;
    
    if(UseIBS) {  
      auto ibs = vtkSmartPointer<ttkCinemaDarkroomIBS>::New();  
      ibs->SetInputConnection(0, sssao->GetOutputPort());

      ibs->SetInputArrayToProcess(0,0,0,0,"Diffuse");
      ibs->SetInputArrayToProcess(1,0,0,0,"Depth");
      ibs->SetInputArrayToProcess(2,0,0,0,"SSSAO");

      ibs->SetStrength(this->Strength);
      ibs->SetLuminance(this->Luminance);
      ibs->SetAmbient(this->Ambient);
      ibs->Update();

      lastUsed = ibs;
    }
  }

  if(UseDoF) {
    auto dof = vtkSmartPointer<ttkCinemaDarkroomSSDoF>::New();

    if(lastUsed == nullptr) {
      dof->SetInputData(inputImage);
    } else {
      dof->SetInputConnection(0, lastUsed->GetOutputPort(0));
    }

    dof->SetInputArrayToProcess(0,0,0,0,"Diffuse");
    dof->SetInputArrayToProcess(1,0,0,0,"Depth");

    dof->SetRadius(this->DepthRadius);
    dof->SetAperture(this->Aperture);
    dof->SetFocalDepth(this->FocalDepth);
    dof->SetMaxBlur(this->MaxBlur);

    lastUsed = dof;
  }

  if(UseFXAA) {
    auto fxaa = vtkSmartPointer<ttkCinemaDarkroomFXAA>::New();

    if(lastUsed == nullptr) {
      fxaa->SetInputData(inputImage);
    } else {
      fxaa->SetInputConnection(0, lastUsed->GetOutputPort(0));
    }

    fxaa->SetInputArrayToProcess(0,0,0,0,"Diffuse");
  
    fxaa->Update();

    lastUsed = fxaa;
  }
  if(lastUsed != nullptr) {
    outputImage->ShallowCopy(lastUsed->GetOutputDataObject(0));
  } else {
    std::cout << "Please select at least one filter to use." << std::endl;
    return -1;
  }
  
  return 1;
}
