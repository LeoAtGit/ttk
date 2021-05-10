#include <ttkCinemaDarkroomRendering.h>

#include <vtkInformation.h>
#include <vtkObjectFactory.h>

#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>

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
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    info->Append(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
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

  auto input = vtkDataObject::GetData(inputVector[0]);
  auto inputAsMB = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  auto outputAsMB = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  if(input->IsA("vtkMultiBlockDataSet"))
    inputAsMB->ShallowCopy(input);
  else
    inputAsMB->SetBlock(0,input);

  ttk::Timer timer;
  std::string msg = "CM -> SSSAO -> IBS";
  if(this->UseSSDoF) msg += " -> SSDoF";
  if(this->UseFXAA) msg += " -> FXAA";
  this->printMsg(msg,0,0,1,ttk::debug::LineMode::REPLACE);

  // build pipeline
  auto cm = vtkSmartPointer<ttkCinemaDarkroomColorMapping>::New();
  auto sssao = vtkSmartPointer<ttkCinemaDarkroomSSSAO>::New();
  auto ibs = vtkSmartPointer<ttkCinemaDarkroomIBS>::New();
  auto dof = vtkSmartPointer<ttkCinemaDarkroomSSDoF>::New();
  auto fxaa = vtkSmartPointer<ttkCinemaDarkroomFXAA>::New();
  ttkAlgorithm* lastUsed;
  std::string lastResult = "";

  // always perform color mapping
  {
    cm->SetInputArrayToProcess(0,this->GetInputArrayInformation(0));
    cm->SyncColorMapWithParaView();
    cm->SetManualColorMap(this->ManualColorMap);
    cm->SetColorMap(-2);
    lastUsed = cm;

    sssao->SetInputConnection(0, lastUsed->GetOutputPort(0));
    sssao->SetInputArrayToProcess(0,0,0,0,"Depth");
    sssao->SetSamples(this->Samples);
    sssao->SetRadius(this->Radius);
    sssao->SetDiffArea(this->DiffArea);
    lastUsed = sssao;

    ibs->SetInputConnection(0, lastUsed->GetOutputPort(0));
    ibs->SetInputArrayToProcess(0,0,0,0,"Diffuse");
    ibs->SetInputArrayToProcess(1,0,0,0,"Depth");
    ibs->SetInputArrayToProcess(2,0,0,0,"SSSAO");
    ibs->SetStrength(this->Strength);
    ibs->SetLuminance(this->Luminance);
    ibs->SetAmbient(this->Ambient);
    lastUsed = ibs;
    lastResult = "IBS";
  }

  if(this->UseSSDoF) {
    dof->SetInputConnection(0, lastUsed->GetOutputPort(0));
    dof->SetInputArrayToProcess(0,0,0,0,lastResult.data());
    dof->SetInputArrayToProcess(1,0,0,0,"Depth");
    dof->SetRadius(this->DepthRadius);
    dof->SetAperture(this->Aperture);
    dof->SetFocalDepth(this->FocalDepth);
    dof->SetMaxBlur(this->MaxBlur);
    lastUsed = dof;
    lastResult = "SSDoF";
  }

  if(this->UseFXAA) {
    fxaa->SetInputConnection(0, lastUsed->GetOutputPort(0));
    fxaa->SetInputArrayToProcess(0,0,0,0,lastResult.data());
    lastUsed = fxaa;
    lastResult = "FXAA";
  }

  const size_t nInputImages = inputAsMB->GetNumberOfBlocks();
  for(size_t i=0; i<nInputImages; i++){
    auto inputImage = vtkDataSet::SafeDownCast( inputAsMB->GetBlock(i) );
    cm->SetInputDataObject( inputImage );
    lastUsed->Update();

    auto outputImage = vtkImageData::SafeDownCast(lastUsed->GetOutputDataObject(0));
    auto array = outputImage->GetPointData()->GetArray(lastResult.data());
    array->SetName("Rendering");

    auto outputImageCopy = vtkSmartPointer<vtkImageData>::New();
    outputImageCopy->ShallowCopy(inputImage);
    outputImageCopy->GetPointData()->AddArray(array);

    outputAsMB->SetBlock(i, outputImageCopy);
  }

  auto output = vtkDataObject::GetData(outputVector);
  if(input->IsA("vtkMultiBlockDataSet"))
    output->ShallowCopy(outputAsMB);
  else
    output->ShallowCopy(outputAsMB->GetBlock(0));

  this->printMsg(msg,1,timer.getElapsedTime(),1);

  return 1;
}
