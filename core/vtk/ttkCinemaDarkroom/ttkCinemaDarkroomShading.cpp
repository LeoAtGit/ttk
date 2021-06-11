#include <ttkCinemaDarkroomShading.h>

#include <vtkInformation.h>
#include <vtkObjectFactory.h>

#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkDataArray.h>
#include <vtkPointData.h>

#include <vtkPassThrough.h>
#include <ttkCinemaDarkroomColorMapping.h>
#include <ttkCinemaDarkroomSSSAO.h>
#include <ttkCinemaDarkroomIBS.h>
#include <ttkCinemaDarkroomSSDoF.h>
#include <ttkCinemaDarkroomFXAA.h>

vtkStandardNewMacro(ttkCinemaDarkroomShading);

ttkCinemaDarkroomShading::ttkCinemaDarkroomShading() : ttkAlgorithm() {
  this->setDebugMsgPrefix("CinemaDarkroomShading");
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkCinemaDarkroomShading::~ttkCinemaDarkroomShading() {
}

int ttkCinemaDarkroomShading::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    info->Append(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

int ttkCinemaDarkroomShading::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkCinemaDarkroomShading::RequestData(vtkInformation *request,
                                      vtkInformationVector **inputVector,
                                      vtkInformationVector *outputVector) {

  auto input = vtkDataObject::GetData(inputVector[0]);
  auto inputAsMB = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  auto outputAsMB = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  if(input->IsA("vtkMultiBlockDataSet"))
    inputAsMB->ShallowCopy(input);
  else
    inputAsMB->SetBlock(0,input);

  const size_t nInputImages = inputAsMB->GetNumberOfBlocks();
  if(nInputImages<1)
    return this->printWrn("Empty Input Image Set");
  auto firstImage = vtkImageData::SafeDownCast(inputAsMB->GetBlock(0));
  if(!firstImage)
    return !this->printErr("Input contains non-vtkImageData Objects.");
  auto firstImagePD = firstImage->GetPointData();

  ttk::Timer timer;
  std::string msg = "";

  // build pipeline
  auto passThrough = vtkSmartPointer<vtkPassThrough>::New();
  auto cm = vtkSmartPointer<ttkCinemaDarkroomColorMapping>::New();
  auto sssao = vtkSmartPointer<ttkCinemaDarkroomSSSAO>::New();
  auto ibs = vtkSmartPointer<ttkCinemaDarkroomIBS>::New();
  auto dof = vtkSmartPointer<ttkCinemaDarkroomSSDoF>::New();
  auto fxaa = vtkSmartPointer<ttkCinemaDarkroomFXAA>::New();
  vtkAlgorithm* lastUsed = passThrough;
  std::string lastResult = "";

  // initialize pipeline
  {
    if(!firstImagePD->HasArray("Diffuse")){
      cm->SetInputConnection(0, lastUsed->GetOutputPort(0));
      cm->SetInputArrayToProcess(0,this->GetInputArrayInformation(0));

      if(this->ColorMap==-3){
        cm->SetColorMap(-2);
        cm->SyncColorMapsWithParaView();
      } else
        cm->SetColorMap(this->ColorMap);

      cm->SetColorMapData(this->ColorMapData);
      cm->SetSingleColor(this->SingleColor);
      cm->SetNANColor(this->NANColor);
      cm->SetScalarRange(this->ScalarRange);

      lastUsed = cm;
      msg += " -> CM";
    }

    if(!firstImagePD->HasArray("SSSAO")){
      sssao->SetInputConnection(0, lastUsed->GetOutputPort(0));
      sssao->SetInputArrayToProcess(0,0,0,0,"Depth");
      sssao->SetSamples(this->Samples);
      sssao->SetRadius(this->Radius);
      sssao->SetDiffArea(this->DiffArea);
      lastUsed = sssao;
      msg += " -> SSSAO";
    }

    ibs->SetInputConnection(0, lastUsed->GetOutputPort(0));
    ibs->SetInputArrayToProcess(0,0,0,0,"Diffuse");
    ibs->SetInputArrayToProcess(1,0,0,0,"Depth");
    ibs->SetInputArrayToProcess(2,0,0,0,"SSSAO");
    ibs->SetStrength(this->Strength);
    ibs->SetLuminance(this->Luminance);
    ibs->SetAmbient(this->Ambient);
    lastUsed = ibs;
    lastResult = "IBS";
    msg += " -> IBS";
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
    msg += " -> SSDoF";
  }

  if(this->UseFXAA) {
    fxaa->SetInputConnection(0, lastUsed->GetOutputPort(0));
    fxaa->SetInputArrayToProcess(0,0,0,0,lastResult.data());
    lastUsed = fxaa;
    lastResult = "FXAA";
    msg += " -> FXAA";
  }
  msg = "["+std::to_string(nInputImages) + "x] " + msg.substr(4);

  this->printMsg(msg,0,0,1,ttk::debug::LineMode::REPLACE);

  // process images
  for(size_t i=0; i<nInputImages; i++){
    auto inputImage = vtkDataSet::SafeDownCast( inputAsMB->GetBlock(i) );
    passThrough->SetInputDataObject( inputImage );
    lastUsed->Update();

    auto outputImage = vtkImageData::SafeDownCast(lastUsed->GetOutputDataObject(0));
    if(!outputImage)
      return 0;

    auto array = outputImage->GetPointData()->GetArray(lastResult.data());
    if(!array)
      return 0;
    array->SetName("Shading");

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
