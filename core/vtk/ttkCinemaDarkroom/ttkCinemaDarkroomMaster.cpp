#include <ttkCinemaDarkroomMaster.h>

#include <vtkInformation.h>
#include <vtkObjectFactory.h>

#include <vtkDataArray.h>
#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointData.h>

#include <ttkCinemaDarkroomColorMapping.h>
#include <ttkCinemaDarkroomFXAA.h>
#include <ttkCinemaDarkroomIBS.h>
#include <ttkCinemaDarkroomPBR.h>
#include <ttkCinemaDarkroomPhong.h>
#include <ttkCinemaDarkroomSSDoF.h>
#include <ttkCinemaDarkroomSSSAO.h>
#include <vtkPassThrough.h>

vtkStandardNewMacro(ttkCinemaDarkroomMaster);

ttkCinemaDarkroomMaster::ttkCinemaDarkroomMaster() : ttkAlgorithm() {
  this->setDebugMsgPrefix("CinemaDarkroomMaster");
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkCinemaDarkroomMaster::~ttkCinemaDarkroomMaster() {
}

int ttkCinemaDarkroomMaster::FillInputPortInformation(int port,
                                                      vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    info->Append(
      vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

int ttkCinemaDarkroomMaster::FillOutputPortInformation(int port,
                                                       vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

int ttkCinemaDarkroomMaster::RequestData(vtkInformation *request,
                                         vtkInformationVector **inputVector,
                                         vtkInformationVector *outputVector) {

  auto input = vtkDataObject::GetData(inputVector[0]);
  auto inputAsMB = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  if(input->IsA("vtkMultiBlockDataSet"))
    inputAsMB->ShallowCopy(input);
  else
    inputAsMB->SetBlock(0, input);

  const size_t nInputImages = inputAsMB->GetNumberOfBlocks();
  if(nInputImages < 1)
    return this->printWrn("Empty Input Image Set");
  auto firstImage = vtkImageData::SafeDownCast(inputAsMB->GetBlock(0));
  if(!firstImage)
    return !this->printErr("Input contains non-vtkImageData Objects.");
  // auto firstImagePD = firstImage->GetPointData();

  ttk::Timer timer;
  std::string msg = "";

  // build pipeline
  auto passThrough = vtkSmartPointer<vtkPassThrough>::New();
  auto cm = vtkSmartPointer<ttkCinemaDarkroomColorMapping>::New();
  auto sssao = vtkSmartPointer<ttkCinemaDarkroomSSSAO>::New();
  auto ibs = vtkSmartPointer<ttkCinemaDarkroomIBS>::New();
  auto pbr = vtkSmartPointer<ttkCinemaDarkroomPBR>::New();
  auto phong = vtkSmartPointer<ttkCinemaDarkroomPhong>::New();
  // auto dof = vtkSmartPointer<ttkCinemaDarkroomSSDoF>::New();
  auto fxaa = vtkSmartPointer<ttkCinemaDarkroomFXAA>::New();
  vtkAlgorithm *lastUsed = passThrough;
  lastUsed->SetInputDataObject(0, inputAsMB);

  // initialize pipeline
  {
    auto scalarArrayInfo = this->GetInputArrayInformation(1);
    const std::string scalarArrayName
      = scalarArrayInfo->Get(vtkDataObject::FIELD_NAME());

    if(scalarArrayName.compare("Albedo") != 0) {
      cm->SetInputConnection(0, lastUsed->GetOutputPort(0));
      cm->SetDebugLevel(this->debugLevel_);
      cm->SetInputArrayToProcess(0, scalarArrayInfo);

      if(this->ColorMap == -3) {
        cm->SetColorMap(-2);
        cm->SyncColorMapsWithParaView();
      } else
        cm->SetColorMap(this->ColorMap);

      cm->SetColorMapData(this->ColorMapData);
      cm->SetSingleColor(this->SingleColor);
      cm->SetNANColor(this->NANColor);
      cm->SetTransparentNAN(this->TransparentNAN);
      cm->SetScalarRange(this->ScalarRange);

      lastUsed = cm;
      msg += "CM";
    }

    sssao->SetDebugLevel(this->debugLevel_);
    sssao->SetInputConnection(0, lastUsed->GetOutputPort(0));
    sssao->SetInputArrayToProcess(0, this->GetInputArrayInformation(0));
    sssao->SetSamples(this->Samples);
    sssao->SetRadius(this->Radius);
    sssao->SetDiffArea(this->DiffArea);
    sssao->SetUseMSAA(this->UseMSAA);
    lastUsed = sssao;
    msg += " -> SSSAO";

    if(this->Shader == SHADER::IBS) {
      ibs->SetDebugLevel(this->debugLevel_);
      ibs->SetInputConnection(0, lastUsed->GetOutputPort(0));
      ibs->SetInputArrayToProcess(1, 0, 0, 0, "Albedo");
      ibs->SetInputArrayToProcess(0, this->GetInputArrayInformation(0));
      ibs->SetInputArrayToProcess(2, 0, 0, 0, "SSSAO");
      ibs->SetStrength(this->Strength);
      ibs->SetDiffuse(this->Diffuse);
      ibs->SetAmbient(this->Ambient);
      ibs->SetUseMSAA(this->UseMSAA);
      lastUsed = ibs;
      msg += " -> IBS";
    } else if(this->Shader == SHADER::PBR) {
      pbr->SetDebugLevel(this->debugLevel_);
      pbr->SetInputConnection(0, lastUsed->GetOutputPort(0));
      pbr->SetInputArrayToProcess(1, 0, 0, 0, "Albedo");
      pbr->SetInputArrayToProcess(0, this->GetInputArrayInformation(0));
      pbr->SetInputArrayToProcess(2, 0, 0, 0, "SSSAO");
      pbr->SetDiffuse(this->Diffuse);
      pbr->SetAmbient(this->Ambient);
      pbr->SetAO(this->AO);
      pbr->SetMetallic(this->Metallic);
      pbr->SetRoughness(this->Roughness);
      pbr->SetUseMSAA(this->UseMSAA);
      lastUsed = pbr;
      msg += " -> PBR";
    } else if(this->Shader == SHADER::PHONG) {
      phong->SetDebugLevel(this->debugLevel_);
      phong->SetInputConnection(0, lastUsed->GetOutputPort(0));
      phong->SetInputArrayToProcess(1, 0, 0, 0, "Albedo");
      phong->SetInputArrayToProcess(0, this->GetInputArrayInformation(0));
      phong->SetInputArrayToProcess(2, 0, 0, 0, "SSSAO");
      phong->SetDiffuse(this->Diffuse);
      phong->SetAmbient(this->Ambient);
      phong->SetExponent(this->Exponent);
      phong->SetSpecular(this->Specular);
      phong->SetUseMSAA(this->UseMSAA);
      lastUsed = phong;
      msg += " -> Phong";
    }

    if(this->UseFXAA) {
      fxaa->SetDebugLevel(this->debugLevel_);
      fxaa->SetInputConnection(0, lastUsed->GetOutputPort(0));
      fxaa->SetInputArrayToProcess(0, 0, 0, 0, "Shading");
      lastUsed = fxaa;
      msg += " -> FXAA";
    }
  }

  // passThrough->SetInputDataObject( inputAsMB );
  this->printMsg(msg);
  lastUsed->Update();

  auto output = vtkDataObject::GetData(outputVector);
  if(input->IsA("vtkMultiBlockDataSet"))
    output->ShallowCopy(lastUsed->GetOutputDataObject(0));
  else
    output->ShallowCopy(
      static_cast<vtkMultiBlockDataSet *>(lastUsed->GetOutputDataObject(0))
        ->GetBlock(0));

  // }

  // if(this->UseSSDoF) {
  //   dof->SetInputConnection(0, lastUsed->GetOutputPort(0));
  //   dof->SetInputArrayToProcess(0,0,0,0,lastResult.data());
  //   dof->SetInputArrayToProcess(1,0,0,0,"Depth");
  //   dof->SetRadius(this->DepthRadius);
  //   dof->SetAperture(this->Aperture);
  //   dof->SetFocalDepth(this->FocalDepth);
  //   dof->SetMaxBlur(this->MaxBlur);
  //   lastUsed = dof;
  //   lastResult = "SSDoF";
  //   msg += " -> SSDoF";
  // }

  // if(this->UseFXAA) {
  //   fxaa->SetInputConnection(0, lastUsed->GetOutputPort(0));
  //   fxaa->SetInputArrayToProcess(0,0,0,0,lastResult.data());
  //   lastUsed = fxaa;
  //   lastResult = "FXAA";
  //   msg += " -> FXAA";
  // }
  // msg = "["+std::to_string(nInputImages) + "x] " + msg.substr(4);

  // this->printMsg(msg,0,0,1,ttk::debug::LineMode::REPLACE);

  // // // process images
  // // for(size_t i=0; i<nInputImages; i++){
  // //   auto inputImage = vtkDataSet::SafeDownCast( inputAsMB->GetBlock(i) );
  // //   passThrough->SetInputDataObject( inputImage );
  // //   lastUsed->Update();

  // //   auto outputImage =
  // vtkImageData::SafeDownCast(lastUsed->GetOutputDataObject(0));
  // //   if(!outputImage)
  // //     return 0;

  // //   auto array = outputImage->GetPointData()->GetArray(lastResult.data());
  // //   if(!array)
  // //     return 0;
  // //   array->SetName("Shading");

  // //   auto outputImageCopy = vtkSmartPointer<vtkImageData>::New();
  // //   outputImageCopy->ShallowCopy(inputImage);
  // //   outputImageCopy->GetPointData()->AddArray(array);

  // //   outputAsMB->SetBlock(i, outputImageCopy);
  // // }

  // auto output = vtkDataObject::GetData(outputVector);
  // if(input->IsA("vtkMultiBlockDataSet"))
  //   output->ShallowCopy(outputAsMB);
  // else
  //   output->ShallowCopy(outputAsMB->GetBlock(0));

  // this->printMsg(msg,1,timer.getElapsedTime(),1);

  return 1;
}
