#include <ttkCinemaDarkroomShader.h>

#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointData.h>

#include <vtkCamera.h>
#include <vtkOpenGLPolyDataMapper.h>
#include <vtkOpenGLRenderWindow.h>
#include <vtkOpenGLRenderer.h>
#include <vtkPlaneSource.h>
#include <vtkRenderWindow.h>
#include <vtkWindowToImageFilter.h>

#include <vtkFloatArray.h>
#include <vtkUnsignedCharArray.h>

#include <vtkOpenGLTexture.h>
#include <vtkProperty.h>
#include <vtkShaderProperty.h>
#include <vtkTextureObject.h>

#include <vtkFramebufferPass.h>
#include <vtkRenderStepsPass.h>

#include <boost/algorithm/string.hpp>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkCinemaDarkroomShader);

const std::string mainColorCode(R"(
  void main(){
    vec2 uv = vPos.xy;
    vec4 color = compute(uv,texture2D( tex1, uv ));
    gl_FragColor = color;
  }
)");

const std::string mainColorMSAACode(R"(
  void main(){
    vec2 pixelSize = 1.0/RESOLUTION;

    vec2 uv = vPos.xy;
    vec2 eps = vec2(1,-1)*0.15;
    vec4 centerColor = texture2D( tex1, uv );

    vec4 color = (
       compute(uv+eps.xx*pixelSize, centerColor)
      +compute(uv+eps.xy*pixelSize, centerColor)
      +compute(uv+eps.yy*pixelSize, centerColor)
      +compute(uv+eps.yx*pixelSize, centerColor)
    );
    color = color/4.0;

    gl_FragColor = color;
  }
)");

const std::string mainScalarCode(R"(
  void main(){
    vec2 uv = vPos.xy;

    float scalar = compute(uv);

    gl_FragDepth = scalar;
  }
)");

const std::string mainScalarMSAACode(R"(
  void main(){
    vec2 uv = vPos.xy;

    vec2 pixelSizeHalf = 0.5/RESOLUTION;
    vec2 eps = vec2(0.5,-0.5);

    float scalar = (
       compute(uv+eps.xx*pixelSizeHalf)
      +compute(uv+eps.xy*pixelSizeHalf)
      +compute(uv+eps.yy*pixelSizeHalf)
      +compute(uv+eps.yx*pixelSizeHalf)
    )/4.0;

    gl_FragDepth = scalar;
  }
)");

const std::string readDepthCode(R"(
  float readDepth( const in vec2 coord ){
    return texture2D( tex0, coord ).r;
  }
)");

const std::string computeNormalCode(R"(
  vec3 computeNormal(in vec2 uv, in float depth){
    vec2 pixelSize = 1./RESOLUTION;
    vec3 eps = vec3( pixelSize.x, pixelSize.y, 0 );
    float depthN = readDepth(uv.xy + eps.zy);
    float depthS = readDepth(uv.xy - eps.zy);
    float depthE = readDepth(uv.xy + eps.xz);
    float depthW = readDepth(uv.xy - eps.xz);

    // vec3 dx = vec3(2.0*eps.xz,depthE-depthW);
    // vec3 dy = vec3(2.0*eps.zy,depthN-depthS);

    vec3 dx = vec3(eps.xz, abs(depth-depthW) < abs(depth-depthE)
      ? depthW-depth
      : depth-depthE
    );
    vec3 dy = vec3(eps.zy, abs(depth-depthN) < abs(depth-depthS)
      ? depth-depthN
      : depthS-depth
    );

    return normalize(cross(dx, dy));
  }
)");

const std::string computeHalfVectorCode(R"(
  vec3 computeHalfVector(
    in vec3 toLight,
    in vec3 toView
  ){
    return normalize(toLight + toView);
  }
)");

const std::string packFloatCode(R"(
  float shift_right (float v, float amt) {
      v = floor(v) + 0.5;
      return floor(v / exp2(amt));
  }
  float shift_left (float v, float amt) {
      return floor(v * exp2(amt) + 0.5);
  }
  float mask_last (float v, float bits) {
      return mod(v, shift_left(1.0, bits));
  }
  float extract_bits (float num, float from, float to) {
      from = floor(from + 0.5); to = floor(to + 0.5);
      return mask_last(shift_right(num, from), to - from);
  }
  vec4 encode_float(float val) {
      if (val == 0.0) return vec4(0, 0, 0, 0);
      float sign = val > 0.0 ? 0.0 : 1.0;
      val = abs(val);
      float exponent = floor(log2(val));
      float biased_exponent = exponent + 127.0;
      float fraction = ((val / exp2(exponent)) - 1.0) * 8388608.0;
      float t = biased_exponent / 2.0;
      float last_bit_of_biased_exponent = fract(t) * 2.0;
      float remaining_bits_of_biased_exponent = floor(t);
      float byte4 = extract_bits(fraction, 0.0, 8.0) / 255.0;
      float byte3 = extract_bits(fraction, 8.0, 16.0) / 255.0;
      float byte2 = (last_bit_of_biased_exponent * 128.0 + extract_bits(fraction, 16.0, 23.0)) / 255.0;
      float byte1 = (sign * 128.0 + remaining_bits_of_biased_exponent) / 255.0;
      return vec4(byte4, byte3, byte2, byte1);
  }
)");

ttkCinemaDarkroomShader::ttkCinemaDarkroomShader() {
  this->CreateFullScreenQuad();
  this->CreateRenderer();

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkCinemaDarkroomShader::~ttkCinemaDarkroomShader() {
}

int ttkCinemaDarkroomShader::FillInputPortInformation(int port,
                                                      vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkImageData");
    info->Append(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

int ttkCinemaDarkroomShader::FillOutputPortInformation(int port,
                                                       vtkInformation *info) {
  if(port == 0) {
    info->Set(ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT(), 0);
    return 1;
  }
  return 0;
}

std::string ttkCinemaDarkroomShader::GetVertexShaderCode() {
  return std::string(R"(

//VTK::System::Dec  // always start with this line in your VS

attribute vec4 vertexMC;
varying vec4 vPos;

void main () {
    vPos = vertexMC/2. + vec4(0.5,0.5,0.5,0);
    gl_Position = vertexMC;
}

  )");
}

std::string ttkCinemaDarkroomShader::GetFragmentShaderCode() {
  return std::string(R"(

//VTK::System::Dec // always start with these lines in your FS
//VTK::Output::Dec // always start with these lines in your FS

varying vec4 vPos;

void main(void) {
    gl_FragData[0] = vec4(1,0,0,1);
}

  )");
}

std::string
  ttkCinemaDarkroomShader::PerformReplacements(const std::string &source) {
  std::string result = source;

  for(int i=this->Replacements.size()-1; i>=0; i--){
    const auto& it = this->Replacements[i];
    boost::replace_all(result, it.first, it.second);
  }

  return result;
}

int ttkCinemaDarkroomShader::AddReplacement(const std::string &name,
                                            const std::vector<double> &values,
                                            const bool &isInt) {
  std::string text;
  if(values.size() == 0)
    return 0;

  if(values.size() > 1) {
    if(isInt)
      text += "i";

    text += "vec" + std::to_string(values.size()) + "(";
  }

  if(isInt)
    text += std::to_string((int)values[0]);
  else
    text += std::to_string(values[0]);

  for(size_t i = 1; i < values.size(); i++)
    if(isInt)
      text += "," + std::to_string((int)values[i]);
    else
      text += "," + std::to_string(values[i]);

  if(values.size() > 1)
    text += ")";

  return this->AddReplacementText(name, text);
}

int ttkCinemaDarkroomShader::AddReplacementText(const std::string &name, const std::string &text){
  this->Replacements.resize(this->Replacements.size()+1);
  auto& pair = this->Replacements.back();
  pair.first = name;
  pair.second = text;

  return 1;
}

int ttkCinemaDarkroomShader::CreateFullScreenQuad() {
  auto ps = vtkSmartPointer<vtkPlaneSource>::New();
  ps->SetOrigin(-1, -1, 0);
  ps->SetPoint1(1, -1, 0);
  ps->SetPoint2(-1, 1, 0);
  ps->Update();

  this->FullScreenQuad = vtkSmartPointer<vtkPolyData>::New();
  this->FullScreenQuad->ShallowCopy(ps->GetOutput());
  this->FullScreenQuadActor = vtkSmartPointer<vtkActor>::New();

  auto mapper = vtkSmartPointer<vtkOpenGLPolyDataMapper>::New();
  mapper->SetInputData(this->FullScreenQuad);
  this->FullScreenQuadActor->SetMapper(mapper);

  return 1;
}

int ttkCinemaDarkroomShader::CreateRenderer() {
  // Renderer
  this->Renderer = vtkSmartPointer<vtkOpenGLRenderer>::New();
  this->Renderer->AddActor(this->FullScreenQuadActor);
  this->Renderer->SetBackground(0, 0, 0);

  // Camera
  auto camera = vtkSmartPointer<vtkCamera>::New();
  camera->SetParallelProjection(true);
  camera->SetClippingRange(0, 2);
  camera->SetPosition(0, 0, 1);
  camera->SetFocalPoint(0, 0, 0);
  camera->SetParallelScale(
    1); // Will be ignored because quad positions are fixed
  this->Renderer->SetActiveCamera(camera);

  this->RenderWindow = vtkSmartPointer<vtkRenderWindow>::New();

  return 1;
}

int ttkCinemaDarkroomShader::InitRenderer(vtkImageData *outputImage) {
  int dim[3];
  outputImage->GetDimensions(dim);

  // Window
  auto size = this->RenderWindow->GetSize();
  if(size[0] != dim[0] || size[1] != dim[1]) {
    ttk::Timer timer;
    this->printMsg("Initializing Renderer (" + std::to_string(dim[0]) + "x"
                     + std::to_string(dim[1]) + ")",
                   0, 0, 1, ttk::debug::LineMode::REPLACE);

    this->RenderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    this->RenderWindow->AddRenderer(this->Renderer);
    this->RenderWindow->SetMultiSamples(0); // Disable AA
    this->RenderWindow->OffScreenRenderingOn();

    auto windowAsOGL = vtkOpenGLRenderWindow::SafeDownCast(this->RenderWindow);
    windowAsOGL->SetSize(dim[0], dim[1]);
    windowAsOGL->Initialize();

    this->printMsg("Initializing Renderer (" + std::to_string(dim[0]) + "x"
                     + std::to_string(dim[1]) + ")",
                   1, timer.getElapsedTime(), 1);
  }

  return 1;
}

int ttkCinemaDarkroomShader::RegisterTextures(vtkImageData *image){
  return 1;
}

int ttkCinemaDarkroomShader::RegisterReplacements(){
  this->Replacements.clear();
  auto size = this->RenderWindow->GetSize();
  this->AddReplacement("RESOLUTION", {(double)size[0], (double)size[1]});
  this->AddReplacementText("READ_DEPTH", readDepthCode);
  this->AddReplacementText("COMPUTE_NORMAL", computeNormalCode);
  this->AddReplacementText("COMPUTE_HALF_VECTOR", computeHalfVectorCode);

  this->AddReplacementText("MAIN_COLOR", this->UseMSAA ? mainColorMSAACode : mainColorCode);
  this->AddReplacementText("MAIN_COLOR_MSAA", mainColorMSAACode);
  this->AddReplacementText("MAIN_SCALAR", this->UseMSAA ? mainScalarMSAACode : mainScalarCode);
  this->AddReplacementText("MAIN_SCALAR_MSAA", mainScalarMSAACode);
  return 1;
}

int ttkCinemaDarkroomShader::AddTexture(vtkImageData *image,
                                        int arrayIdx,
                                        int textureIdx) {
  int dim[3];
  image->GetDimensions(dim);

  auto inputArray = this->GetInputArrayToProcess(arrayIdx, image);
  if(!inputArray || this->GetInputArrayAssociation(arrayIdx, image) != 0) {
    this->printErr("Unable to retrieve input point data array "
                   + std::to_string(arrayIdx) + ".");
    return 0;
  }

  std::string textureName = "tex" + std::to_string(textureIdx);

  auto properties = this->FullScreenQuadActor->GetProperty();
  if(!properties->GetTexture(textureName.data())){
    auto texture = vtkSmartPointer<vtkOpenGLTexture>::New();
    texture->InterpolateOn();
    properties->SetTexture(textureName.data(), texture);
  }

  auto texture = static_cast<vtkOpenGLTexture*>(properties->GetTexture(textureName.data()));
  auto textureObj = vtkSmartPointer<vtkTextureObject>::New();
  textureObj->SetContext(
    vtkOpenGLRenderWindow::SafeDownCast(this->RenderWindow));
  textureObj->SetWrapT(vtkTextureObject::ClampToEdge);
  textureObj->SetWrapS(vtkTextureObject::ClampToEdge);
  // textureObj->SetLinearMagnification(inputArray->GetNumberOfComponents()!=4);
  textureObj->SetLinearMagnification(true);
  textureObj->Create2DFromRaw(
    dim[0], dim[1], inputArray->GetNumberOfComponents(),
    inputArray->GetDataType(), ttkUtils::GetVoidPointer(inputArray));

  texture->SetTextureObject(textureObj);

  return 1;
}

int ttkCinemaDarkroomShader::Render(vtkImageData *image) {
  ttk::Timer timer;
  int dim[3];
  image->GetDimensions(dim);
  const int nPixels = dim[0]*dim[1];

  this->printMsg(
    "Rendering (" + std::to_string(dim[0]) + "x" + std::to_string(dim[1]) + ")",
    0, 0, 1, ttk::debug::LineMode::REPLACE, ttk::debug::Priority::DETAIL);

  this->RenderWindow->Render();

  vtkSmartPointer<vtkDataArray> outputArray;
  if(this->FloatOutput){
    outputArray = vtkSmartPointer<vtkFloatArray>::New();
    outputArray->SetNumberOfComponents(1);
    outputArray->SetNumberOfTuples(nPixels);
    this->RenderWindow->GetZbufferData(
      0, 0, dim[0] - 1, dim[1] - 1, vtkFloatArray::SafeDownCast(outputArray));
  } else {
    outputArray = vtkSmartPointer<vtkUnsignedCharArray>::New();
    outputArray->SetNumberOfComponents(4);
    outputArray->SetNumberOfTuples(nPixels);
    this->RenderWindow->GetRGBACharPixelData(
      0, 0, dim[0] - 1, dim[1] - 1, 1, vtkUnsignedCharArray::SafeDownCast(outputArray));
  }
  outputArray->SetName(this->GetOutputName().data());

  image->GetPointData()->AddArray(outputArray);
  image->GetPointData()->SetActiveScalars(outputArray->GetName());

  this->printMsg(
    "Rendering (" + std::to_string(dim[0]) + "x" + std::to_string(dim[1]) + ")",
    1, timer.getElapsedTime(), 1, ttk::debug::LineMode::NEW, ttk::debug::Priority::DETAIL);

  return 1;
}

int ttkCinemaDarkroomShader::RequestData(vtkInformation *request,
                                      vtkInformationVector **inputVector,
                                      vtkInformationVector *outputVector) {

  // fetch input
  auto input = vtkDataObject::GetData(inputVector[0]);
  auto inputAsMB = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  if(input->IsA("vtkMultiBlockDataSet")){
    inputAsMB->ShallowCopy(input);
  } else if(input->IsA("vtkImageData")) {
    inputAsMB->SetBlock(0,input);
  }

  const size_t nImages = inputAsMB->GetNumberOfBlocks();
  if(nImages<1)
    return this->printWrn("Empty Input Image Set");

  // initialize renderer based on first image
  auto firstImage = vtkImageData::SafeDownCast(inputAsMB->GetBlock(0));
  this->InitRenderer(firstImage);

  // prepare shaders
  this->RegisterReplacements();
  this->FullScreenQuadActor->GetShaderProperty()->SetVertexShaderCode(
    this->PerformReplacements(this->GetVertexShaderCode()).data()
  );
  this->FullScreenQuadActor->GetShaderProperty()->SetFragmentShaderCode(
    this->PerformReplacements(this->GetFragmentShaderCode()).data()
  );

  // compute output
  ttk::Timer timer;
  int dim[3];
  firstImage->GetDimensions(dim);
  const std::string msg = "Rendering "+std::to_string(nImages)+" (" + std::to_string(dim[0]) + "x" + std::to_string(dim[1]) + ") Images";
  auto debugLevel = this->debugLevel_==static_cast<int>(ttk::debug::Priority::DETAIL) ? ttk::debug::Priority::VERBOSE : static_cast<ttk::debug::Priority>(this->debugLevel_);

  this->printMsg(msg, 0, 0, 1, ttk::debug::LineMode::REPLACE, debugLevel);
  auto outputAsMB = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  for(size_t i=0; i<nImages; i++){
    // init output
    auto outputImage = vtkSmartPointer<vtkImageData>::New();
    outputImage->ShallowCopy(inputAsMB->GetBlock(i));
    outputAsMB->SetBlock(i, outputImage);

    // register textures
    this->RegisterTextures(outputImage);

    // render image
    this->Render(outputImage);

    // print progress
    this->printMsg(msg, ((float)i+1)/((float)nImages), timer.getElapsedTime(), 1, ttk::debug::LineMode::REPLACE, debugLevel);
  }
  this->printMsg(msg, 1, timer.getElapsedTime(), 1, ttk::debug::LineMode::NEW, debugLevel);

  auto output = vtkDataObject::GetData(outputVector);
  if(output->IsA("vtkMultiBlockDataSet"))
    output->ShallowCopy(outputAsMB);
  else
    output->ShallowCopy(outputAsMB->GetBlock(0));

  return 1;
}