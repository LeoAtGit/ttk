#include <ttkCinemaDarkroomIBS.h>

#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>

vtkStandardNewMacro(ttkCinemaDarkroomIBS);

ttkCinemaDarkroomIBS::ttkCinemaDarkroomIBS() : ttkCinemaDarkroomShader() {
  this->setDebugMsgPrefix("CinemaDarkroomIBS");
}

ttkCinemaDarkroomIBS::~ttkCinemaDarkroomIBS() {
}

std::string ttkCinemaDarkroomIBS::GetFragmentShaderCode() {
  return std::string(R"(
//VTK::System::Dec // always start with these lines in your FS
//VTK::Output::Dec // always start with these lines in your FS

varying vec4 vPos;

uniform sampler2D tex0; // depth
uniform sampler2D tex1; // albedo
uniform sampler2D tex2; // ao

READ_DEPTH
COMPUTE_NORMAL

vec4 compute(const in vec2 sampleUV, const in vec2 pixelUV){

    float depth = readDepth(sampleUV);
    vec3 albedo = texture2D( tex1, sampleUV ).rgb;
    float ao = texture2D( tex2, sampleUV ).r;

    // Compute Luminance
    vec3 lumcoeff = vec3( 0.299, 0.587, 0.114 );
    vec3 luminance = vec3( dot( albedo, lumcoeff ) );

    // Silhouette Effect
    vec2 pixelSize = 1./RESOLUTION;
    vec3 eps = 1.5*vec3( pixelSize.x, pixelSize.y, 0 );
    float depthN = readDepth(sampleUV + eps.zy);
    float depthE = readDepth(sampleUV + eps.xz);
    float depthS = readDepth(sampleUV - eps.zy);
    float depthW = readDepth(sampleUV - eps.xz);
    float dxdz = (depthE-depthW);
    float dydz = (depthN-depthS);
    vec3 n = normalize( vec3(dxdz, dydz, 1./STRENGTH) );

    vec3 lightPos = vec3(0,0,1);
    float lightInt = 1.0*dot(n,normalize(lightPos));
    vec3 outputColor = vec3( albedo * mix( vec3(ao), vec3(1.0), luminance * DIFFUSE ) );
    outputColor = outputColor*AMBIENT + outputColor*lightInt;

    return vec4(outputColor, depth>0.99 ? 0.0 : 1.0);
}

MAIN_COLOR

  )");
}

int ttkCinemaDarkroomIBS::RegisterReplacements() {
  ttkCinemaDarkroomShader::RegisterReplacements();

  this->AddReplacement("STRENGTH", {this->Strength});
  this->AddReplacement("DIFFUSE", {this->Diffuse});
  this->AddReplacement("AMBIENT", {this->Ambient});
  return 1;
}

int ttkCinemaDarkroomIBS::RegisterTextures(vtkImageData *image) {
  if(!this->AddTexture(image, 0))
    return !this->printErr("Unable to add texture 0.");
  if(!this->AddTexture(image, 1))
    return !this->printErr("Unable to add texture 1.");
  if(!this->AddTexture(image, 2))
    return !this->printErr("Unable to add texture 2.");
  return 1;
}