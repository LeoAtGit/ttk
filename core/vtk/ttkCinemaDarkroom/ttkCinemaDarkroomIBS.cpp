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

void main() {
    float depth = readDepth(vPos.xy);
    vec3 albedo = texture2D( tex1, vPos.xy ).rgb;
    float ao = texture2D( tex2, vPos.xy ).r;

    // Compute Luminance
    vec3 lumcoeff = vec3( 0.299, 0.587, 0.114 );
    vec3 luminance = vec3( dot( albedo, lumcoeff ) );

    // Silhouette Effect
    vec2 pixelSize = 1./RESOLUTION;
    vec3 eps = 1.5*vec3( pixelSize.x, pixelSize.y, 0 );
    float depthN = readDepth(vPos.xy + eps.zy);
    float depthE = readDepth(vPos.xy + eps.xz);
    float depthS = readDepth(vPos.xy - eps.zy);
    float depthW = readDepth(vPos.xy - eps.xz);
    float dxdz = (depthE-depthW);
    float dydz = (depthN-depthS);
    vec3 n = normalize( vec3(dxdz, dydz, 1./STRENGTH) );

    vec3 lightPos = vec3(0,0,1);
    float lightInt = 1.0*dot(n,normalize(lightPos));
    vec3 outputColor = vec3( albedo * mix( vec3(ao), vec3(1.0), luminance * LUMINANCE ) );
    outputColor = outputColor*AMBIENT + outputColor*lightInt;

    gl_FragColor = vec4(outputColor, depth>0.99 ? 0.0 : 1.0);
}
  )");
}

int ttkCinemaDarkroomIBS::RegisterReplacements() {
  ttkCinemaDarkroomShader::RegisterReplacements();

  this->AddReplacement("STRENGTH", {this->Strength});
  this->AddReplacement("LUMINANCE", {this->Luminance});
  this->AddReplacement("AMBIENT", {this->Ambient});
  return 1;
}

int ttkCinemaDarkroomIBS::RegisterTextures(vtkImageData *image) {
  if(!this->AddTexture(image, 0, 0))
    return 0;
  if(!this->AddTexture(image, 1, 1))
    return 0;
  if(!this->AddTexture(image, 2, 2))
    return 0;
  return 1;
}