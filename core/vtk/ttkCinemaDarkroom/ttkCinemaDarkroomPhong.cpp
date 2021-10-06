#include <ttkCinemaDarkroomPhong.h>

#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>

vtkStandardNewMacro(ttkCinemaDarkroomPhong);

ttkCinemaDarkroomPhong::ttkCinemaDarkroomPhong() : ttkCinemaDarkroomShader() {
  this->setDebugMsgPrefix("CinemaDarkroomPhong");
}

ttkCinemaDarkroomPhong::~ttkCinemaDarkroomPhong() {
}

std::string ttkCinemaDarkroomPhong::GetFragmentShaderCode() {
  return std::string(R"(
//VTK::System::Dec // always start with these lines in your FS
//VTK::Output::Dec // always start with these lines in your FS

#extension GL_OES_standard_derivatives : enable

varying vec4 vPos;

uniform sampler2D tex0; // depth
uniform sampler2D tex1; // albedo
uniform sampler2D tex2; // ao

READ_DEPTH
COMPUTE_NORMAL

vec4 compute(const in vec2 sampleUV, const in vec2 pixelUV){

    vec4 pixelAlbedoRGBA = texture2D( tex1, pixelUV );
    vec4 sampleAlbedoRGBA = texture2D( tex1, sampleUV );

    vec3 albedo = sampleAlbedoRGBA.rgb;
    float alpha = sampleAlbedoRGBA.a;
    float ao = texture2D( tex2, sampleUV ).r;
    float depth = readDepth(sampleUV);
    vec3 normal = computeNormal(sampleUV, depth);

    vec3 lightDir = normalize(vec3(1,1,1));
    vec3 viewDir = vec3(0,0,1);

    vec3 halfDir = normalize(lightDir + viewDir);
    float lambertian = max(dot(lightDir, normal), 0.0);
    float specAngle = max(dot(halfDir, normal), 0.0);
    float specular = lambertian > 0.0 ? pow(specAngle, EXPONENT) : 0.0;
    vec3 ambientColor = mix(vec3(0),albedo, ao);
    vec3 diffuseColor = albedo * lambertian;

    vec3 color = ambientColor*AMBIENT + diffuseColor*DIFFUSE + specular*SPECULAR;

    return alpha<1.0 || depth==1.0
      ? vec4(pixelAlbedoRGBA.rgb,floor(alpha))
      : vec4(color,1.0)
    ;
}

MAIN_COLOR

  )");
}

int ttkCinemaDarkroomPhong::RegisterReplacements() {
  ttkCinemaDarkroomShader::RegisterReplacements();

  this->AddReplacement("EXPONENT", {this->Exponent});
  this->AddReplacement("AMBIENT", {this->Ambient});
  this->AddReplacement("DIFFUSE", {this->Diffuse});
  this->AddReplacement("SPECULAR", {this->Specular});
  return 1;
}

int ttkCinemaDarkroomPhong::RegisterTextures(vtkImageData *image) {
  if(!this->AddTexture(image, 0))
    return 0;
  if(!this->AddTexture(image, 1))
    return 0;
  if(!this->AddTexture(image, 2))
    return 0;
  return 1;
}