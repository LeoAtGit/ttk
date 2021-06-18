#include <ttkCinemaDarkroomPBR.h>

#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>

vtkStandardNewMacro(ttkCinemaDarkroomPBR);

ttkCinemaDarkroomPBR::ttkCinemaDarkroomPBR() : ttkCinemaDarkroomShader() {
  this->setDebugMsgPrefix("CinemaDarkroomPBR");
}

ttkCinemaDarkroomPBR::~ttkCinemaDarkroomPBR() {
}

std::string ttkCinemaDarkroomPBR::GetFragmentShaderCode() {
  return std::string(R"(
//VTK::System::Dec // always start with these lines in your FS
//VTK::Output::Dec // always start with these lines in your FS

varying vec4 vPos;

uniform sampler2D tex0; // albedo
uniform sampler2D tex1; // depth
uniform sampler2D tex2; // ao

READ_DEPTH
COMPUTE_NORMAL
COMPUTE_HALF_VECTOR

#define PI              3.14159
#define ONE_OVER_PI     0.31831

/**
 * GGX/Trowbridge-Reitz NDF
 *
 * Calculates the specular highlighting from surface roughness.
 *
 * Roughness lies on the range [0.0, 1.0], with lower values
 * producing a smoother, "glossier", surface. Higher values
 * produce a rougher surface with the specular lighting distributed
 * over a larger surface area.
 *
 * See it graphed at:
 * https://www.desmos.com/calculator/pjzk3yafzs
 */
float CalculateNDF(
  in vec3  surfNorm,
  in vec3  halfVector,
  in float roughness
){
  float a = (roughness * roughness);
  float halfAngle = dot(surfNorm, halfVector);

  return (a / (PI * pow((pow(halfAngle, 2.0) * (a - 1.0) + 1.0), 2.0)));
}

/**
 * GGX/Schlick-Beckmann microfacet geometric attenuation.
 *
 * The attenuation is modified by the roughness (input as k)
 * and approximates the influence/amount of microfacets in the surface.
 * A microfacet is a sub-pixel structure that affects light
 * reflection/occlusion.
 */
float CalculateAttenuation(
    in vec3  surfNorm,
    in vec3  vector,
    in float k
){
  float d = max(dot(surfNorm, vector), 0.0);
  return (d / ((d * (1.0 - k)) + k));
}

/**
 * GGX/Schlick-Beckmann attenuation for analytical light sources.
 */
float CalculateAttenuationAnalytical(
    in vec3  surfNorm,
    in vec3  toLight,
    in vec3  toView,
    in float roughness
){
  float k = pow((roughness + 1.0), 2.0) * 0.125;

  float lightAtten = CalculateAttenuation(surfNorm, toLight, k);
  float viewAtten  = CalculateAttenuation(surfNorm, toView, k);

  return (lightAtten * viewAtten);
}

/**
 * GGX/Schlick-Beckmann attenuation for IBL light sources.
 * Uses Disney modification of k to reduce hotness.
 */
float CalculateAttenuationIBL(
    in float roughness,
    in float normDotLight,          // Clamped to [0.0, 1.0]
    in float normDotView            // Clamped to [0.0, 1.0]
){
    float k = pow(roughness, 2.0) * 0.5;

    float lightAtten = (normDotLight / ((normDotLight * (1.0 - k)) + k));
    float viewAtten  = (normDotView / ((normDotView * (1.0 - k)) + k));

    return (lightAtten * viewAtten);
}

/**
 * Calculates the Fresnel reflectivity.
 * The metalic parameter controls the fresnel incident value (fresnel0).
 */
vec3 CalculateFresnel(
    in vec3 surfNorm,
    in vec3 toView,
    in vec3 fresnel0
){
  float d = max(dot(surfNorm, toView), 0.0);
  float p = ((-5.55473 * d) - 6.98316) * d;

  return fresnel0 + ((1.0 - fresnel0) * pow(1.0 - d, 5.0));
}

/**
 * Standard Lambertian diffuse lighting.
 */
vec3 CalculateDiffuse(
    in vec3 albedo
){
    return (albedo * ONE_OVER_PI);
}

/**
 * Cook-Torrance BRDF for analytical light sources.
 */
vec3 CalculateSpecularAnalytical(
    in    vec3  surfNorm,            // Surface normal
    in    vec3  toLight,             // Normalized vector pointing to light source
    in    vec3  toView,              // Normalized vector point to the view/camera
    in    vec3  fresnel0,            // Fresnel incidence value
    inout vec3  sfresnel,            // Final fresnel value used a kS
    in    float roughness            // Roughness parameter (microfacet contribution)
){
    vec3 halfVector = computeHalfVector(toLight, toView);

    float ndf      = CalculateNDF(surfNorm, halfVector, roughness);
    float geoAtten = CalculateAttenuationAnalytical(surfNorm, toLight, toView, roughness);

    sfresnel = CalculateFresnel(surfNorm, toView, fresnel0);

    vec3  numerator   = (sfresnel * ndf * geoAtten);
    float denominator = 4.0 * dot(surfNorm, toLight) * dot(surfNorm, toView);

    return (numerator / denominator);
}

/**
 * Calculates the total light contribution for the analytical light source.
 */
vec3 CalculateLightingAnalytical(
    in vec3  surfNorm,
    in vec3  toLight,
    in vec3  toView,
    in vec3  albedo,
    in float roughness
){
    vec3 fresnel0 = mix(vec3(0.04), albedo, METALLIC);
    vec3 ks       = vec3(0.0);
    vec3 diffuse  = CalculateDiffuse(albedo);
    vec3 specular = CalculateSpecularAnalytical(surfNorm, toLight, toView, fresnel0, ks, roughness);
    vec3 kd       = (1.0 - ks);

    float angle = clamp(dot(surfNorm, toLight), 0.0, 1.0);

    return ((kd * diffuse) + specular) * angle;
}

vec4 compute(in vec2 uv, in vec4 centerColor){

  vec4 albedoRGBA = texture2D( tex1, uv );
  float alpha = albedoRGBA.a;
  vec3 albedo = albedoRGBA.rgb;
  float ao = texture2D( tex2, uv ).r;
  float depth = readDepth(uv);

  vec3 normal = computeNormal(uv, depth);
  vec3 lightDir = normalize(vec3(1,1,1));
  // vec3 viewDir = normalize(vec3(0,0,3)-vec3(uv*2.0-1.0, -depth) );
  vec3 viewDir = normalize(vec3(0,0,1));

  vec3 ambientColor = albedo;

  vec3 aoColor = albedo*ao;

  vec3 diffuseColor = CalculateLightingAnalytical(
    normal,
    lightDir,
    viewDir,
    albedo,
    ROUGHNESS
  );

  vec3 color = ambientColor*AMBIENT + diffuseColor*DIFFUSE + aoColor*AO;

  return alpha<1.0 || depth==1.0
    ? vec4(centerColor.rgb,floor(alpha))
    : vec4(color,1.0)
  ;
}

MAIN_COLOR

  )");
}

int ttkCinemaDarkroomPBR::RegisterReplacements() {
  ttkCinemaDarkroomShader::RegisterReplacements();

  this->AddReplacement("AMBIENT", {this->Ambient});
  this->AddReplacement("DIFFUSE", {this->Diffuse});
  this->AddReplacement("AO", {this->AO});
  this->AddReplacement("ROUGHNESS", {this->Roughness});
  this->AddReplacement("METALLIC", {this->Metallic});
  return 1;
}

int ttkCinemaDarkroomPBR::RegisterTextures(vtkImageData *image) {
  if(!this->AddTexture(image, 0, 0))
    return 0;
  if(!this->AddTexture(image, 1, 1))
    return 0;
  if(!this->AddTexture(image, 2, 2))
    return 0;
  return 1;
}