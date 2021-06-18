#include <ttkCinemaDarkroomSSDoF.h>

#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>

vtkStandardNewMacro(ttkCinemaDarkroomSSDoF);

ttkCinemaDarkroomSSDoF::ttkCinemaDarkroomSSDoF() : ttkCinemaDarkroomShader() {
  this->setDebugMsgPrefix("CinemaDarkroomSSDoF");
}

ttkCinemaDarkroomSSDoF::~ttkCinemaDarkroomSSDoF() {
}

std::string ttkCinemaDarkroomSSDoF::GetFragmentShaderCode() {
  return std::string(R"(
//VTK::System::Dec // always start with these lines in your FS
//VTK::Output::Dec // always start with these lines in your FS

varying vec4 vPos;

uniform sampler2D tex0; // depth
uniform sampler2D tex1; // color

READ_DEPTH

float computeCircleOfConfusion(
  const in vec2 coord
){
  float s2 = readDepth(coord);
  float c = APERTURE * abs(s2-FOCAL_DEPTH);
  return clamp(c, 0.0, MAX_BLUR);
}

const vec2 poissonDisk[32] = vec2[](
  vec2( -0.94201624,  -0.39906216 ),
  vec2(  0.94558609,  -0.76890725 ),
  vec2( -0.094184101, -0.92938870 ),
  vec2(  0.34495938,   0.29387760 ),
  vec2( -0.91588581,   0.45771432 ),
  vec2( -0.81544232,  -0.87912464 ),
  vec2( -0.38277543,   0.27676845 ),
  vec2(  0.97484398,   0.75648379 ),
  vec2(  0.44323325,  -0.97511554 ),
  vec2(  0.53742981,  -0.47373420 ),
  vec2( -0.26496911,  -0.41893023 ),
  vec2(  0.79197514,   0.19090188 ),
  vec2( -0.24188840,   0.99706507 ),
  vec2( -0.81409955,   0.91437590 ),
  vec2(  0.19984126,   0.78641367 ),
  vec2(  0.14383161,  -0.14100790 ),
  vec2( -0.44201624,  -0.29906216 ),
  vec2(  0.94558609,  -0.46890725 ),
  vec2( -0.194184101, -0.42938870 ),
  vec2(  0.24495938,   0.99387760 ),
  vec2( -0.31588581,   0.45771432 ),
  vec2( -0.81544232,  -0.87912464 ),
  vec2( -0.08277543,   0.87676845 ),
  vec2(  0.57484398,   0.55648379 ),
  vec2(  0.74323325,  -0.27511554 ),
  vec2(  0.44298431,  -0.47373420 ),
  vec2( -0.21196911,  -0.22893023 ),
  vec2(  0.79197514,   0.12020188 ),
  vec2( -0.11184840,   0.99706507 ),
  vec2( -0.4309955,   0.111437590 ),
  vec2(  0.12344126,   0.78641367 ),
  vec2(  0.2183161,   -0.89100790 )
);

vec4 SSDoF(
    const in vec2 coord
){
    float bleedingBias = 0.02;
    float bleedingMult = 30.0;

    float centerDepth = readDepth(coord);
    float centerCoC = computeCircleOfConfusion(coord);

    vec4 color = vec4(0);
    float totalWeight = 0.0;

    vec2 adjustedRadius = vec2(
        RESOLUTION.y/RESOLUTION.x,
        1.0
    )*RADIUS;

    for(int i=0; i<32; i++){
        vec2 offset = poissonDisk[i] * adjustedRadius;

        vec2 sampleCoords = coord + offset * centerCoC;
        float sampleCoC = computeCircleOfConfusion(sampleCoords);

        vec4 samplePixel = texture2D(tex1, sampleCoords);
        float sampleDepth = readDepth(sampleCoords);

        float weight = sampleDepth < centerDepth ? sampleCoC * bleedingMult : 1.0;
        weight = (centerCoC > sampleCoC + bleedingBias) ? weight : 1.0;
        weight = clamp(weight,0.0,1.0);

        color += samplePixel*weight;
        totalWeight += weight;
    }

    return color / totalWeight;
}

void main() {
    gl_FragColor = SSDoF(vPos.xy);
}
  )");
}

int ttkCinemaDarkroomSSDoF::RegisterReplacements() {
  ttkCinemaDarkroomShader::RegisterReplacements();

  this->AddReplacement("RADIUS", {this->Radius});
  this->AddReplacement("MAX_BLUR", {this->MaxBlur});
  this->AddReplacement("APERTURE", {this->Aperture});
  this->AddReplacement("FOCAL_DEPTH", {this->FocalDepth});
  return 1;
}

int ttkCinemaDarkroomSSDoF::RegisterTextures(vtkImageData *image) {
  if(!this->AddTexture(image, 0, 0))
    return 0;
  if(!this->AddTexture(image, 1, 1))
    return 0;
  return 1;
}