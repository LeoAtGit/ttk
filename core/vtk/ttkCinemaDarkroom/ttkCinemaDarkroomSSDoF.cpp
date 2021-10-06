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
POISSON_DISC

float computeCircleOfConfusion(
  const in vec2 coord
){
  float s2 = readDepth(coord);
  float c = APERTURE * abs(s2-FOCAL_DEPTH);
  return clamp(c, 0.0, MAX_BLUR);
}

vec4 compute(const in vec2 sampleUV, const in vec2 pixelUV){
    float bleedingBias = 0.02;
    float bleedingMult = 30.0;

    float centerDepth = readDepth(sampleUV);
    float centerCoC = computeCircleOfConfusion(sampleUV);

    vec4 color = vec4(0);
    float totalWeight = 0.0;

    vec2 adjustedRadius = vec2(
        RESOLUTION.y/RESOLUTION.x,
        1.0
    )*RADIUS;

    for(int i=0; i<32; i++){
        vec2 offset = poissonDisc[i] * adjustedRadius;

        vec2 sampleCoords = sampleUV + offset * centerCoC;
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

MAIN_COLOR

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
  if(!this->AddTexture(image, 0))
    return 0;
  if(!this->AddTexture(image, 1))
    return 0;
  return 1;
}