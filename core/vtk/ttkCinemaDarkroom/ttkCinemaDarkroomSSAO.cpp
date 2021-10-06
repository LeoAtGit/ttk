#include <ttkCinemaDarkroomSSAO.h>

#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>

vtkStandardNewMacro(ttkCinemaDarkroomSSAO);

ttkCinemaDarkroomSSAO::ttkCinemaDarkroomSSAO() : ttkCinemaDarkroomShader() {
  this->setDebugMsgPrefix("CinemaDarkroomSSAO");
  this->SetOutputName("SSAO");
  this->SetFloatOutput(true);
}

ttkCinemaDarkroomSSAO::~ttkCinemaDarkroomSSAO() {
}

std::string ttkCinemaDarkroomSSAO::GetFragmentShaderCode() {
  return std::string(R"(
//VTK::System::Dec // always start with these lines in your FS
//VTK::Output::Dec // always start with these lines in your FS

uniform sampler2D tex0;
varying vec4 vPos;

READ_DEPTH

COMPUTE_NORMAL
POISSON_DISC

float compute(const in vec2 sampleUV, const in vec2 pixelUV){
    float centerDepth = readDepth(sampleUV);
    vec3 pos = vec3( sampleUV, centerDepth);
    vec3 n = computeNormal(sampleUV, centerDepth);

    float occlusion = 0.0;

    vec2 aspect = vec2(RESOLUTION.y/RESOLUTION.x, 1) * RADIUS;

    for (int i = 0; i < 32; ++i){
        // get sample
        vec2 sampleTexCoord = sampleUV + poissonDisc[i] * aspect;
        vec3 samplePos      = vec3(sampleTexCoord, readDepth(sampleTexCoord));
        vec3 sampleDir      = normalize(pos - samplePos);

        // distance between SURFACE-POSITION and SAMPLE-POSITION
        float VPdistSP = distance(pos, samplePos);

        // angle between SURFACE-NORMAL and SAMPLE-DIRECTION (vector from SURFACE-POSITION to SAMPLE-POSITION)
        float dotNS = max(dot(sampleDir, n), 0.0);

        // occlusion factor
        float a = 1.0 - smoothstep(DIFF_AREA, DIFF_AREA * 2.0, VPdistSP);

        // aggregate
        occlusion += a*dotNS;
    }

    return centerDepth>0.99 ? 1.0 : 1.-occlusion/32.0;
}

MAIN_SCALAR

  )");
}

int ttkCinemaDarkroomSSAO::RegisterReplacements() {
  ttkCinemaDarkroomShader::RegisterReplacements();

  this->AddReplacement("RADIUS", {this->Radius});
  this->AddReplacement("DIFF_AREA", {this->DiffArea});
  return 1;
}

int ttkCinemaDarkroomSSAO::RegisterTextures(vtkImageData *image) {
  if(!this->AddTexture(image, 0))
    return 0;
  return 1;
}