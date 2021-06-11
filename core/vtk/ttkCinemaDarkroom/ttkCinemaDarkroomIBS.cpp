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

#extension GL_OES_standard_derivatives : enable

varying vec4 vPos;

uniform sampler2D tex0; // color
uniform sampler2D tex1; // depth
uniform sampler2D tex2; // ao

float readDepth( const in vec2 coord ){
    return texture2D( tex1, coord ).r;
}

void main() {
    vec3 diffuse = texture2D( tex0, vPos.xy ).rgb;
    float ao = texture2D( tex2, vPos.xy ).r;
    float depth = readDepth(vPos.xy);

    // Compute Luminance
    vec3 lumcoeff = vec3( 0.299, 0.587, 0.114 );
    vec3 luminance = vec3( dot( diffuse, lumcoeff ) );

    // Silhouette Effect
    vec2 pixelSize = 1./cResolution;
    vec3 eps = vec3( pixelSize.x, pixelSize.y, 0 );
    float depthN = readDepth(vPos.xy + eps.zy);
    float depthS = readDepth(vPos.xy - eps.zy);
    float depthE = readDepth(vPos.xy + eps.xz);
    float depthW = readDepth(vPos.xy - eps.xz);

    // vec3 dx = vec3(2.0*eps.xz,depthE-depthW);
    // vec3 dy = vec3(2.0*eps.zy,depthN-depthS);

    vec3 dx = vec3(eps.xz, abs(depth-depthW) < abs(depth-depthE)
      ? depthW-depth
      : depth-depthE
    );
    vec3 dy = vec3(eps.zy, abs(depth-depthN) < abs(depth-depthS)
      ? depthN-depth
      : depth-depthS
    );

    // float dxdz = abs(depthE-depthW);
    // float dydz = abs(depthN-depthS);
    // float dxdz = abs(dFdx(depth));
    // float dydz = abs(dFdy(depth));
    // float dxdz = (depthE-depthW)/2.0;
    // float dydz = (depthN-depthS)/2.0;
    // vec3 normal = normalize( vec3(-dxdz, -dydz, 1.0) );

    // float hDepth = abs(depthE-depth) < abs(depthW-depth) ? depthE : depthW;
    // float vDepth = abs(depthN-depth) < abs(depthS-depth) ? depthN : depthS;
    // float dxdz = abs(depthE-depth) < abs(depthW-depth)
    //   ? (depthW-depth)
    //   : (depth-depthW);
    // float dydz = abs(depthN-depth) < abs(depthS-depth)
    //   ? (depthN-depth)
    //   : (depth-depthS);
    // // vec3 normal = normalize( vec3(-dxdz, -dydz, 1.0) );
    // vec3 normal = normalize(cross( vec3(dxdz,0,0), vec3(dydz,0,0)));
    // gl_FragColor = vec4(normal/2.0+0.5, depth>0.99 ? 0.0 : 1.0);


    // vec3 dx = abs(depthE-depth) < abs(depthW-depth)
    //   ? vec3( eps.xz,depth-depthE)
    //   : vec3(-eps.xz,depth-depthW);

    // vec3 dy = abs(depthN-depth) < abs(depthS-depth)
    //   ? vec3( eps.zy,depth-depthN)
    //   : vec3(-eps.zy,depth-depthS);

    vec3 normal = normalize(cross( dx, dy));
    gl_FragColor = vec4(normal/2.0+0.5, depth>0.99 ? 0.0 : 1.0);

    // vec3 pos = vec3(vPos.xy,depth);
    // vec3 normal = normalize(cross(dFdx(pos), dFdy(pos)));;
    // gl_FragColor = vec4(normal/2.0+0.5, depth>0.99 ? 0.0 : 1.0);

    // vec3 lightPos = vec3(0,-1,0);
    // float lighting = max(0,dot(lightPos,normal));

    // gl_FragColor = vec4(vec3(lighting), depth>0.99 ? 0.0 : 1.0);

    // vec3 n = normalize( vec3(dxdz, dydz, 1./cStrength) );

    // float lightInt = 1.0*dot(n,normalize(lightPos));

    // vec3 outputColor = vec3( diffuse * mix( vec3(ao), vec3(0.0), luminance * cLuminance ) );

    // outputColor = outputColor*cAmbient + outputColor*lightInt;

    // gl_FragColor = vec4(outputColor, depth>0.99 ? 0.0 : 1.0);
}
  )");
}

int ttkCinemaDarkroomIBS::RequestData(vtkInformation *request,
                                      vtkInformationVector **inputVector,
                                      vtkInformationVector *outputVector) {

  auto inputImage = vtkImageData::GetData(inputVector[0]);
  auto outputImage = vtkImageData::GetData(outputVector);
  outputImage->ShallowCopy(inputImage);

  this->InitRenderer(outputImage);

  this->AddReplacement("cStrength", {this->Strength});
  this->AddReplacement("cLuminance", {this->Luminance});
  this->AddReplacement("cAmbient", {this->Ambient});

  if(!this->AddTexture(outputImage, 0, 0))
    return 0;
  if(!this->AddTexture(outputImage, 1, 1))
    return 0;
  if(!this->AddTexture(outputImage, 2, 2))
    return 0;

  this->Render(outputImage, "IBS");

  return 1;
}