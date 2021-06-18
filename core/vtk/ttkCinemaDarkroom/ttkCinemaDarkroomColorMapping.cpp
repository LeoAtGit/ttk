#include <ttkCinemaDarkroomColorMapping.h>

#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkUnsignedCharArray.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

#include <vtkPythonInterpreter.h>

#include <math.h>

vtkStandardNewMacro(ttkCinemaDarkroomColorMapping);

ttkCinemaDarkroomColorMapping::ttkCinemaDarkroomColorMapping() {
  this->setDebugMsgPrefix("CinemaDarkroomColorMapping");

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkCinemaDarkroomColorMapping::~ttkCinemaDarkroomColorMapping() {
}

template <typename DT>
int mapScalarsToColors(unsigned char *colors,
                      const std::vector<double> &colorMap,
                      const double *nanColor,
                      const bool transparentNAN,
                      const DT *array,
                      const double range[2],
                      const size_t nPixels,
                      const int threadNumber) {
  const size_t nKeys = colorMap.size() / 4;
  const double rangeD = range[1] - range[0];

  unsigned char nanAlpha = transparentNAN ? 0 : 255;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber)
#endif
  for(size_t i = 0; i < nPixels; i++) {
    size_t idx = i * 4;

    const double value = (double)array[i];
    if(isnan(value)) {
      colors[idx + 0] = 255.0 * nanColor[0];
      colors[idx + 1] = 255.0 * nanColor[1];
      colors[idx + 2] = 255.0 * nanColor[2];
      colors[idx + 3] = nanAlpha;
      continue;
    }

    const double normalizedValue = (value - range[0]) / rangeD;
    size_t ki = 0;
    for(size_t k = 1; k < nKeys; k++) {
      if(normalizedValue <= colorMap[k * 4]) {
        ki = k - 1;
        break;
      }
    }

    double lambda = (normalizedValue - colorMap[ki * 4]) / (colorMap[(ki + 1) * 4] - colorMap[ki * 4]);
    lambda = std::min(std::max(lambda,0.0),1.0);
    double lambdaInv = 1 - lambda;
    size_t idx2 = ki * 4;

    colors[idx + 0]
      = 255.0 * (lambdaInv * colorMap[idx2 + 1] + lambda * colorMap[idx2 + 5]);
    colors[idx + 1]
      = 255.0 * (lambdaInv * colorMap[idx2 + 2] + lambda * colorMap[idx2 + 6]);
    colors[idx + 2]
      = 255.0 * (lambdaInv * colorMap[idx2 + 3] + lambda * colorMap[idx2 + 7]);
    colors[idx + 3] = 255;
  }

  return 1;
};

int ttkCinemaDarkroomColorMapping::SyncColorMapsWithParaView(){
  std::string code(R"(
from paraview.simple import GetColorTransferFunction
from paraview.simple import GetSources
from paraview.simple import FindSource

sources = GetSources()

for sourceName in sources:
  source = FindSource(sourceName[0])
  className = source.__class__.__name__
  if (className=='TTKDarkroomColorMapping' or className=='TTKDarkroomShading') and (source.ColorMap=='ParaView' or source.ColorMap==-3):
    lut = GetColorTransferFunction(source.Scalars[1])
    cm = map(str,lut.RGBPoints)
    source.ColorMapData = ','.join(cm)
    source.NANColor = lut.NanColor
    source.ScalarRange = [0,1]
)");

  vtkPythonInterpreter::RunSimpleString(code.data());

  return 1;
}

int ttkCinemaDarkroomColorMapping::Render(vtkImageData *image) {
  ttk::Timer timer;
  int dim[3];
  image->GetDimensions(dim);
  const int nPixels = dim[0]*dim[1];

  auto scalarArray = this->GetInputArrayToProcess(0, image);
  if(!scalarArray || this->GetInputArrayAssociation(0, image) != 0
     || scalarArray->GetNumberOfComponents() != 1)
    return !this->printErr("Unable to retrieve point scalar array.");

  const std::string msg = "Mapping "+std::string(scalarArray->GetName())+" (" + std::to_string(dim[0]) + "x" + std::to_string(dim[1]) + ")";
  this->printMsg(msg, 0, 0, this->threadNumber_,ttk::debug::LineMode::REPLACE, ttk::debug::Priority::DETAIL);


  auto colorArray = vtkSmartPointer<vtkUnsignedCharArray>::New();
  colorArray->SetName((std::string("Albedo_")+scalarArray->GetName()).data());
  colorArray->SetNumberOfComponents(4);
  colorArray->SetNumberOfTuples(nPixels);
  image->GetPointData()->AddArray(colorArray);
  auto colorArrayData = ttkUtils::GetPointer<unsigned char>(colorArray);

  std::vector<double> colorMapData;
  const std::vector<double> *colorMap{nullptr};

  if(this->ColorMap == -3)
    this->SyncColorMapsWithParaView();

  if(this->ColorMap == -1) { // Solid
    colorMapData.resize(8);
    colorMapData[0] = 0;
    colorMapData[1] = this->SingleColor[0];
    colorMapData[2] = this->SingleColor[1];
    colorMapData[3] = this->SingleColor[2];
    colorMapData[4] = 1;
    colorMapData[5] = this->SingleColor[0];
    colorMapData[6] = this->SingleColor[1];
    colorMapData[7] = this->SingleColor[2];
    colorMap = &colorMapData;
  } else if(this->ColorMap == -2 || this->ColorMap == -3) {
    int status = ttkUtils::stringListToDoubleVector(
      this->ColorMapData, colorMapData);
    if(!status || colorMapData.size() < 8 || colorMapData.size() % 4 != 0)
      return !this->printErr("Invalid manual color map input.");

    colorMap = &colorMapData;
  } else {
    if(this->ColorMap < 0 || this->ColorMap >= (int)this->ColorMaps.size())
      return !this->printErr("Invalid color map index: "
                     + std::to_string(this->ColorMap));
    colorMap = &this->ColorMaps[this->ColorMap];
  }

  switch(scalarArray->GetDataType()) {
    vtkTemplateMacro(mapScalarsToColors<VTK_TT>(
      colorArrayData, *colorMap, this->NANColor, this->TransparentNAN,
      ttkUtils::GetPointer<VTK_TT>(scalarArray),
      this->ScalarRange,
      nPixels,
      this->threadNumber_));
  }

  this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_, ttk::debug::LineMode::NEW, ttk::debug::Priority::DETAIL);

  return 1;
}
