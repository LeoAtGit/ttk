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
int mapScalarsToColor(unsigned char *color,
                      const std::vector<double> &colorMap,
                      const double *nanColor,
                      const DT *array,
                      const double range[2],
                      const size_t nPixels,
                      const int threadNumber) {
  const size_t nKeys = colorMap.size() / 4;
  const double valueDelta = range[1] - range[0];

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber)
#endif
  for(size_t i = 0; i < nPixels; i++) {

    const double value = (double)array[i];
    if(isnan(value)) {
      size_t idx = i * 3;
      color[idx + 0] = 255.0 * nanColor[0];
      color[idx + 1] = 255.0 * nanColor[1];
      color[idx + 2] = 255.0 * nanColor[2];
      continue;
    }

    const double normalizedValue
      = std::max(0.0, std::min(1.0, (value - range[0]) / valueDelta));

    size_t ki = 0;
    for(size_t k = 1; k < nKeys; k++) {
      if(normalizedValue <= colorMap[k * 4]) {
        ki = k - 1;
        break;
      }
    }

    double lambda = (normalizedValue - colorMap[ki * 4])
                    / (colorMap[(ki + 1) * 4] - colorMap[ki * 4]);
    double lambdaInv = 1 - lambda;

    size_t idx = i * 3;
    size_t idx2 = ki * 4;
    color[idx + 0]
      = 255.0 * (lambdaInv * colorMap[idx2 + 1] + lambda * colorMap[idx2 + 5]);
    color[idx + 1]
      = 255.0 * (lambdaInv * colorMap[idx2 + 2] + lambda * colorMap[idx2 + 6]);
    color[idx + 2]
      = 255.0 * (lambdaInv * colorMap[idx2 + 3] + lambda * colorMap[idx2 + 7]);
  }

  return 1;
};

int ttkCinemaDarkroomColorMapping::SyncColorMapWithParaView(){
  std::string code(R"(
from paraview.simple import GetColorTransferFunction
from paraview.simple import FindSource

lut = GetColorTransferFunction('$FIELD')

cm = FindSource('TTKDarkroomColorMapping1')

if cm==None:
  cm = FindSource('TTKDarkroomRendering1')

cm.ManualColorMap = ','.join(map(str,lut.RGBPoints))

)");

  auto arrayInfo = this->GetInputArrayInformation(0);
  auto name = arrayInfo->Get(vtkDataObject::FIELD_NAME());

  code.replace(code.find("$FIELD"), 6, name);

  vtkPythonInterpreter::RunSimpleString(code.data());

  return 1;
}

int ttkCinemaDarkroomColorMapping::RequestData(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {

  ttk::Timer timer;

  auto input = vtkImageData::GetData(inputVector[0]);
  auto output = vtkImageData::GetData(outputVector);
  output->ShallowCopy(input);

  auto scalarArray = this->GetInputArrayToProcess(0, output);
  if(!scalarArray || this->GetInputArrayAssociation(0, output) != 0
     || scalarArray->GetNumberOfComponents() != 1)
    return !this->printErr("Unable to retrieve point scalar array.");

  std::string msg = "Mapping "+std::string(scalarArray->GetName());
  this->printMsg(msg, 0, 0, this->threadNumber_,ttk::debug::LineMode::REPLACE);

  size_t nPixels = scalarArray->GetNumberOfTuples();

  auto colorArray = vtkSmartPointer<vtkUnsignedCharArray>::New();
  colorArray->SetName("Diffuse");
  colorArray->SetNumberOfComponents(3);
  colorArray->SetNumberOfTuples(nPixels);
  output->GetPointData()->AddArray(colorArray);
  auto colorArrayData = ttkUtils::GetPointer<unsigned char>(colorArray);

  std::vector<double> manualColorMap;
  const std::vector<double> *colorMap{nullptr};

  double range[2];
  if(this->AutoRange && this->ColorMap != -3){
    scalarArray->GetRange(range);
  } else {
    range[0] = this->ValueRange[0];
    range[1] = this->ValueRange[1];
  }

  if(this->ColorMap == -3)
    this->SyncColorMapWithParaView();

  if(this->ColorMap == -1) { // Solid
    manualColorMap.resize(8);
    manualColorMap[0] = 0;
    manualColorMap[1] = this->SingleColor[0];
    manualColorMap[2] = this->SingleColor[1];
    manualColorMap[3] = this->SingleColor[2];
    manualColorMap[4] = 1;
    manualColorMap[5] = this->SingleColor[0];
    manualColorMap[6] = this->SingleColor[1];
    manualColorMap[7] = this->SingleColor[2];
    colorMap = &manualColorMap;
  } else if(this->ColorMap == -2 || this->ColorMap == -3) {
    int status = ttkUtils::stringListToDoubleVector(
      this->ManualColorMap, manualColorMap);
    if(!status || manualColorMap.size() < 8 || manualColorMap.size() % 4 != 0)
      return !this->printErr("Invalid manual color map input.");

    // normalize color map
    range[0] = manualColorMap[0];
    range[1] = manualColorMap[manualColorMap.size()-4];
    const auto d = range[1]-range[0];
    for(size_t i=0; i<manualColorMap.size(); i+=4)
      manualColorMap[i] = (manualColorMap[i]-range[0])/d;

    colorMap = &manualColorMap;
  } else {
    if(this->ColorMap < 0 || this->ColorMap >= (int)this->ColorMaps.size())
      return !this->printErr("Invalid color map index: "
                     + std::to_string(this->ColorMap));
    colorMap = &this->ColorMaps[this->ColorMap];
  }

  switch(scalarArray->GetDataType()) {
    vtkTemplateMacro(mapScalarsToColor<VTK_TT>(
      colorArrayData, *colorMap, this->NANColor,
      ttkUtils::GetPointer<VTK_TT>(scalarArray), range,
      nPixels, this->threadNumber_));
  }

  this->printMsg(msg, 1, timer.getElapsedTime(), this->threadNumber_);

  return 1;
}