#include <ttkEventDataConverter.h>

#include <vtkDataObject.h> // For port info
#include <vtkObjectFactory.h> // for new macro

#include <sstream>
#include <string>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkStringArray.h>
#include <vtkUnsignedCharArray.h>

#include <ttkUtils.h>

using namespace std;

vtkStandardNewMacro(ttkEventDataConverter);

ttkEventDataConverter::ttkEventDataConverter() {
  this->SetFilepath("");

  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
}

ttkEventDataConverter::~ttkEventDataConverter() {
}

// see ttkAlgorithm::FillInputPortInformation for details about this method
int ttkEventDataConverter::FillInputPortInformation(int port,
                                                    vtkInformation *info) {
  return 0;
}

// see ttkAlgorithm::FillOutputPortInformation for details about this method
int ttkEventDataConverter::FillOutputPortInformation(int port,
                                                     vtkInformation *info) {
  switch(port) {
    case 0:
      info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
      break;
    default:
      return 0;
  }
  return 1;
}

int ttkEventDataConverter::RequestData(vtkInformation *request,
                                       vtkInformationVector **inputVector,
                                       vtkInformationVector *outputVector) {
  this->printMsg("VTK Layer infos go here");
  this->printMsg(ttk::debug::Separator::L1);

  std::string path = this->GetFilepath();
  this->printMsg("path: " + path);
  if(path.empty()) {
    return 0;
  }

  // Get the output
  auto output = vtkUnstructuredGrid::GetData(outputVector);

  // TODO: compute buffer sizes
  int nPoints = this->getNumberOfPoints(path);
  this->printMsg("nPoints: " + std::to_string(nPoints));

  // TODO: initialize buffers
  auto points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(nPoints);
  output->SetPoints(points);

  auto categoryIndex = vtkSmartPointer<vtkUnsignedCharArray>::New();
  categoryIndex->SetName("CategoryIndex");
  categoryIndex->SetNumberOfComponents(1);
  categoryIndex->SetNumberOfTuples(nPoints);

  auto yearIndex = vtkSmartPointer<vtkUnsignedCharArray>::New();
  yearIndex->SetName("YearIndex");
  yearIndex->SetNumberOfComponents(1);
  yearIndex->SetNumberOfTuples(nPoints);

  std::vector<string> categoryDictionary;
  std::vector<string> yearDictionary;

  // TODO: fill buffers
  int status = this->parsePointCoords(
    path, (float *)ttkUtils::GetVoidPointer(points),
    (unsigned char *)ttkUtils::GetVoidPointer(categoryIndex),
    (unsigned char *)ttkUtils::GetVoidPointer(yearIndex),
    categoryDictionary,
    yearDictionary);

  // finalize output
  {
    auto cellOffsets = vtkSmartPointer<vtkIdTypeArray>::New();
    cellOffsets->SetNumberOfTuples(nPoints + 1);
    auto cellOffsetsData = (vtkIdType *)ttkUtils::GetVoidPointer(cellOffsets);
    for(int i = 0; i < nPoints + 1; i++)
      cellOffsetsData[i] = i;

    auto cellConnectivity = vtkSmartPointer<vtkIdTypeArray>::New();
    cellConnectivity->SetNumberOfTuples(nPoints);
    auto cellConnectivityData
      = (vtkIdType *)ttkUtils::GetVoidPointer(cellConnectivity);
    for(int i = 0; i < nPoints; i++)
      cellConnectivityData[i] = i;

    auto cellArray = vtkSmartPointer<vtkCellArray>::New();
    cellArray->SetData(cellOffsets, cellConnectivity);
    output->SetCells(VTK_VERTEX, cellArray);

    // create filed data for the category 
    auto categoryDictionaryArray = vtkSmartPointer<vtkStringArray>::New();
    categoryDictionaryArray->SetName("CategoryDictionary");
    categoryDictionaryArray->SetNumberOfComponents(1);
    categoryDictionaryArray->SetNumberOfTuples(categoryDictionary.size());

    for(int i = 0; i < categoryDictionary.size(); i++)
      categoryDictionaryArray->SetValue(i, categoryDictionary[i]);

    // create filed data for the year of each event
    auto yearDictionaryArray = vtkSmartPointer<vtkStringArray>::New();
    yearDictionaryArray->SetName("YearDictionary");
    yearDictionaryArray->SetNumberOfComponents(1);
    yearDictionaryArray->SetNumberOfTuples(yearDictionary.size());

    for(int i = 0; i < yearDictionary.size(); i++)
      yearDictionaryArray->SetValue(i, yearDictionary[i]);
    
    // add property and filed d ata to the output

    output->GetPointData()->AddArray(categoryIndex);
    output->GetPointData()->AddArray(yearIndex);

    output->GetFieldData()->AddArray(categoryDictionaryArray);
    output->GetFieldData()->AddArray(yearDictionaryArray);

    output->Print(std::cout);
  }

  return 1;
}
