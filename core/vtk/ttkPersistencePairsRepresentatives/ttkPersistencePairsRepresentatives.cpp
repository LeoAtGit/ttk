#include <ttkPersistencePairsRepresentatives.h>

#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkPersistencePairsRepresentatives);

ttkPersistencePairsRepresentatives::ttkPersistencePairsRepresentatives() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkPersistencePairsRepresentatives::~ttkPersistencePairsRepresentatives() {
}

int ttkPersistencePairsRepresentatives::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

int ttkPersistencePairsRepresentatives::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");
    return 1;
  }
  return 0;
}

int ttkPersistencePairsRepresentatives::RequestData(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {
  // Fetch Input Data
  auto inputDataSet = vtkUnstructuredGrid::GetData(inputVector[0]);
  if(!inputDataSet)
    return 0;

  // get persistence diagram
  auto coords0 = inputDataSet->GetPoints()->GetData();

  // get and format persistence diagrams
  using dataType = double;
  using dT
    = std::tuple<int, ttk::CriticalType, int, ttk::CriticalType, dataType, int,
                 dataType, float, float, float, dataType, float, float, float>;
  using mT = std::tuple<int, int, double>;
  std::vector<dT> inputDiagram;
  using NodeType = ttk::CriticalType;
  auto LocalMax = ttk::CriticalType::Local_maximum;
  auto LocalMin = ttk::CriticalType::Local_minimum;

  const double spacing
    = 0; // 2d tracking -> height spacing parameter for 3d display
  auto status = getDiagram<double>(inputDiagram, inputDataSet, spacing, 0);
  if(status < 0) {
    this->printErr("Could not extract diagram from first input data-set");
    return 0;
  }

  // PD Output
  {
    const int nComponents = inputDiagram.size();
    auto outputComponents = vtkPolyData::GetData(outputVector, 0); // 1);

    // points
    {
      auto persistenceArray = vtkSmartPointer<vtkFloatArray>::New();
      persistenceArray->SetName("Persistence");
      persistenceArray->SetNumberOfTuples(nComponents);
      auto persistenceArrayData = ttkUtils::GetPointer<float>(persistenceArray);

      auto valueArray = vtkSmartPointer<vtkFloatArray>::New();
      valueArray->SetName("Value");
      valueArray->SetNumberOfTuples(nComponents);
      auto valueArrayData = ttkUtils::GetPointer<float>(valueArray);

      auto idArray = vtkSmartPointer<vtkIntArray>::New();
      idArray->SetName("PairId");
      idArray->SetNumberOfTuples(nComponents);
      auto idArrayData = ttkUtils::GetPointer<int>(idArray);

      auto points = vtkSmartPointer<vtkPoints>::New();
      points->SetDataTypeToFloat();
      points->SetNumberOfPoints(nComponents);
      auto pointsData = ttkUtils::GetPointer<float>(points->GetData());
      for(int i = 0, j = 0; i < nComponents; i++) {
        const auto &t = inputDiagram[i];

        NodeType type1 = std::get<1>(t);
        NodeType type2 = std::get<3>(t);
        bool t11Min = type1 == LocalMin;
        bool t11Max = type1 == LocalMax;
        bool t12Min = type2 == LocalMin;
        bool t12Max = type2 == LocalMax;
        bool t1Max = t11Max || t12Max;
        bool t1Min = !t1Max && (t11Min || t12Min);

        float x1 = (float)(t1Max   ? std::get<11>(t)
                           : t1Min ? std::get<7>(t)
                                   : 0);
        float y1 = (float)(t1Max   ? std::get<12>(t)
                           : t1Min ? std::get<8>(t)
                                   : 0);
        float z1 = (float)(t1Max   ? std::get<13>(t)
                           : t1Min ? std::get<9>(t)
                                   : 0);

        float v1 = (float)(t1Max   ? std::get<10>(t)
                           : t1Min ? std::get<6>(t)
                                   : 0);

        pointsData[j++] = x1;
        pointsData[j++] = y1;
        pointsData[j++] = z1;

        persistenceArrayData[i] = std::get<4>(t);
        valueArrayData[i] = v1;
        idArrayData[i] = i;
      }

      outputComponents->SetPoints(points);
      auto pd = outputComponents->GetPointData();
      pd->AddArray(persistenceArray);
      pd->AddArray(valueArray);
      pd->AddArray(idArray);
    }

    // cells
    {
      auto connectivityArray = vtkSmartPointer<vtkIntArray>::New();
      connectivityArray->SetNumberOfTuples(nComponents);
      auto connectivityArrayData = ttkUtils::GetPointer<int>(connectivityArray);
      for(int i = 0; i < nComponents; i++)
        connectivityArrayData[i] = i;

      auto offsetArray = vtkSmartPointer<vtkIntArray>::New();
      offsetArray->SetNumberOfTuples(nComponents + 1);
      auto offsetArrayData = ttkUtils::GetPointer<int>(offsetArray);
      for(int i = 0; i <= nComponents; i++)
        offsetArrayData[i] = i;

      auto cellArray = vtkSmartPointer<vtkCellArray>::New();
      cellArray->SetData(offsetArray, connectivityArray);

      outputComponents->SetVerts(cellArray);
    }
  }

  // return success
  return 1;
}

template <typename dataType>
int ttkPersistencePairsRepresentatives::getDiagram(
  std::vector<std::tuple<int,
                         ttk::CriticalType,
                         int,
                         ttk::CriticalType,
                         dataType,
                         int,
                         dataType,
                         float,
                         float,
                         float,
                         dataType,
                         float,
                         float,
                         float>> &diagram,
  vtkUnstructuredGrid *CTPersistenceDiagram_,
  const double spacing,
  const int diagramNumber) {
  auto pointData = CTPersistenceDiagram_->GetPointData();
  auto cellData = CTPersistenceDiagram_->GetCellData();

  if(pointData == nullptr || cellData == nullptr) {
    return -1;
  }

  auto vertexIdentifierScalars = ttkSimplexIdTypeArray::SafeDownCast(
    pointData->GetArray(ttk::VertexScalarFieldName));
  auto nodeTypeScalars
    = vtkIntArray::SafeDownCast(pointData->GetArray("CriticalType"));
  auto pairIdentifierScalars
    = ttkSimplexIdTypeArray::SafeDownCast(cellData->GetArray("PairIdentifier"));
  auto extremumIndexScalars
    = vtkIntArray::SafeDownCast(cellData->GetArray("PairType"));
  auto persistenceScalars
    = vtkDoubleArray::SafeDownCast(cellData->GetArray("Persistence"));
  auto birthScalars
    = vtkDoubleArray::SafeDownCast(pointData->GetArray("Birth"));
  auto deathScalars
    = vtkDoubleArray::SafeDownCast(pointData->GetArray("Death"));

  vtkPoints *points = CTPersistenceDiagram_->GetPoints();
  if(!pairIdentifierScalars)
    return -2;

  auto pairingsSize = (int)pairIdentifierScalars->GetNumberOfTuples();

  // Continuous indexing (no gap in indices)
  for(int pairIndex = 0; pairIndex < pairingsSize; ++pairIndex) {
    const float indexOfPair = pairIndex;
    if(*pairIdentifierScalars->GetTuple(pairIndex) != -1) // except diagonal
      pairIdentifierScalars->SetTuple(pairIndex, &indexOfPair);
  }

  float s{0.0};

  if(!deathScalars != !birthScalars)
    return -2;
  bool is2D = !deathScalars && !birthScalars;
  bool is3D = !is2D;

  if(pairingsSize < 1 || !vertexIdentifierScalars || !nodeTypeScalars
     || !persistenceScalars || !extremumIndexScalars || !points)
    return -2;

  diagram.resize((unsigned long)pairingsSize);
  int nbNonCompact = 0;

  for(int i = 0; i < pairingsSize; ++i) {

    int vertexId1 = vertexIdentifierScalars->GetValue(2 * i);
    int vertexId2 = vertexIdentifierScalars->GetValue(2 * i + 1);
    int nodeType1 = nodeTypeScalars->GetValue(2 * i);
    int nodeType2 = nodeTypeScalars->GetValue(2 * i + 1);

    int pairIdentifier = pairIdentifierScalars->GetValue(i);
    int pairType = extremumIndexScalars->GetValue(i);
    double persistence = persistenceScalars->GetValue(i);

    int index1 = 2 * i;
    double *coords1 = points->GetPoint(index1);
    auto x1 = (float)coords1[0];
    auto y1 = (float)coords1[1];
    auto z1 = (float)coords1[2];

    int index2 = index1 + 1;
    double *coords2 = points->GetPoint(index2);
    auto x2 = (float)coords2[0];
    auto y2 = (float)coords2[1];
    auto z2 = (float)coords2[2];

    dataType value1 = (!birthScalars) ? (dataType)x1
                                      : (dataType)birthScalars->GetValue(2 * i);
    dataType value2 = (!deathScalars)
                        ? (dataType)y2
                        : (dataType)deathScalars->GetValue(2 * i + 1);

    if(pairIdentifier != -1 && pairIdentifier < pairingsSize)
      diagram.at(pairIdentifier) = std::make_tuple(
        vertexId1, (ttk::CriticalType)nodeType1, vertexId2,
        (ttk::CriticalType)nodeType2, (dataType)persistence, pairType, value1,
        x1, y1, z1 + s, value2, x2, y2, z2 + s);

    if(pairIdentifier >= pairingsSize) {
      nbNonCompact++;
      if(nbNonCompact == 0) {
        std::stringstream msg;
        msg << "Diagram pair identifiers "
            << "must be compact (not exceed the diagram size). " << std::endl;
        this->printWrn(msg.str());
      }
    }
  }

  if(nbNonCompact > 0) {
    std::stringstream msg;
    msg << "Missed " << nbNonCompact << " pairs due to non-compactness."
        << std::endl;
    this->printWrn(msg.str());
  }

  using diagramTuple
    = std::tuple<int, ttk::CriticalType, int, ttk::CriticalType, dataType, int,
                 dataType, float, float, float, dataType, float, float, float>;
  sort(diagram.begin(), diagram.end(),
       [](const diagramTuple &a, const diagramTuple &b) -> bool {
         return std::get<6>(a) < std::get<6>(b);
       });

  return 1;
}
