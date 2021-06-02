#include <ttkCorrespondenceByPersistencePairs.h>

#include <vtkObjectFactory.h>
#include <vtkInformation.h>
#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointSet.h>
#include <vtkPointData.h>
#include <vtkIntArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>

#include <ttkUtils.h>
#include <ttkMacros.h>

vtkStandardNewMacro(ttkCorrespondenceByPersistencePairs);

ttkCorrespondenceByPersistencePairs::ttkCorrespondenceByPersistencePairs() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkCorrespondenceByPersistencePairs::~ttkCorrespondenceByPersistencePairs() {
}

int ttkCorrespondenceByPersistencePairs::Correlate(
    vtkImageData *correspondences,
    vtkDataObject* inputDataObjects0,
    vtkDataObject* inputDataObjects1
){
  int status = 0;
    
    // unpack input
  vtkUnstructuredGrid *p0 = vtkUnstructuredGrid::SafeDownCast(inputDataObjects0); // set of persist. pairs from previous timestep
  vtkUnstructuredGrid *p1 = vtkUnstructuredGrid::SafeDownCast(inputDataObjects1); // PPs from current timestep
  if(!p0 || !p1)
    return !this->printErr("Input data objects need to be vtkPointSets.");

  // get point coordinates
  auto coords0 = p0->GetPoints()->GetData();
  auto coords1 = p1->GetPoints()->GetData();

  // get birth & death scalars
  //auto birth0 = this->GetInputArrayToProcess(0, p0);
  //auto death0 = this->GetInputArrayToProcess(1, p0);
  //auto birth1 = this->GetInputArrayToProcess(0, p1);
  //auto death1 = this->GetInputArrayToProcess(1, p1);

  // get and format persistence diagrams
  using dataType = double;
  using dT = std::tuple<int, ttk::CriticalType, int, ttk::CriticalType, dataType,
      int, dataType, float, float, float, dataType, float, float, float>;
  using mT = std::tuple<int, int, double>;

  std::vector<dT> CTDiagram0;
  std::vector<dT> CTDiagram1;

  const double spacing = 0; // 2d tracking -> height spacing parameter for 3d display
  status = getDiagram<double>(CTDiagram0, p0, spacing, 0);
  if (status < 0) {
    this->printErr("Could not extract diagram from first input data-set");
    return 0;
  }
  status = getDiagram<double>(CTDiagram1, p1, spacing, 1);
  if(status < 0) {
    this->printErr("Could not extract diagram from second input data-set");
    return 0;
  }
  if(coords0->GetDataType() != coords1->GetDataType())
    return !this->printErr("Input diagrams need to have the same data type.");

  // compute correspondence matrix dimensions
  //const int nPairs0 = coords0->GetNumberOfTuples() / 2;
  //const int nPairs1 = coords1->GetNumberOfTuples() / 2;
  const int nFeatures0 = CTDiagram0.size();
  const int nFeatures1 = CTDiagram1.size();

  // initialize correspondence matrix i.e., distance matrix
  correspondences->SetDimensions(nFeatures0, nFeatures1, 1);
  correspondences->AllocateScalars(VTK_DOUBLE, 1); // matching output = double

  auto correspondencesArray = correspondences->GetPointData()->GetArray(0);
  correspondencesArray->SetName("LiftedWassersteinDistance");
  auto correspondanceMatrix = ttkUtils::GetPointer<float>(correspondencesArray);
  for (int i = 0; i < nFeatures0; ++i)
      for (int j = 0; j < nFeatures1; ++j) 
      {
          correspondanceMatrix[j * nFeatures0 + i] = 0;
      }

  // get metric parameters
  const std::string algorithm = DistanceAlgorithm;
  const std::string wasserstein = WassersteinMetric;
  const double alpha = Alpha;
  const int pvAlgorithm = PVAlgorithm;
  const double px = PX;
  const double py = PY;
  const double pz = PZ;
  const double ps = PS;
  const double pe = PE;
  
  // compute correspondences in basecode
  //status = 0;
  std::vector<mT> matchings;
  switch (coords0->GetDataType()) {
  vtkTemplateMacro(
      status = this->computeDistanceMatrix<double>(
      //ttkUtils::GetPointer<float>(correspondencesArray),
      //ttkUtils::GetPointer<const VTK_TT>(coords0),
      //ttkUtils::GetPointer<const VTK_TT>(coords1),
      CTDiagram0,
      CTDiagram1,
      matchings,
      //nFeatures0,
      //nFeatures1,
      px, py, pz, ps, pe, algorithm, wasserstein, alpha, pvAlgorithm));
  }
  if (status < 0) 
  {
      this->printErr("Error computing distance matrix.");
      return -1;
  }
    //));
    //}

  // build matrix
  auto matchingsSize = matchings.size();
  if (matchingsSize > 0)
  {
    for (int i = 0; i < matchingsSize; ++i) {
      vtkIdType ids[2];
      mT t = matchings.at((unsigned long)i);
      auto n1 = (int)std::get<0>(t); // diagram 0
      auto n2 = (int)std::get<1>(t); // diagram 1
      if (n1 >= nFeatures0 || n2 >= nFeatures1) 
      {
          this->printErr("Invalid indexing: feature index > feature number.");
          continue;
      }

      correspondanceMatrix[n2 * nFeatures0 + n1]++;
    }
  }
  
  // add index label maps
  int a=0;
  auto fd = correspondences->GetFieldData();
  for (auto& it : std::vector<vtkUnstructuredGrid*> ({p0, p1})) // vtkPointSett
  {
    auto labels = this->GetInputArrayToProcess(0, it);
    if(!labels)
      return !this->printErr("Unable to retrieve labels.");
    auto array = vtkSmartPointer<vtkDataArray>::Take( labels->NewInstance() );
    array->ShallowCopy(labels);
    array->SetName(("IndexLabelMap"+std::to_string(a++)).data());
    fd->AddArray(array);
  }

  return 1;
}

template <typename dataType>
int ttkCorrespondenceByPersistencePairs::getDiagram(
  std::vector<std::tuple<int, ttk::CriticalType, int, ttk::CriticalType, dataType,
      int, dataType, float, float, float, dataType, float, float, float> > &diagram,
  vtkUnstructuredGrid *CTPersistenceDiagram_,
  const double spacing,
  const int diagramNumber)
{
  auto pointData = CTPersistenceDiagram_->GetPointData();
  auto cellData = CTPersistenceDiagram_->GetCellData();

  if(pointData == nullptr || cellData == nullptr) {
    return -1;
  }

  auto vertexIdentifierScalars = ttkSimplexIdTypeArray::SafeDownCast(
    pointData->GetArray(ttk::VertexScalarFieldName));

  auto nodeTypeScalars = vtkIntArray::SafeDownCast(pointData->GetArray("CriticalType"));

  auto pairIdentifierScalars = ttkSimplexIdTypeArray::SafeDownCast(cellData->GetArray("PairIdentifier"));

  auto extremumIndexScalars = vtkIntArray::SafeDownCast(cellData->GetArray("PairType"));

  auto persistenceScalars = vtkDoubleArray::SafeDownCast(cellData->GetArray("Persistence"));

  auto birthScalars = vtkDoubleArray::SafeDownCast(pointData->GetArray("Birth"));

  auto deathScalars = vtkDoubleArray::SafeDownCast(pointData->GetArray("Death"));

  vtkPoints *points = CTPersistenceDiagram_->GetPoints();
  if (!pairIdentifierScalars) return -2;

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
  //if(Is3D && !is3D)
  //  Is3D = false;
  //if(!Is3D && diagramNumber == 1)
  //  s = (float)spacing;

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
      diagram.at(pairIdentifier)
        = std::make_tuple(vertexId1, (BNodeType)nodeType1, vertexId2,
                          (BNodeType)nodeType2, (dataType)persistence, pairType,
                          value1, x1, y1, z1 + s, value2, x2, y2, z2 + s);

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

  sort(diagram.begin(), diagram.end(),
       [](const diagramTuple &a, const diagramTuple &b) -> bool {
         return std::get<6>(a) < std::get<6>(b);
       });

  return 1;
}
