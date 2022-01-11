#include <ttkCorrespondenceByJacobiSet.h>

#include <vtkInformation.h>
#include <vtkObjectFactory.h>

#include <vtkFloatArray.h>
#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointSet.h>

#include <vtkPointData.h>
#include <vtkStringArray.h>

#include <ttkUtils.h>
#include <ttkMacros.h>

#include <vtkUnstructuredGrid.h>

// TTK Base Includes
#include <ttkJacobiSet.h>
#include <ttkPointMerger.h>
#include <ttkConnectedComponents.h>

vtkStandardNewMacro(ttkCorrespondenceByJacobiSet);

ttkCorrespondenceByJacobiSet::ttkCorrespondenceByJacobiSet() {
  this->setDebugMsgPrefix("CorrespondenceByJacobiSet");

  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

ttkCorrespondenceByJacobiSet::~ttkCorrespondenceByJacobiSet() {
}

template<typename DT>
int computeStackedArray(
  DT* sd,
  const DT* d0,
  const DT* d1,
  const int n
){
  for(int i=0; i<n; i++){
    sd[i] = d0[i];
  }
  for(int i=0, j=n; i<n; i++,j++){
    sd[j] = d1[i];
  }
  return 1;
}

int ttkCorrespondenceByJacobiSet::ComputeCorrespondences(
  vtkImageData *correspondenceMatrix,
  vtkDataObject *inputDataObjects0,
  vtkDataObject *inputDataObjects1) {

  // unpack input
  auto inputsAsMB0 = vtkMultiBlockDataSet::SafeDownCast(inputDataObjects0);
  auto inputsAsMB1 = vtkMultiBlockDataSet::SafeDownCast(inputDataObjects1);
  if(!inputsAsMB0 || !inputsAsMB1)
    return !this->printErr("Unable to retrieve input data objects.");

  auto image0 = vtkImageData::SafeDownCast(inputsAsMB0->GetBlock(0));
  auto image1 = vtkImageData::SafeDownCast(inputsAsMB1->GetBlock(0));
  if(!image0 || !image1)
    return !this->printErr("Unable to retrieve input grid data objects.");

  auto points0 = vtkPointSet::SafeDownCast(inputsAsMB0->GetBlock(1));
  auto points1 = vtkPointSet::SafeDownCast(inputsAsMB1->GetBlock(1));
  if(!points0 || !points1)
    return !this->printErr("Unable to retrieve input critical point data objects.");

  // check if input images are two dimensional
  int dim[3];
  {
    int dim_[3];
    image0->GetDimensions(dim);
    image1->GetDimensions(dim_);
    if(dim[0]!=dim_[0] || dim[1]!=dim_[1] || dim[2]!=dim_[2] || dim[2]!=1)
      return !this->printErr("Input grids need to be two dimensional and must have same dimension.");
  }

  // // retrieve scalar and order arrays
  // auto scalar0 = this->GetInputArrayToProcess(0,image0);
  // auto scalar1 = this->GetInputArrayToProcess(0,image1);
  // if(!scalar0 || !scalar1)
  //   return !this->printErr("Unable to retrieve scalar arrays.");

  // const int nTuplesPerImage = scalar0->GetNumberOfTuples();

  // auto stackedScalarArray = vtkSmartPointer<vtkDataArray>::Take(scalar0->NewInstance());
  // stackedScalarArray->SetName("Scalars");
  // stackedScalarArray->SetNumberOfTuples(2*nTuplesPerImage);
  // ttkTypeMacroA(
  //   stackedScalarArray->GetDataType(),
  //   computeStackedArray<T0>(
  //     ttkUtils::GetPointer<T0>(stackedScalarArray),
  //     ttkUtils::GetPointer<T0>(scalar0),
  //     ttkUtils::GetPointer<T0>(scalar1),
  //     nTuplesPerImage
  //   )
  // );

  // retrieve order arrays
  auto order0 = this->GetOrderArray(image0, 0);
  auto order1 = this->GetOrderArray(image1, 0);
  if(!order0 || !order1)
    return !this->printErr("Unable to retrieve order arrays.");
  const int nTuplesPerImage = order0->GetNumberOfTuples();

  // initialize data arrays of stacked image data object
  auto stackedOrderArray = vtkSmartPointer<vtkDataArray>::Take(order0->NewInstance());
  stackedOrderArray->SetName("ORDER");
  stackedOrderArray->SetNumberOfTuples(nTuplesPerImage*2);
  auto stackedOrderArrayData = ttkUtils::GetPointer<ttk::SimplexId>(stackedOrderArray);

  auto timeArray = vtkSmartPointer<vtkDataArray>::Take(order0->NewInstance());
  timeArray->SetName("TIME");
  timeArray->SetNumberOfTuples(nTuplesPerImage*2);
  auto timeArrayData = ttkUtils::GetPointer<ttk::SimplexId>(timeArray);
  {
    auto order0Data = ttkUtils::GetPointer<ttk::SimplexId>(order0);
    auto order1Data = ttkUtils::GetPointer<ttk::SimplexId>(order1);
    for(int i=0; i<nTuplesPerImage; i++){
      timeArrayData[i] = 0;
      stackedOrderArrayData[i] = order0Data[i];
    }
    for(int i=0,j=nTuplesPerImage; i<nTuplesPerImage; i++,j++){
      timeArrayData[j] = 1;
      stackedOrderArrayData[j] = order1Data[i];
    }
  }

  vtkSmartPointer<vtkDataArray> stackedVertexIdentifiers;
  {
    // stack vertex identifiers
    auto vertexIdentifiers0 = this->GetInputArrayToProcess(1,image0);
    auto vertexIdentifiers1 = this->GetInputArrayToProcess(1,image1);
    if(!vertexIdentifiers0 || !vertexIdentifiers1)
      return !this->printErr("Unable to retrieve vertex identifier arrays from input grid.");

    stackedVertexIdentifiers = vtkSmartPointer<vtkDataArray>::Take(vertexIdentifiers1->NewInstance());
    stackedVertexIdentifiers->SetName(vertexIdentifiers0->GetName());
    stackedVertexIdentifiers->SetNumberOfTuples(nTuplesPerImage*2);
    ttkTypeMacroA(
      stackedVertexIdentifiers->GetDataType(),
      computeStackedArray<T0>(
        ttkUtils::GetPointer<T0>(stackedVertexIdentifiers),
        ttkUtils::GetPointer<T0>(vertexIdentifiers0),
        ttkUtils::GetPointer<T0>(vertexIdentifiers1),
        nTuplesPerImage
      )
    );
  }

  // build stacked image data object
  auto stackedImage = vtkSmartPointer<vtkImageData>::New();
  stackedImage->SetDimensions(dim[0], dim[1], 2);

  auto stackedImagePD = stackedImage->GetPointData();
  stackedImagePD->AddArray(stackedOrderArray);
  stackedImagePD->AddArray(timeArray);
  stackedImagePD->AddArray(stackedVertexIdentifiers);


  // compute jacobi set on stacked image
  auto jacobiSetFilter = vtkSmartPointer<ttkJacobiSet>::New();
  jacobiSetFilter->SetInputDataObject(stackedImage);
  jacobiSetFilter->SetDebugLevel(this->debugLevel_);
  jacobiSetFilter->SetVertexScalars(true);
  jacobiSetFilter->SetInputArrayToProcess(0,0,0,0,"ORDER");
  jacobiSetFilter->SetInputArrayToProcess(1,0,0,0,"TIME");
  jacobiSetFilter->Update();

  auto jacobiSet = vtkUnstructuredGrid::SafeDownCast(jacobiSetFilter->GetOutputDataObject(0));

  // count number of temporal cells
  int nTemporalCells = 0;
  {
    auto pointCoords = ttkUtils::GetPointer<float>(jacobiSet->GetPoints()->GetData());

    int nCells = jacobiSet->GetNumberOfCells();
    for(int i=0; i<nCells; i++){
      auto pointIds = jacobiSet->GetCell(i)->GetPointIds();
      if(pointCoords[pointIds->GetId(0)*3+2]!=pointCoords[pointIds->GetId(1)*3+2])
        nTemporalCells++;
    }
  }

  this->printMsg("Number of Temporal Jacobi Edges: " + std::to_string(nTemporalCells));

  // build new edge set
  {
    ttk::Timer t;
    const std::string msg = "Extracting Temporal Edges";
    this->printMsg(msg,0,0,1,ttk::debug::LineMode::REPLACE);

    auto temporalEdges = vtkSmartPointer<vtkUnstructuredGrid>::New();
    temporalEdges->AllocateExact(nTemporalCells,nTemporalCells*2);
    auto pointCoords = ttkUtils::GetPointer<float>(jacobiSet->GetPoints()->GetData());

    auto points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(nTemporalCells*2);

    int nCells = jacobiSet->GetNumberOfCells();
    for(int i=0,j=0; i<nCells; i++){
      const auto pointIds = jacobiSet->GetCell(i)->GetPointIds();
      const auto u = pointIds->GetId(0);
      const auto v = pointIds->GetId(1);

      if(pointCoords[u*3+2]!=pointCoords[v*3+2]){
        points->SetPoint(j, &pointCoords[u*3]);
        points->SetPoint(j+1, &pointCoords[v*3]);
        vtkIdType ids[2] = {j,j+1};
        temporalEdges->InsertNextCell(VTK_LINE,2,ids);
        j+=2;
      }
    }
    temporalEdges->SetPoints(points);

    this->printMsg(msg,0.33,t.getElapsedTime(),1,ttk::debug::LineMode::REPLACE);

    auto pointMerger = vtkSmartPointer<ttkPointMerger>::New();
    pointMerger->SetInputDataObject(temporalEdges);
    pointMerger->SetBoundaryOnly(false);
    pointMerger->SetDistanceThreshold(0.01);

    auto connectedComponents = vtkSmartPointer<ttkConnectedComponents>::New();
    connectedComponents->SetInputConnection(pointMerger->GetOutputPort());
    connectedComponents->SetInputArrayToProcess(0,0,0,0,"NONE");
    connectedComponents->SetUseSeedIdAsComponentId(false);
    connectedComponents->Update();

    this->printMsg(msg,1,t.getElapsedTime(),1);

    auto components = vtkUnstructuredGrid::SafeDownCast(connectedComponents->GetOutputDataObject(0));
    if(!components)
      return !this->printErr("Unable to merge points and compute connected components of temporal Jacobi edges.");

    auto componentIds = ttkUtils::GetPointer<int>(components->GetPointData()->GetArray("ComponentId"));
    if(!componentIds)
      return !this->printErr("Unable to retrieve componentIds from temporal Jacobi edges");
  }

  // iterate over features and find for each critical points its corresponding connected component






  correspondenceMatrix->ShallowCopy(stackedImage);

  // this->setSosOffsetsU(combinedOrderArrayData);
  // this->setSosOffsetsV(stackedOrderArrayData);
  // auto triangulation = ttkAlgorithm::GetTriangulation(stackedGrid);
  // if(!triangulation)
  //   return !this->printErr("Unable to derive triangulation of stacked image data object.");
  // this->preconditionTriangulation(triangulation);

  // std::vector<std::pair<ttk::SimplexId, char>> jacobiSet{};
  // // std::vector<char> isPareto{};


  // int status = this->execute<ttk::SimplexId,ttk::SimplexId,ttk::ImplicitTriangulation>(
  //   jacobiSet,
  //   combinedOrderArrayData,
  //   timeArrayData,
  //   *static_cast<ttk::ImplicitTriangulation*>(triangulation->getData())
  //   // &isPareto
  // );

  // if(status!=0)
  //   return 0;



  // vtkNew<vtkPolyData> output;
  // output->AllocateExact(0,0, jacobiSet.size(), jacobiSet.size()*2, 0,0, 0,0);

  // vtkNew<vtkPoints> pts;
  // output->SetPoints(pts);
  // vtkIdType q=0;
  // for(auto p: jacobiSet){
  //   for(int i=0; i<2; i++){
  //     ttk::SimplexId v;
  //     triangulation->getEdgeVertex(p.first,i,v);
  //     float pos[3];
  //     triangulation->getVertexPoint(v,pos[0],pos[1],pos[2]);
  //     pts->InsertNextPoint(pos);
  //   }
  //   const vtkIdType ids[2]{q,q+1};
  //   output->InsertNextCell(VTK_LINE, 2, ids);
  //   q+=2;
  // }

  // vtkNew<vtkPolyDataWriter> writer;
  // writer->SetInputDataObject(output);
  // writer->SetFileName("/home/jones/external/data/ttk-data/test.vtp");
  // writer->Update();

  // pts->InsertNextPoint(origin);
  // pts->InsertNextPoint(p0);
  // pts->InsertNextPoint(p1);
  // output->AllocateExact(
  //   0, // vtkIdType numVerts,
  //   0, // vtkIdType vertConnSize,
  //   jacobiSet.size(), // vtkIdType numLines,
  //   jacobiSet.size(), // vtkIdType lineConnSize,
  //   // vtkIdType numPolys,
  //   // vtkIdType polyConnSize,
  //   // vtkIdType numStrips,
  //   // vtkIdType stripConnSize
  // )


  // const auto uComponent = this->GetInputArrayToProcess(0, input);
  // const auto vComponent = this->GetInputArrayToProcess(1, input);


  // const int nPoints0 = p0->GetNumberOfPoints();
  // const int nPoints1 = p1->GetNumberOfPoints();

  // // get point coordinates

  // auto temp = vtkSmartPointer<vtkFloatArray>::New();
  // auto coords0 = nPoints0 > 0 ? p0->GetPoints()->GetData() : temp;
  // auto coords1 = nPoints1 > 0 ? p1->GetPoints()->GetData() : temp;

  // if(coords0->GetDataType() != coords1->GetDataType())
  //   return !this->printErr("Input vtkPointSet need to have same precision.");

  // // initialize correspondence matrix i.e., distance matrix
  // correspondenceMatrix->SetDimensions(nPoints0, nPoints1, 1);
  // correspondenceMatrix->AllocateScalars(coords0->GetDataType(), 1);
  // auto matrixData = correspondenceMatrix->GetPointData()->GetArray(0);
  // matrixData->SetName("Distance");

  // // compute distance matrix
  // int status = 0;
  // switch(coords0->GetDataType()) {
  //   vtkTemplateMacro(status = this->computeDistanceMatrix<VTK_TT>(
  //                     ttkUtils::GetPointer<VTK_TT>(matrixData),
  //                     ttkUtils::GetPointer<const VTK_TT>(coords0),
  //                     ttkUtils::GetPointer<const VTK_TT>(coords1), nPoints0,
  //                     nPoints1));
  // }
  // if(!status)
  //   return 0;

  // status = ttkCorrespondenceAlgorithm::AddIndexLabelMaps(
  //   correspondenceMatrix, this->GetInputArrayToProcess(0, p0),
  //   this->GetInputArrayToProcess(0, p1));
  // if(!status)
  //   return 0;

  return 1;
}
