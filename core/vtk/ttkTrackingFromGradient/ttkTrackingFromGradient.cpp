#include <ttkTrackingFromGradient.h>

#include <vtkInformation.h>
#include <vtkDataSet.h>
#include <vtkMultiBlockDataSet.h>

#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkImageData.h>
#include <vtkIntArray.h>
#include <vtkSignedCharArray.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkTrackingFromGradient);

ttkTrackingFromGradient::ttkTrackingFromGradient() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

ttkTrackingFromGradient::~ttkTrackingFromGradient() {}

int ttkTrackingFromGradient::FillInputPortInformation(int port, vtkInformation *info) {
  if(port == 0 || port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    info->Append(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

int ttkTrackingFromGradient::FillOutputPortInformation(int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

vtkIntArray* GetVertexIdArray(vtkDataSet* input){
  return vtkIntArray::SafeDownCast( input->GetPointData()->GetArray("ttkVertexScalarField") );
}

int ttkTrackingFromGradient::RequestData(
  vtkInformation *request,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector
){
  auto output = vtkMultiBlockDataSet::GetData(outputVector);

  auto inputDomain = vtkDataObject::GetData(inputVector[0]);
  auto inputCriticalPoints = vtkDataObject::GetData(inputVector[1]);
  if(!inputDomain || !inputCriticalPoints)
    return !this->printErr("Unable to retrieve input data objects.");

  auto inputDomainAsMB = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  auto inputCriticalPointsAsMB = vtkSmartPointer<vtkMultiBlockDataSet>::New();

  auto inputDomainAsDS = vtkDataSet::SafeDownCast(inputDomain);
  auto inputCriticalPointsAsDS = vtkDataSet::SafeDownCast(inputCriticalPoints);

  if(inputDomainAsDS && inputCriticalPointsAsDS){
    inputDomainAsMB->SetBlock(0, this->PreviousDomain);
    inputDomainAsMB->SetBlock(1, inputDomainAsDS);

    inputCriticalPointsAsMB->SetBlock(0, this->PreviousCriticalPoints);
    inputCriticalPointsAsMB->SetBlock(1, inputCriticalPointsAsDS);

    // Get iteration information
    auto iterationInformation = vtkDoubleArray::SafeDownCast(
      inputDomainAsDS->GetFieldData()->GetArray("_ttk_IterationInfo")
    );

    // terminate if first element of streaming
    if(!this->PreviousDomain || (iterationInformation && iterationInformation->GetValue(0)<1)){
      this->printMsg("Initializing Streamed Tracking.");

      this->PreviousDomain = vtkSmartPointer<vtkDataSet>::Take( inputDomainAsDS->NewInstance() );
      this->PreviousDomain->ShallowCopy(inputDomainAsDS);

      this->PreviousCriticalPoints = vtkSmartPointer<vtkDataSet>::Take( inputCriticalPointsAsDS->NewInstance() );
      this->PreviousCriticalPoints->ShallowCopy(inputCriticalPointsAsDS);

      return 1;
    }
  } else if(inputDomain->IsA("vtkMultiBlockDataSet") && inputCriticalPoints->IsA("vtkMultiBlockDataSet")) {
    inputDomainAsMB->ShallowCopy(inputDomain);
    inputCriticalPointsAsMB->ShallowCopy(inputCriticalPoints);
  } else {
    return !this->printErr("Unsupported input data object types: "+std::string(inputDomain->GetClassName()) +" "+std::string(inputCriticalPoints->GetClassName()));
  }

  if(this->GetInputArrayAssociation(0, inputVector) != 0)
    return !this->printErr("Scalars must be a point data array.");

  const int nSteps = inputDomainAsMB->GetNumberOfBlocks();
  if(nSteps != inputCriticalPointsAsMB->GetNumberOfBlocks())
    return !this->printErr("Number of domains and critical point sets must be equal.");

  for(int b=1; b<nSteps; b++){

    auto domain0 = vtkDataSet::SafeDownCast(inputDomainAsMB->GetBlock(b-1));
    auto domain1 = vtkDataSet::SafeDownCast(inputDomainAsMB->GetBlock(b));
    auto criticalPoints0 = vtkDataSet::SafeDownCast(inputCriticalPointsAsMB->GetBlock(b-1));
    auto criticalPoints1 = vtkDataSet::SafeDownCast(inputCriticalPointsAsMB->GetBlock(b));

    int nFeatures0 = criticalPoints0->GetNumberOfElements(0);
    int nFeatures1 = criticalPoints1->GetNumberOfElements(0);

    // allocate correlation matrices
    auto correlations = vtkSmartPointer<vtkImageData>::New();
    correlations->SetDimensions(nFeatures0,nFeatures1,1);
    correlations->AllocateScalars(VTK_INT,1);

    auto forward = correlations->GetPointData()->GetArray(0);
    forward->SetName("Forward");

    auto backward = vtkSmartPointer<vtkIntArray>::New();
    backward->DeepCopy(forward);
    backward->SetName("Backward");
    correlations->GetPointData()->AddArray(backward);

    auto orderArray0 = ttkAlgorithm::GetOrderArray(domain0, 0);
    auto orderArray1 = ttkAlgorithm::GetOrderArray(domain1, 0);

    int status = 0;

    auto triangulation = this->GetTriangulation(domain0);

    this->preconditionTriangulation(triangulation);

    ttkVtkTemplateMacro(
      orderArray1->GetDataType(),
      triangulation->getType(),
      (
        status = this->computeCorrelations<ttk::SimplexId, TTK_TT>(
          ttkUtils::GetPointer<int>(forward),

          ttkUtils::GetPointer<ttk::SimplexId>(orderArray1),
          static_cast<TTK_TT*>(triangulation->getData()),

          ttkUtils::GetPointer<ttk::SimplexId>(GetVertexIdArray(criticalPoints0)),
          ttkUtils::GetPointer<ttk::SimplexId>(GetVertexIdArray(criticalPoints1)),
          nFeatures0,
          nFeatures1,
          std::greater<ttk::SimplexId>{},
          [](ttk::SimplexId i, ttk::SimplexId j, ttk::SimplexId n, ttk::SimplexId m){return j*n+i;}
        )
      )
    );
    if(!status)
      return 0;

    ttkVtkTemplateMacro(
      orderArray0->GetDataType(),
      triangulation->getType(),
      (
        status = this->computeCorrelations<ttk::SimplexId, TTK_TT>(
          ttkUtils::GetPointer<int>(backward),

          ttkUtils::GetPointer<ttk::SimplexId>(orderArray0),
          static_cast<TTK_TT*>(triangulation->getData()),

          ttkUtils::GetPointer<ttk::SimplexId>(GetVertexIdArray(criticalPoints1)),
          ttkUtils::GetPointer<ttk::SimplexId>(GetVertexIdArray(criticalPoints0)),
          nFeatures1,
          nFeatures0,
          std::greater<ttk::SimplexId>{},
          [](ttk::SimplexId i, ttk::SimplexId j, ttk::SimplexId n, ttk::SimplexId m){return i*m+j;}
        )
      )
    );
    if(!status)
      return 0;

    output->SetBlock(b-1, correlations);
  }

  // if streaming mode
  if(inputDomainAsDS && inputCriticalPointsAsDS){
    this->PreviousDomain = vtkSmartPointer<vtkDataSet>::Take( inputDomainAsDS->NewInstance() );
    this->PreviousDomain->ShallowCopy(inputDomainAsDS);

    this->PreviousCriticalPoints = vtkSmartPointer<vtkDataSet>::Take( inputCriticalPointsAsDS->NewInstance() );
    this->PreviousCriticalPoints->ShallowCopy(inputCriticalPointsAsDS);
  }

  return 1;
}