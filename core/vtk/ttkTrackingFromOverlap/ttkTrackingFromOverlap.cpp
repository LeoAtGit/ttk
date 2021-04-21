#include <ttkTrackingFromOverlap.h>

#include <vtkInformation.h>

#include <vtkMultiBlockDataSet.h>
#include <vtkImageData.h>
#include <vtkSmartPointer.h>

#include <vtkDoubleArray.h>
#include <vtkPointData.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkTrackingFromOverlap);

ttkTrackingFromOverlap::ttkTrackingFromOverlap() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkTrackingFromOverlap::~ttkTrackingFromOverlap() {
}

int ttkTrackingFromOverlap::FillInputPortInformation(int port,
                                                    vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    info->Append(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

int ttkTrackingFromOverlap::FillOutputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
    return 1;
  }
  return 0;
}

int GetNumberOfLabels(vtkDataArray* labels){
  if(labels->GetNumberOfTuples()<1)
    return 0;
  double labelRange[2];
  labels->GetRange(labelRange);
  return labelRange[1]>=0 ? labelRange[1]+1 : 0;
}

int ttkTrackingFromOverlap::Correlate(vtkImageData* correlations, vtkDataArray* labels0, vtkDataArray* labels1){
  // validate arrays
  if(!labels0 || !labels1)
    return !this->printErr("Unable to retrieve labels.");

  if(labels0->GetNumberOfComponents()!=1 || labels1->GetNumberOfComponents()!=1)
    return !this->printErr("Labels must have exactly one component.");

  if(labels0->GetNumberOfTuples()!=labels1->GetNumberOfTuples())
    return !this->printErr("Labels must have same number of values.");

  if(labels0->GetDataType()!=labels1->GetDataType())
    return !this->printErr("Labels must have same data type.");

  const int nVertices = labels0->GetNumberOfTuples();
  const int nLabels0 = GetNumberOfLabels(labels0);
  const int nLabels1 = GetNumberOfLabels(labels1);

  correlations->SetDimensions(nLabels0,nLabels1,1);
  correlations->AllocateScalars(VTK_INT,1);
  auto correlationsArray = correlations->GetPointData()->GetArray(0);
  correlationsArray->SetName("Overlap");

  int status = 0;
  switch(labels0->GetDataType()) {
    vtkTemplateMacro(
      status = this->computeAdjacencyMatrix<VTK_TT>(
        ttkUtils::GetPointer<int>(correlationsArray),
        ttkUtils::GetPointer<const VTK_TT>(labels0),
        ttkUtils::GetPointer<const VTK_TT>(labels1),
        nVertices,
        nLabels0,
        nLabels1));
  }
  if(!status)
    return 0;

  return 1;
}

int ttkTrackingFromOverlap::RequestData(vtkInformation *request,
                                       vtkInformationVector **inputVector,
                                       vtkInformationVector *outputVector) {

  auto output = vtkMultiBlockDataSet::GetData(outputVector);

  auto input = vtkDataObject::GetData(inputVector[0]);
  if(!input)
    return !this->printErr("Unable to retrieve input data object.");

  auto inputAsDS = vtkDataSet::SafeDownCast(input);
  auto inputAsMB = vtkSmartPointer<vtkMultiBlockDataSet>::New();
  if(inputAsDS){

    inputAsMB->SetBlock(0, this->PreviousInput);
    inputAsMB->SetBlock(1, inputAsDS);

    // Get iteration information
    auto iterationInformation = vtkDoubleArray::SafeDownCast(
      inputAsDS->GetFieldData()->GetArray("_ttk_IterationInfo")
    );

    // terminate if first element of streaming
    if(!this->PreviousInput || (iterationInformation && iterationInformation->GetValue(0)<1)){
      this->printMsg("Initializing Tracking with first vtkDataSet.");
      this->PreviousInput = vtkSmartPointer<vtkDataSet>::Take( inputAsDS->NewInstance() );
      this->PreviousInput->ShallowCopy(inputAsDS);

      return 1;
    }
  } else if(input->IsA("vtkMultiBlockDataSet")) {
    inputAsMB->ShallowCopy(input);
  } else {
    return !this->printErr("Unsupported input data object type: "+std::string(input->GetClassName()));
  }

  if(this->GetInputArrayAssociation(0, inputVector) != 0)
    return !this->printErr("Feature labels must be point data.");

  const int nSteps = inputAsMB->GetNumberOfBlocks();
  for(int b=1; b<nSteps; b++){
    auto correlations = vtkSmartPointer<vtkImageData>::New();

    if(!this->Correlate(
      correlations,
      this->GetInputArrayToProcess(0, inputAsMB->GetBlock(b-1)),
      this->GetInputArrayToProcess(0, inputAsMB->GetBlock(b+0))
    ))
      return 0;

    output->SetBlock(b-1, correlations);
  }

  if(inputAsDS){
    this->PreviousInput = vtkSmartPointer<vtkDataSet>::Take( inputAsDS->NewInstance() );
    this->PreviousInput->ShallowCopy(input);
  }

  return 1;
}
