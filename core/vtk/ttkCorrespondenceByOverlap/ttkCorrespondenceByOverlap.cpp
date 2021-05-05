#include <ttkCorrespondenceByOverlap.h>

#include <vtkObjectFactory.h>

#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkIntArray.h>
#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointData.h>

#include <ttkUtils.h>

vtkStandardNewMacro(ttkCorrespondenceByOverlap);

ttkCorrespondenceByOverlap::ttkCorrespondenceByOverlap() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

ttkCorrespondenceByOverlap::~ttkCorrespondenceByOverlap() {
}

int ttkCorrespondenceByOverlap::Correlate(vtkImageData *correspondences,
                                      vtkDataObject *inputDataObjects0,
                                      vtkDataObject *inputDataObjects1
                                      ) {

  if(this->GetInputArrayAssociation(0, inputDataObjects0) != 0)
    return !this->printErr("Labels must be point data.");

  // get label arrays
  auto labels0 = this->GetInputArrayToProcess(0, inputDataObjects0);
  auto labels1 = this->GetInputArrayToProcess(0, inputDataObjects1);

  // validate arrays
  if(!labels0 || !labels1)
    return !this->printErr("Unable to retrieve labels.");

  if(labels0->GetNumberOfComponents() != 1
     || labels1->GetNumberOfComponents() != 1)
    return !this->printErr("Labels must have exactly one component.");

  if(labels0->GetNumberOfTuples() != labels1->GetNumberOfTuples())
    return !this->printErr("Labels must have same number of values.");

  if(labels0->GetDataType() != labels1->GetDataType())
    return !this->printErr("Labels must have same data type.");

  const int nVertices = labels0->GetNumberOfTuples();

  // extract unique labels from volume
  ttk::CorrespondenceByOverlap::LabelIndexMap labelsIndexMap0;
  ttk::CorrespondenceByOverlap::LabelIndexMap labelsIndexMap1;
  int status = 0;
  for(auto& it : std::vector<std::pair<vtkDataArray*,LabelIndexMap*>>({{labels0,&labelsIndexMap0},{labels1,&labelsIndexMap1}})){
    switch(labels0->GetDataType()) {
      vtkTemplateMacro(status = this->computeLabelIndexMap<VTK_TT>(
                         *it.second,
                         ttkUtils::GetPointer<const VTK_TT>(it.first),
                         nVertices
                       )
      );
      if(!status)
        return 0;
    }
  }

  const int nLabels0 = labelsIndexMap0.size();
  const int nLabels1 = labelsIndexMap1.size();

  // initialize correspondence matrix
  correspondences->SetDimensions(nLabels0, nLabels1, 1);
  correspondences->AllocateScalars(VTK_INT, 1);
  auto correspondencesArray = correspondences->GetPointData()->GetArray(0);
  correspondencesArray->SetName("Overlap");

  // compute overlaps
  switch(labels0->GetDataType()) {
    vtkTemplateMacro(status = this->computeAdjacencyMatrix<VTK_TT>(
                      ttkUtils::GetPointer<int>(correspondencesArray),
                      ttkUtils::GetPointer<const VTK_TT>(labels0),
                      ttkUtils::GetPointer<const VTK_TT>(labels1),
                      nVertices,
                      labelsIndexMap0,
                      labelsIndexMap1
    ));
  }
  if(!status)
    return 0;

  // add index label maps
  int a=0;
  for(const auto it : std::vector<LabelIndexMap*>({&labelsIndexMap0,&labelsIndexMap1})){
    auto array = vtkSmartPointer<vtkIntArray>::New();
    array->SetName(std::string("IndexLabelMap"+std::to_string(a++)).data());
    array->SetNumberOfTuples(it->size());
    auto arrayData = ttkUtils::GetPointer<int>(array);
    for(const auto& it2: *it)
      arrayData[it2.second] = it2.first;

    correspondences->GetFieldData()->AddArray(array);
  }

  return 1;
}