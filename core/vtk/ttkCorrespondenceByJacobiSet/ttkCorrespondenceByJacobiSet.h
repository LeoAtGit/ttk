#pragma once

// VTK Module
#include <ttkCorrespondenceByJacobiSetModule.h>

// VTK Includes
#include <ttkCorrespondenceAlgorithm.h>

class TTKCORRESPONDENCEBYJACOBISET_EXPORT ttkCorrespondenceByJacobiSet
  : public ttkCorrespondenceAlgorithm {

public:
  static ttkCorrespondenceByJacobiSet *New();
  vtkTypeMacro(ttkCorrespondenceByJacobiSet, ttkCorrespondenceAlgorithm);

protected:
  ttkCorrespondenceByJacobiSet();
  ~ttkCorrespondenceByJacobiSet();

  int ComputeCorrespondences(vtkImageData *correspondenceMatrix,
                             vtkDataObject *inputDataObjects0,
                             vtkDataObject *inputDataObjects1) override;
};
