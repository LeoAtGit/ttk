#pragma once

// VTK Module
#include <ttkCorrespondenceByOverlapModule.h>

// VTK Includes
#include <ttkCorrespondenceAlgorithm.h>

// TTK Base Includes
#include <CorrespondenceByOverlap.h>

class TTKCORRESPONDENCEBYOVERLAP_EXPORT ttkCorrespondenceByOverlap
  : public ttkCorrespondenceAlgorithm,
    protected ttk::CorrespondenceByOverlap {
public:
  static ttkCorrespondenceByOverlap *New();
  vtkTypeMacro(ttkCorrespondenceByOverlap, ttkCorrespondenceAlgorithm);

protected:
  ttkCorrespondenceByOverlap();
  ~ttkCorrespondenceByOverlap();

  int ComputeCorrespondences(vtkImageData *correspondenceMatrix,
                             vtkDataObject *inputDataObjects0,
                             vtkDataObject *inputDataObjects1) override;
};
