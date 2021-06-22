#pragma once

// VTK Module
#include <ttkCorrespondenceByGradientModule.h>

// VTK Includes
#include <ttkCorrespondenceAlgorithm.h>

// TTK Base Includes
#include <CorrespondenceByGradient.h>

class TTKCORRESPONDENCEBYGRADIENT_EXPORT ttkCorrespondenceByGradient
  : public ttkCorrespondenceAlgorithm,
    protected ttk::CorrespondenceByGradient {

public:
  static ttkCorrespondenceByGradient *New();
  vtkTypeMacro(ttkCorrespondenceByGradient, ttkCorrespondenceAlgorithm);

protected:
  ttkCorrespondenceByGradient();
  ~ttkCorrespondenceByGradient();

  int ComputeCorrespondences(vtkImageData *correspondenceMatrix,
                             vtkDataObject *inputDataObjects0,
                             vtkDataObject *inputDataObjects1) override;
};