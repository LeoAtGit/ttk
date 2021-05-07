#pragma once

// VTK Module
#include <ttkCorrespondenceByMTSModule.h>

// VTK Includes
#include <ttkCorrespondenceAlgorithm.h>

// TTK Base Includes
#include <CorrespondenceByMTS.h>

class TTKCORRESPONDENCEBYMTS_EXPORT ttkCorrespondenceByMTS
  : public ttkCorrespondenceAlgorithm,
    protected ttk::CorrespondenceByMTS {

public:
  static ttkCorrespondenceByMTS *New();
  vtkTypeMacro(ttkCorrespondenceByMTS, ttkCorrespondenceAlgorithm);

protected:
  ttkCorrespondenceByMTS();
  ~ttkCorrespondenceByMTS() override;

  int Correlate(vtkImageData *correspondences,
                vtkDataObject *inputDataObjects0,
                vtkDataObject *inputDataObjects1
  ) override;
};
