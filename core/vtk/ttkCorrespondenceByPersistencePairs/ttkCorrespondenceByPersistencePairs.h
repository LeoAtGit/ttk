#pragma once

// VTK Module
#include <ttkCorrespondenceByPersistencePairsModule.h>

// VTK Includes
#include <ttkCorrespondenceAlgorithm.h>

// TTK Base Includes
#include <CorrespondenceByPersistencePairs.h>

class TTKCORRESPONDENCEBYPERSISTENCEPAIRS_EXPORT ttkCorrespondenceByPersistencePairs
  : public ttkCorrespondenceAlgorithm,
    protected ttk::CorrespondenceByPersistencePairs {

public:
  static ttkCorrespondenceByPersistencePairs *New();
  vtkTypeMacro(ttkCorrespondenceByPersistencePairs, ttkCorrespondenceAlgorithm);

protected:
  ttkCorrespondenceByPersistencePairs();
  ~ttkCorrespondenceByPersistencePairs();

  int Correlate(vtkImageData *correspondences,
                vtkDataObject *inputDataObjects0,
                vtkDataObject *inputDataObjects1
  ) override;
};
