#pragma once

#include <CinemaImagingNative.h>

class vtkMultiBlockDataSet;
class vtkPointSet;

namespace ttk {
  class ttkCinemaImagingNative : public CinemaImagingNative {
  public:
    ttkCinemaImagingNative();
    ~ttkCinemaImagingNative();

    int RenderVTKObject(vtkMultiBlockDataSet *images,

                        vtkPointSet *object,
                        vtkPointSet *cameras) const;
  };
}; // namespace ttk
