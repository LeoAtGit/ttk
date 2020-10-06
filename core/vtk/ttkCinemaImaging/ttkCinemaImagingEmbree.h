/// \ingroup base
/// \class ttk::EmbreeRenderer
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.05.2020
///
/// \brief TTK %EmbreeRenderer processing package.
///
/// %EmbreeRenderer is a TTK processing package that

#pragma once

#include <Debug.h>

class vtkMultiBlockDataSet;
class vtkCellArray;
class vtkPointSet;

namespace ttk {
  class ttkCinemaImagingEmbree : virtual public Debug {
  public:
    ttkCinemaImagingEmbree();
    ~ttkCinemaImagingEmbree();

    int RenderImages(
        vtkMultiBlockDataSet* outputImages,

        vtkCellArray* inputObjectCells,
        vtkPointSet* inputObject,
        vtkPointSet* inputGrid
    ) const;
  };
};