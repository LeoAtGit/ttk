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

class vtkRenderer;
class vtkRenderWindow;
class vtkPointSet;
class vtkMultiBlockDataSet;
class vtkRenderPassCollection;
class vtkCamera;
class vtkAbstractArray;

namespace ttk {
  class ttkCinemaImagingVTK : virtual public Debug {
  public:
    ttkCinemaImagingVTK();
    ~ttkCinemaImagingVTK();

    int RenderImages(
        vtkMultiBlockDataSet* outputImages,

        vtkPointSet* inputObject,
        vtkPointSet* inputGrid,
        const double defaultResolution[2],
        const int    projectionMode,
        const double defaultCamAngle,
        const double defaultCamFocus[3],
        const double defaultCamNearFar[2],
        const double defaultCamHeight,
        const double defaultCamDir[3],
        const double defaultCamUp[3],
        const bool   generateScalarImages
    ) const;

  protected:
    int setupRenderer(
        vtkRenderer *renderer,
        vtkPointSet *object,
        vtkCamera *camera
    ) const;

    int setupWindow(
        vtkRenderWindow *window,
        vtkRenderer *renderer,
        const double resolution[2]
    ) const;

    int addValuePass(
        vtkPointSet *object,
        int fieldType,
        vtkRenderPassCollection *valuePassCollection,
        std::vector<std::string> &valuePassNames
    ) const;

    int addGobalArray(
        std::vector<vtkAbstractArray*> &globalArrays,
        std::string name,
        size_t nValues,
        const double *values
    ) const;
  };
};