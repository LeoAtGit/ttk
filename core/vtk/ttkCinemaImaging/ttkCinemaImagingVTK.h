#pragma once

#include <Debug.h>

class vtkRenderer;
class vtkRenderWindow;
class vtkPointSet;
class vtkMultiBlockDataSet;
class vtkRenderPassCollection;
class vtkCamera;

namespace ttk {
  class ttkCinemaImagingVTK : virtual public Debug {
  public:
    ttkCinemaImagingVTK();
    ~ttkCinemaImagingVTK();

    int RenderVTKObject(vtkMultiBlockDataSet *images,
                        vtkPointSet *object,
                        vtkPointSet *cameras) const;

  protected:
    int setupRenderer(vtkRenderer *renderer,
                      vtkPointSet *object,
                      vtkCamera *camera) const;

    int setupWindow(vtkRenderWindow *window, vtkRenderer *renderer) const;
    int addValuePass(vtkPointSet *object,
                     int fieldType,
                     vtkRenderPassCollection *valuePassCollection,
                     std::vector<std::string> &valuePassNames) const;
  };
}; // namespace ttk
