/// \ingroup vtk
/// \class ttkCinemaImaging
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 1.9.2018
///
/// \brief TTK VTK-filter that generates images of a vtkDataSet.
///
/// This filter takes images of a vtkDataObject from positions specified on a
/// vtkPointSet. Each image will be a block of a vtkMultiBlockDataSet where
/// block order corresponds to point order. Each sample point can optionally
/// have vtkDoubleArrays to override the default rendering parameters, i.e, the
/// resolution, focus, clipping planes, and viewport height.
///
/// VTK wrapping code for the @CinemaImaging package.
///
/// \param Input vtkDataObject that will be depicted (vtkDataObject)
/// \param Input vtkPointSet that records the camera sampling locations
/// (vtkPointSet) \param Output vtkMultiBlockDataSet that represents a list of
/// images (vtkMultiBlockDataSet)

#pragma once

// VTK Module
#include <ttkCinemaImagingModule.h>

// TTK includes
#include <ttkAlgorithm.h>

class vtkPolyData;
class vtkMultiBlockDataSet;
class vtkPointSet;
class vtkFieldData;
class vtkImageData;
class vtkPointData;

class TTKCINEMAIMAGING_EXPORT ttkCinemaImaging : public ttkAlgorithm {

private:
  int Backend{0};

  int Resolution[2]{256, 256};
  bool GenerateScalarImages{true};

  int CamProjectionMode{0};
  bool CamFocusAuto{true};
  double CamFocus[3]{0, 0, 0};
  bool CamNearFarAuto{true};
  double CamNearFar[2]{0, 1};

  bool CamHeightAuto{true};
  double CamHeight{1};
  double CamAngle{90};

public:
  static ttkCinemaImaging *New();
  vtkTypeMacro(ttkCinemaImaging, ttkAlgorithm);

  // Backend
  vtkSetMacro(Backend, int);
  vtkGetMacro(Backend, int);

  // General Settings
  vtkSetVector2Macro(Resolution, int);
  vtkGetVector2Macro(Resolution, int);
  vtkSetMacro(GenerateScalarImages, bool);
  vtkGetMacro(GenerateScalarImages, bool);

  // Camera
  vtkSetMacro(CamProjectionMode, int);
  vtkGetMacro(CamProjectionMode, int);

  vtkSetMacro(CamFocusAuto, bool);
  vtkGetMacro(CamFocusAuto, bool);
  vtkSetVector3Macro(CamFocus, double);
  vtkGetVector3Macro(CamFocus, double);

  vtkSetMacro(CamNearFarAuto, bool);
  vtkGetMacro(CamNearFarAuto, bool);
  vtkSetVector2Macro(CamNearFar, double);
  vtkGetVector2Macro(CamNearFar, double);

  // Orthographic
  vtkSetMacro(CamHeightAuto, bool);
  vtkGetMacro(CamHeightAuto, bool);
  vtkSetMacro(CamHeight, double);
  vtkGetMacro(CamHeight, double);

  // Perspective
  vtkSetMacro(CamAngle, double);
  vtkGetMacro(CamAngle, double);

  static int addFieldDataArray(
    vtkFieldData* fd,
    vtkDataArray* array,
    int tupelIdx,
    std::string name = ""
  );

  static int addAllFieldDataArrays(
    vtkPointSet* inputGrid,
    vtkImageData* image,
    int tupelIdx
  );

  static int computeDirFromFocus(
    vtkPointSet* inputGrid
  );

  static int ensureGridData(
    vtkPointData* fd,
    std::string name,
    int nTuples,
    const std::vector<double>& defaultValues
  );

protected:
  ttkCinemaImaging();
  ~ttkCinemaImaging();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
