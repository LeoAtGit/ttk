/// \ingroup vtk
/// \class ttkCinemaDarkroomCamera
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.11.2020
///
/// \brief This source generates a Cinema Darkroom Camera.
///
/// \param Output vtkUnstructuredGrid.
///
/// This source generates a single vertex with point data to represent a camera
/// that can be used for Cinema Darkroom rendering. The source can also be
/// synchronized with the camera of an active view port.
///
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkAlgorithm.h>
#include <ttkCinemaDarkroomModule.h>

class vtkPointSet;
class vtkPolyData;

class TTKCINEMADARKROOM_EXPORT ttkCinemaDarkroomCamera : public ttkAlgorithm {

  int Resolution[2]{256, 256};

  double Position[3]{0, 0, 0};
  double Direction[3]{0, 0, -1};
  double Up[3]{0, 1, 0};

  double NearFar[2]{0, 1};
  bool Orthographic{true};
  double Angle{90}; // perspective only
  double Scale{1}; // orthographic only

public:
  struct Calibration {
    struct Parameter {
      vtkDataArray *Array;
      std::vector<double> Get(const int idx = 0);
    };

    int nCameras{-1};
    bool Orthographic{false};

    Parameter Resolution;
    Parameter Position;
    Parameter Direction;
    Parameter Up;
    Parameter NearFar;
    Parameter Angle;
    Parameter Scale;

    Calibration(vtkPointSet *cameras);
  };

  static ttkCinemaDarkroomCamera *New();
  vtkTypeMacro(ttkCinemaDarkroomCamera, ttkAlgorithm);

  vtkSetVector2Macro(Resolution, int);
  vtkGetVector2Macro(Resolution, int);

  vtkSetVector3Macro(Position, double);
  vtkGetVector3Macro(Position, double);
  vtkSetVector3Macro(Direction, double);
  vtkGetVector3Macro(Direction, double);
  vtkSetVector3Macro(Up, double);
  vtkGetVector3Macro(Up, double);

  vtkSetVector2Macro(NearFar, double);
  vtkGetVector2Macro(NearFar, double);
  vtkSetMacro(Orthographic, bool);
  vtkGetMacro(Orthographic, bool);
  vtkSetMacro(Angle, double);
  vtkGetMacro(Angle, double);
  vtkSetMacro(Scale, double);
  vtkGetMacro(Scale, double);

  int SyncParaViewCamera();
  int SyncDarkroomCamera();

protected:
  // int computeFrustumGeometry();

  ttkCinemaDarkroomCamera();
  ~ttkCinemaDarkroomCamera() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int ComputeCalibrations(vtkPolyData *cameras,
                          vtkPointSet *positions = nullptr,
                          vtkDataObject *nearFarBoundingBoxObject = nullptr,
                          vtkDataObject *focusObject = nullptr);
  int ComputeFrustums(vtkPointSet *cameras, vtkPolyData *frustums);

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
