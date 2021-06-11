/// \ingroup vtk
/// \class ttkCinemaDarkroomShading
/// \author Rosty Hnatyshyn <rostyslav.hnatyshyn@gmail.com> and Jonas Lukasczyk <jl@jluk.de>
/// \date 01.06.2021
///
/// \brief Image Based Shading
///
/// \param Input vtkImageData.
/// \param Output vtkImageData.
///
/// \b Related \b Publication:
/// "Cinema Darkroom: A Deferred Rendering Framework for Large-Scale Datasets".
/// J. Lukasczyk, C. Garth, M. Larsen, W. Engelke, I. Hotz, D. Rogers, J.
/// Ahrens, and R. Maciejewski. IEEE 10th Symposium on Large Data Analysis and
/// Visualization (LDAV), 2020.
///
/// \sa ttkCinemaDarkroomShader

#pragma once

// VTK Module
#include <ttkCinemaDarkroomModule.h>
#include <ttkAlgorithm.h>

class TTKCINEMADARKROOM_EXPORT ttkCinemaDarkroomShading
  : public ttkAlgorithm {
private:
  // CM
  double ScalarRange[2]{0, 1};
  int ColorMap{0};
  std::string ColorMapData{""};
  double SingleColor[3]{0, 0, 0};
  double NANColor[3]{0, 0, 0};

  //SSSAO
  int Samples{1};
  double Radius{0.05};
  double DiffArea{0.5};

  //IBS
  double Strength{1.0};
  double Luminance{1.0};
  double Ambient{0.2};

  //DoF
  bool UseSSDoF{false};
  double DepthRadius{0.05};
  double Aperture{1};
  double FocalDepth{1};
  double MaxBlur{1};

  //FXAA
  bool UseFXAA{true};

public:
  // CM
  vtkSetVector2Macro(ScalarRange, double);
  vtkGetVector2Macro(ScalarRange, double);
  vtkSetMacro(ColorMap, int);
  vtkGetMacro(ColorMap, int);
  vtkSetMacro(ColorMapData, const std::string &);
  vtkGetMacro(ColorMapData, std::string);
  vtkSetVector3Macro(NANColor, double);
  vtkGetVector3Macro(NANColor, double);
  vtkSetVector3Macro(SingleColor, double);
  vtkGetVector3Macro(SingleColor, double);

  //SSSAO
  vtkSetMacro(Samples, int);
  vtkGetMacro(Samples, int);
  vtkSetMacro(Radius, double);
  vtkGetMacro(Radius, double);
  vtkSetMacro(DiffArea, double);
  vtkGetMacro(DiffArea, double);

  //IBS
  vtkSetMacro(Strength, double);
  vtkGetMacro(Strength, double);
  vtkSetMacro(Luminance, double);
  vtkGetMacro(Luminance, double);
  vtkSetMacro(Ambient, double);
  vtkGetMacro(Ambient, double);

  //DoF
  vtkSetMacro(UseSSDoF, int);
  vtkGetMacro(UseSSDoF,int);
  vtkSetMacro(DepthRadius, double);
  vtkGetMacro(DepthRadius, double);
  vtkSetMacro(Aperture, double);
  vtkGetMacro(Aperture, double);
  vtkSetMacro(FocalDepth, double);
  vtkGetMacro(FocalDepth, double);
  vtkSetMacro(MaxBlur, double);
  vtkGetMacro(MaxBlur, double);

  //FXAA
  vtkSetMacro(UseFXAA, int);
  vtkGetMacro(UseFXAA, int);

  static ttkCinemaDarkroomShading *New();
  vtkTypeMacro(ttkCinemaDarkroomShading, ttkAlgorithm);

protected:
  ttkCinemaDarkroomShading();
  ~ttkCinemaDarkroomShading() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
