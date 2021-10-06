/// \ingroup vtk
/// \class ttkCinemaDarkroomMaster
/// \author Rosty Hnatyshyn <rostyslav.hnatyshyn@gmail.com> and Jonas Lukasczyk
/// <jl@jluk.de> \date 01.06.2021
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
#include <ttkAlgorithm.h>
#include <ttkCinemaDarkroomModule.h>

#include <ttkMacros.h>

class TTKCINEMADARKROOM_EXPORT ttkCinemaDarkroomMaster : public ttkAlgorithm {

public:
  enum class SHADER { IBS = 0, PBR = 1, PHONG = 2 };

private:
  SHADER Shader{SHADER::PBR};

  bool UseMSAA{false};

  // CM
  double ScalarRange[2]{0, 1};
  int ColorMap{0};
  std::string ColorMapData{""};
  double SingleColor[3]{0, 0, 0};
  double NANColor[3]{0, 0, 0};
  bool TransparentNAN{true};

  // SSSAO
  int Samples{32};
  double Radius{0.05};
  double DiffArea{0.5};

  // General Shader Parameters
  double Ambient{0.2};
  double Diffuse{1.0};
  double AO{1.0};

  // IBS
  double Strength{1.0};

  // PBR
  double Metallic{0.1};
  double Roughness{0.1};

  // Phong
  double Exponent{16.0};
  double Specular{1.0};

  // //DoF
  // bool UseSSDoF{false};
  // double DepthRadius{0.05};
  // double Aperture{1};
  // double FocalDepth{1};
  // double MaxBlur{1};

  // FXAA
  bool UseFXAA{true};

public:
  ttkSetEnumMacro(Shader, SHADER);
  vtkGetEnumMacro(Shader, SHADER);

  vtkSetMacro(UseMSAA, bool);
  vtkGetMacro(UseMSAA, bool);

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
  vtkSetMacro(TransparentNAN, bool);
  vtkGetMacro(TransparentNAN, bool);

  // SSSAO
  vtkSetMacro(Samples, int);
  vtkGetMacro(Samples, int);
  vtkSetMacro(Radius, double);
  vtkGetMacro(Radius, double);
  vtkSetMacro(DiffArea, double);
  vtkGetMacro(DiffArea, double);

  // General Shader Parameters
  vtkSetMacro(Ambient, double);
  vtkGetMacro(Ambient, double);
  vtkSetMacro(Diffuse, double);
  vtkGetMacro(Diffuse, double);
  vtkSetMacro(AO, double);
  vtkGetMacro(AO, double);

  // IBS
  vtkSetMacro(Strength, double);
  vtkGetMacro(Strength, double);

  // PBR
  vtkSetMacro(Metallic, double);
  vtkGetMacro(Metallic, double);
  vtkSetMacro(Roughness, double);
  vtkGetMacro(Roughness, double);

  // Phong
  vtkSetMacro(Exponent, double);
  vtkGetMacro(Exponent, double);
  vtkSetMacro(Specular, double);
  vtkGetMacro(Specular, double);

  // //DoF
  // vtkSetMacro(UseSSDoF, bool);
  // vtkGetMacro(UseSSDoF, bool);
  // vtkSetMacro(DepthRadius, double);
  // vtkGetMacro(DepthRadius, double);
  // vtkSetMacro(Aperture, double);
  // vtkGetMacro(Aperture, double);
  // vtkSetMacro(FocalDepth, double);
  // vtkGetMacro(FocalDepth, double);
  // vtkSetMacro(MaxBlur, double);
  // vtkGetMacro(MaxBlur, double);

  // FXAA
  vtkSetMacro(UseFXAA, bool);
  vtkGetMacro(UseFXAA, bool);

  static ttkCinemaDarkroomMaster *New();
  vtkTypeMacro(ttkCinemaDarkroomMaster, ttkAlgorithm);

protected:
  ttkCinemaDarkroomMaster();
  ~ttkCinemaDarkroomMaster() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
