/// \ingroup vtk
/// \class ttkCinemaDarkroomIBS
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.11.2020
///
/// \brief Image Based Shading
///
/// \param Input vtkImageData.
/// \param Output vtkImageData.
///
/// \b Related \b Publication:
/// "Dynamic Nested Tracking Graphs".
/// Jonas Lukasczyk, Christoph Garth, Gunther H. Weber, Tim Biedert, Ross
/// Maciejewski, Heike Leitte. IEEE Transactions on Visualization and Computer
/// Graphics. 2019
///
/// \sa ttkCinemaDarkroomShader

#pragma once

// VTK Module
#include <ttkCinemaDarkroomModule.h>
#include <ttkAlgorithm.h>

class TTKCINEMADARKROOM_EXPORT ttkCinemaDarkroomRendering
  : public ttkAlgorithm {
private:
  // CM
  std::string ManualColorMap{""};

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
  vtkSetMacro(ManualColorMap, const std::string &);
  vtkGetMacro(ManualColorMap, std::string);

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

  static ttkCinemaDarkroomRendering *New();
  vtkTypeMacro(ttkCinemaDarkroomRendering, ttkAlgorithm);

protected:
  ttkCinemaDarkroomRendering();
  ~ttkCinemaDarkroomRendering() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
