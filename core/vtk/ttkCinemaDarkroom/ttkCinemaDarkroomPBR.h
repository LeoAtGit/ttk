/// \ingroup vtk
/// \class ttkCinemaDarkroomPBR
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
#include <ttkCinemaDarkroomShader.h>

class TTKCINEMADARKROOM_EXPORT ttkCinemaDarkroomPBR
  : public ttkCinemaDarkroomShader {
private:
  double Ambient{0.2};
  double AO{0.5};
  double Diffuse{1.0};
  double Roughness{0.2};
  double Metallic{0.0};

public:
  vtkSetMacro(Ambient, double);
  vtkGetMacro(Ambient, double);
  vtkSetMacro(AO, double);
  vtkGetMacro(AO, double);
  vtkSetMacro(Diffuse, double);
  vtkGetMacro(Diffuse, double);
  vtkSetMacro(Roughness, double);
  vtkGetMacro(Roughness, double);
  vtkSetMacro(Metallic, double);
  vtkGetMacro(Metallic, double);

  static ttkCinemaDarkroomPBR *New();
  vtkTypeMacro(ttkCinemaDarkroomPBR, ttkCinemaDarkroomShader);

protected:
  ttkCinemaDarkroomPBR();
  ~ttkCinemaDarkroomPBR() override;

  std::string GetFragmentShaderCode() override;

  int RegisterReplacements() override;
  int RegisterTextures(vtkImageData* image) override;
};
