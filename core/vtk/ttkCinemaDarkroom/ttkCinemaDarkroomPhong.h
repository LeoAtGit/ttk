/// \ingroup vtk
/// \class ttkCinemaDarkroomPhong
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

class TTKCINEMADARKROOM_EXPORT ttkCinemaDarkroomPhong
  : public ttkCinemaDarkroomShader {
private:
  double Exponent{16.0};
  double Ambient{0.7};
  double Diffuse{0.7};
  double Specular{1.0};

public:
  vtkSetMacro(Exponent, double);
  vtkGetMacro(Exponent, double);
  vtkSetMacro(Ambient, double);
  vtkGetMacro(Ambient, double);
  vtkSetMacro(Diffuse, double);
  vtkGetMacro(Diffuse, double);
  vtkSetMacro(Specular, double);
  vtkGetMacro(Specular, double);

  static ttkCinemaDarkroomPhong *New();
  vtkTypeMacro(ttkCinemaDarkroomPhong, ttkCinemaDarkroomShader);

protected:
  ttkCinemaDarkroomPhong();
  ~ttkCinemaDarkroomPhong() override;

  std::string GetFragmentShaderCode() override;

  int RegisterReplacements() override;
  int RegisterTextures(vtkImageData *image) override;
};
