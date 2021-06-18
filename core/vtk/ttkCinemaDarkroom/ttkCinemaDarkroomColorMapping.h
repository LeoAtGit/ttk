/// \ingroup vtk
/// \class ttkCinemaDarkroomColorMapping
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.11.2020
///
/// \brief Performs color mapping of a scalar field.
///
/// \param Input vtkImageData.
/// \param Output vtkImageData.
///
/// This class maps each value of a scalar point data array to either a solid
/// color, a predefined color map, or a manually defined color map.
///
/// \b Related \b Publication:
/// "Cinema Database Specification - Dietrich Release v1.2".
/// D. Rogers, J. Woodring, J. Ahrens, J. Patchett, and J. Lukasczyk.
/// Technical Report LA-UR-17-25072, Los Alamos National Laboratory,
/// 2018.
///
/// \sa ttkCinemaDarkroomShader

#pragma once

// VTK Module
#include <ttkCinemaDarkroomModule.h>
#include <ttkCinemaDarkroomShader.h>

class TTKCINEMADARKROOM_EXPORT ttkCinemaDarkroomColorMapping
  : public ttkCinemaDarkroomShader {
private:
  static const std::vector<std::vector<double>> ColorMaps;

  double ScalarRange[2]{0, 1};
  int ColorMap{0};
  std::string ColorMapData{""};
  double SingleColor[3]{0, 0, 0};
  double NANColor[3]{0, 0, 0};
  bool TransparentNAN{true};

public:
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

  int SyncColorMapsWithParaView();

  static ttkCinemaDarkroomColorMapping *New();
  vtkTypeMacro(ttkCinemaDarkroomColorMapping, ttkCinemaDarkroomShader);

protected:
  ttkCinemaDarkroomColorMapping();
  ~ttkCinemaDarkroomColorMapping() override;

  int Render(vtkImageData *image) override;
};