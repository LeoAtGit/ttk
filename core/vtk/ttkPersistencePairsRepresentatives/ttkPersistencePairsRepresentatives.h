/// \ingroup vtk
/// \class ttkPersistencePairsRepresentatives
/// \author Maxime Soler
/// \date 09.06.2021
///
/// \brief TTK VTK-filter that computes point representatives from a persistence diagram.
///
/// This filter consumes a ttkPersistenceDiagram, outputs a vtkPointSet a domain-embedded point
/// (x,y,z) for each input persistence pair.
/// The order in which points are laid out must be the same as ttkCorrespondenceByPersistencePairs.
///

#pragma once

// VTK Module
#include <ttkPersistencePairsRepresentativesModule.h>
#include <vtkUnstructuredGrid.h>

// TTK Include
//#include <ConnectedComponents.h>
#include <ttkAlgorithm.h>

class TTKPERSISTENCEPAIRSREPRESENTATIVES_EXPORT ttkPersistencePairsRepresentatives
  : public ttkAlgorithm {

private:
  //bool UseSeedIdAsComponentId{true};

public:
  //vtkSetMacro(UseSeedIdAsComponentId, bool);
  //vtkGetMacro(UseSeedIdAsComponentId, bool);

  static ttkPersistencePairsRepresentatives *New();
  vtkTypeMacro(ttkPersistencePairsRepresentatives, ttkAlgorithm);

  template <typename dataType>
  int getDiagram(
      std::vector<std::tuple<int, ttk::CriticalType, int, ttk::CriticalType, dataType,
          int, dataType, float, float, float, dataType, float, float, float> > &diagram,
      vtkUnstructuredGrid *CTPersistenceDiagram_,
      const double spacing,
      const int diagramNumber);

protected:
  ttkPersistencePairsRepresentatives();
  ~ttkPersistencePairsRepresentatives();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
