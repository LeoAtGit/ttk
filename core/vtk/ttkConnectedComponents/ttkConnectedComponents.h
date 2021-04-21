/// \ingroup vtk
/// \class ttkConnectedComponents
/// \author Jonas Lukasczyk (jl@jluk.de)
/// \date 01.02.2019
///
/// \brief TTK VTK-filter that computes connected components based on a scalar
/// field.
///
/// VTK wrapping code for the @ConnectedComponents package.
///
/// This filter consumes a scalar field containing feature labels and computes
/// for each edge connected group of vertices with the same label a so-called
/// component, where negative labels represent the background. The computed components store the size and center of mass of each component. The module also assigns a unqiue positive integer to each component and maps
/// this id to the segmentation.
///
/// The input data array that contains the feature labels needs to be specified
/// via the standard VTK call SetInputArrayToProcess() with the following
/// parameters:
/// \param idx 0 (FIXED: the first array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the VTK array)
///
/// \sa ttk::ConnectedComponents

#pragma once

// VTK Module
#include <ttkConnectedComponentsModule.h>

// TTK Include
#include <ConnectedComponents.h>
#include <ttkAlgorithm.h>

class TTKCONNECTEDCOMPONENTS_EXPORT ttkConnectedComponents
  : public ttkAlgorithm,
    protected ttk::ConnectedComponents {

public:
  static ttkConnectedComponents *New();
  vtkTypeMacro(ttkConnectedComponents, ttkAlgorithm);

protected:
  ttkConnectedComponents();
  ~ttkConnectedComponents();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
