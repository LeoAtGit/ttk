/// TODO 4: Provide your information and **update** the documentation (in
/// particular regarding the order convention if input arrays need to be
/// specified with the standard VTK call SetInputArrayToProcess()).
///
/// \ingroup vtk
/// \class ttkBranchDecomposition
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// \brief TTK VTK-filter that wraps the ttk::BranchDecomposition module.
///
/// This VTK filter uses the ttk::BranchDecomposition module to compute an averaging of
/// the data values of an input point data array defined on the input
/// vtkDataSet.
///
/// \param Input vtkDataSet.
/// \param Output vtkDataSet.
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutputDataObject()).
///
/// The input data array needs to be specified via the standard VTK call
/// vtkAlgorithm::SetInputArrayToProcess() with the following parameters:
/// \param idx 0 (FIXED: the first array the algorithm requires)
/// \param port 0 (FIXED: first port)
/// \param connection 0 (FIXED: first connection)
/// \param fieldAssociation 0 (FIXED: point data)
/// \param arrayName (DYNAMIC: string identifier of the input array)
///
/// See the corresponding standalone program for a usage example:
///   - standalone/BranchDecomposition/main.cpp
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::BranchDecomposition
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkBranchDecompositionModule.h>

// VTK Includes
#include <ttkAlgorithm.h>

// TTK Base Includes
#include <BranchDecomposition.h>

class TTKBRANCHDECOMPOSITION_EXPORT ttkBranchDecomposition
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
  ,
    protected ttk::BranchDecomposition // and we inherit from the base class
{
private:

public:

  static ttkBranchDecomposition *New();
  vtkTypeMacro(ttkBranchDecomposition, ttkAlgorithm);

protected:
  ttkBranchDecomposition();
  ~ttkBranchDecomposition() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
