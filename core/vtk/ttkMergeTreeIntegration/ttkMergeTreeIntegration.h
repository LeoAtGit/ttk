#pragma once

// VTK Module
#include <ttkMergeTreeIntegrationModule.h>

// Base Includes
#include <MergeTreeIntegration.h>

// VTK Includes
#include <ttkAlgorithm.h>

class TTKMERGETREEINTEGRATION_EXPORT ttkMergeTreeIntegration : public ttkAlgorithm, public ttk::MergeTreeIntegration
{
    public:
        static ttkMergeTreeIntegration *New();
        vtkTypeMacro(ttkMergeTreeIntegration, ttkAlgorithm);

    protected:
        ttkMergeTreeIntegration();
        ~ttkMergeTreeIntegration() override;

        int FillInputPortInformation(int port, vtkInformation* info) override;
        int FillOutputPortInformation(int port, vtkInformation* info) override;
        int RequestData(
            vtkInformation* request,
            vtkInformationVector** inputVector,
            vtkInformationVector* outputVector
        ) override;
};