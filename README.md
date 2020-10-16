# 3DGrowthCurves
A MATLAB repository to facilitate clinical assessment of normal and dysmorphic facial shape
Written by Hraold Matthews (harry.matthews@kuleuven.be) and Peter Claes (peter.claes@kuleuven.be).

[Insert illustration]
## Contents
This repository contains:
1) Code supporting replication of the analyses in [REF].
2) The Facial Assessment Toolbox, for assessing normal and abnormal facial shape.

## Usage
### Dependencies
1) A basic MATLAB installation.
2) The MeshMonk Toolbox https://github.com/TheWebMonks/meshmonk.
3) This from the MATLAB file exchange  https://www.mathworks.com/matlabcentral/fileexchange/10056-scattered-data-interpolation-and-approximation-using-radial-base-functions

### Facial Assessment Toolbox
Please see FacialAssessmentToolbox/demo/PatientAssessmentInstructions.docx for detailed instructions on how to assess new patients with the toolbos

### Manuscript Supporting Code
ManuscriptSupportingCode/AnalysisScripts contains three scripts, illustrating, on simulated data, how to build 3D growth curves from 3D surface scans that are in dense correspondence. Please see S1 Text of the manuscript for further detials on each step.

Script 1 and 2 address the tuning of the kernel bandidth
Scr1_ReconstructionErrorParameterSweep.m - computes reconstruction error for each training observation using different bandwidths
Scr2_fFtBandwidthTuningInterpolants.m - models how the bandwidth should change as a function of age to maninatin acceptable reconstruction error for each model
Scr3_MakeGrowthCurves.m - illustrates how to to produce an age and sex appropriate model for a given age and sex.


