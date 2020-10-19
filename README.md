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
3) This from the MATLAB file exchange  https://www.mathworks.com/matlabcentral/fileexchange/10056-scattered-data-interpolation-and-approximation-using-radial-base-functions.
4) Running a patient assessment requires some crude positioning landmarks be placed on each face. We provide a network for doing this. To run this you will need to have MeVisLab installed (free version is fine). https://www.mevislab.de/download; you may also need to edit your meshes prior to processing. This can be done in meshlab https://www.meshlab.net. this tutorial miht be useful https://www.3printr.com/meshlab-bearbeitung-der-polygon-oberflaechennetze-von-3d-modellen-199419/

### Facial Assessment Toolbox
#### Step 1 Landmarking Scans
Five landmarks need to be roughly positioned on each scan. 

[[https://github.com/harrymatthews50/3DGrowthCurves/blob/master/img/Landmarks.png]]

1. Open the Mevislab network FacialAssessmentToolbox.LandmarkIndicationMeVisLab/IndicateLandmarks_Meshes.mlab.
2.  Press the 'play' icon on the bubble in the bottom left 'SelectObjFiles'.
3. In the first dialog box select all the image files that you want to landmark.
4. In the second dialog box select the folder where ypu want to save 





### Manuscript Supporting Code
ManuscriptSupportingCode/AnalysisScripts contains three scripts, illustrating, on simulated data, how to build 3D growth curves from 3D surface scans that are in dense correspondence. Please see S1 Text of the manuscript for further detials on each step.

Script 1 and 2 address the tuning of the kernel bandidth
Scr1_ReconstructionErrorParameterSweep.m - computes reconstruction error for each training observation using different bandwidths
Scr2_fFtBandwidthTuningInterpolants.m - models how the bandwidth should change as a function of age to maninatin acceptable reconstruction error for each model
Scr3_MakeGrowthCurves.m - illustrates how to to produce an age and sex appropriate model for a given age and sex.


