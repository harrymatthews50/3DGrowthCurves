# 3DGrowthCurves
A MATLAB repository to facilitate clinical assessment of normal and dysmorphic facial shape
Written by Hraold Matthews (harry.matthews@kuleuven.be) and Peter Claes (peter.claes@kuleuven.be).

![alt text](https://github.com/harrymatthews50/3DGrowthCurves/blob/master/img/3DGrowthCurves.png)

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

### Facial Assessment With the Facial Assessment Toolbox
#### Step 1 Landmarking Scans
Five landmarks need to be roughly positioned on each scan. 

![alt text](https://github.com/harrymatthews50/3DGrowthCurves/blob/master/img/Landmarks.png)

1. Open the Mevislab network FacialAssessmentToolbox/LandmarkIndicationMeVisLab/IndicateLandmarks_Meshes.mlab.
2.  Press the 'play' icon on the bubble in the bottom left 'SelectObjFiles'.
3. In the first dialog box select all the image files that you want to landmark.
4. Immediately afterwards a second dialog box wil open. Select the folder where ypu want to save the landmarks. Then the first scan should load. If you cannot see it, try moving the camera with scrollwheels in the viewer window; try pressing the go to home icon (button second frm the top in the viewer toolbar).
5. Indicate five landmarks in order (see above figure; subject's right endocanthion,left endocanthion, pronasale, right chelion, left chelion). These do not have to be precise, this should not take more than a few seconds. 
6. Press the 'play' icon on 'SaveMarkersAndLoadNext'. This will save the landmarks as a '.xml' file is the folder specified in step 3, and will load the next scan.

#### Step 2 Assessing the patient
The script FacialAssessmentToolbox/demo/DemoAssessment.m demonstrates how to perform an assessment in MATLAB.

1.  In the first cell set the paths appropriatels to the locations of the MeshMonk toolbox; the FacialAssessmentToolbox; and an appropriate place where you want to output the results.

2. In the third cell headed 'PatientInfo' set the path to to the '.obj' image file of the patient and the path to the '.xml' file of the landmarks indicated on the subject.

```matlab
% path to patient file
path2patient = '/usr/local/avalok/tmp/hmatth5/Chapter5/Submission1/3DGrowthCurvesPatientAssessmenToolbox/demo/demofaces/demoFace.obj';
% path to patient landmarks
path2patientLandmarks = '/usr/local/avalok/tmp/hmatth5/Projects/3DGrowthCurves/PatientAssessmenToolbox/demo/demofaces/demoFace.xml';
```

```javascript
function fancyAlert(arg) {
  if(arg) {
    $.facebox({div:'#foo'})
  }
}
```




### Manuscript Supporting Code
ManuscriptSupportingCode/AnalysisScripts contains three scripts, illustrating, on simulated data, how to build 3D growth curves from 3D surface scans that are in dense correspondence. Please see S1 Text of the manuscript for further detials on each step.

Script 1 and 2 address the tuning of the kernel bandidth
Scr1_ReconstructionErrorParameterSweep.m - computes reconstruction error for each training observation using different bandwidths
Scr2_fFtBandwidthTuningInterpolants.m - models how the bandwidth should change as a function of age to maninatin acceptable reconstruction error for each model
Scr3_MakeGrowthCurves.m - illustrates how to to produce an age and sex appropriate model for a given age and sex.


