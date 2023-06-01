# 3DGrowthCurves
A MATLAB repository to facilitate clinical assessment of normal and dysmorphic facial shape
Written by Harold Matthews (harry.matthews@kuleuven.be) and Peter Claes (peter.claes@kuleuven.be). If you find the toolbox useful please acknowledge and cite 
us:

Matthews, H.S., Palmer, R.L., Baynam, G.S. et al. Large-scale open-source three-dimensional growth curves for clinical facial assessment and objective description of facial dysmorphism. Sci Rep 11, 12175 (2021). https://doi.org/10.1038/s41598-021-91465-z

![alt text](https://github.com/harrymatthews50/3DGrowthCurves/blob/master/img/3DGrowthCurves.png)

The functionality of this toolbox is now available in the 3D MedX software produced by the 3D Lab of Radboudmc, Nijmegen (https://www.3dmedx.nl)


## Contents
This repository contains:
1) Code supporting replication of the analyses in:
Matthews, H.S., Palmer, R.L., Baynam, G.S. et al. Large-scale open-source three-dimensional growth curves for clinical facial assessment and objective description of facial dysmorphism. Sci Rep 11, 12175 (2021). https://doi.org/10.1038/s41598-021-91465-z

4) The Facial Assessment Toolbox, for assessing normal and abnormal facial shape.

## Dependencies
1) A basic MATLAB installation.
2) The MeshMonk Toolbox https://github.com/TheWebMonks/meshmonk.
3) This from the MATLAB file exchange  https://www.mathworks.com/matlabcentral/fileexchange/10056-scattered-data-interpolation-and-approximation-using-radial-base-functions.
4) Running a patient assessment requires some crude positioning landmarks be placed on each face. We provide a network for doing this. To run this you will need to have MeVisLab installed (free version is fine). https://www.mevislab.de/download; you may also need to edit your meshes prior to processing. This can be done using the python [MeshEditor](https://github.com/harrymatthews50/MeshEditor).

## Facial Assessment With the Facial Assessment Toolbox
### Part 1 Landmarking Scans
Five landmarks need to be roughly positioned on each scan. 

![alt text](https://github.com/harrymatthews50/3DGrowthCurves/blob/master/img/Landmarks.png)

1. Open the Mevislab network FacialAssessmentToolbox/LandmarkIndicationMeVisLab/IndicateLandmarks_Meshes.mlab.
2.  Press the 'play' icon on the bubble in the bottom left 'SelectObjFiles'.
3. In the first dialog box select all the image files that you want to landmark.
4. Immediately afterwards a second dialog box wil open. Select the folder where ypu want to save the landmarks. Then the first scan should load. If you cannot see it, try moving the camera with scrollwheels in the viewer window; try pressing the go to home icon (button second frm the top in the viewer toolbar).
5. Indicate five landmarks in order (see above figure; subject's right endocanthion,left endocanthion, pronasale, right chelion, left chelion). These do not have to be precise, this should not take more than a few seconds. 
6. Press the 'play' icon on 'SaveMarkersAndLoadNext'. This will save the landmarks as a '.xml' file is the folder specified in step 3, and will load the next scan.

### Part 2 Assessing the patient
The script FacialAssessmentToolbox/demo/DemoAssessment.m demonstrates how to perform an assessment in MATLAB.

1. Set up the script
*  In the first cell set the paths appropriatels to the locations of the MeshMonk toolbox; the FacialAssessmentToolbox; and an appropriate place where you want to output the results.

* In Cell 3 set the path to to the '.obj' image file of the patient and the path to the '.xml' file of the landmarks indicated on the subject. And specify the patient's sex and age.

```matlab
%% Cell 3 Edit Patient Info

% path to patient file
path2patient = '/usr/local/avalok/tmp/hmatth5/Chapter5/Submission1/3DGrowthCurvesPatientAssessmenToolbox/demo/demofaces/demoFace.obj';
% path to patient landmarks
path2patientLandmarks = '/usr/local/avalok/tmp/hmatth5/Projects/3DGrowthCurves/PatientAssessmenToolbox/demo/demofaces/demoFace.xml';

% specify patient's age and sex
sex = 'M';
age = 35;
```

2. Process the image with MeshMonk
This determines the locations of 7160 quasi-landmarks on the facial scan of the image. The accuracy of this process can be checked by assessing whether certain landmark points end up on the correct locations on the face.
* Run the script up to cell 7 inclusive. Two windows will open. On eshows the locations of some landmarks on the template face. The other shows some landmarks, automatically found by MeshMonk on the face of the patient. The accuracy of these ladnmarks on the patient indicate how well MeshMonk has established the quasi-landmarks, which determines the accuracy of the comparison of the patient to the Growth Curve.

![alt text](https://github.com/harrymatthews50/3DGrowthCurves/blob/master/img/inspectCorrespondence.png)

3. Assess patient
* Run Cell 8 to produce the facial signature of the patient. 

![alt text](https://github.com/harrymatthews50/3DGrowthCurves/blob/master/img/demoFacesignature_signatureWeight=62.2219.png)


* Run Cell 9 to produce the normal equivalent of the patient. As the demo patient has no major facial abnormalities the normal equivalent is essentially the same as the image of the patient.

![alt text](https://github.com/harrymatthews50/3DGrowthCurves/blob/master/img/NormalEquivalentResults.png)



## Manuscript Supporting Code
ManuscriptSupportingCode/AnalysisScripts contains three scripts, illustrating, on simulated data, how to build 3D growth curves from 3D surface scans that are in dense correspondence. Please see S1 Text of the manuscript for further detials on each step.

This simulated data can be downloaded here: https://drive.google.com/drive/folders/1KzElGVdS-aByku4nVL5Wq2X8Z205M7Ed?usp=sharing

Script 1 and 2 address the tuning of the kernel bandidth
Scr1_ReconstructionErrorParameterSweep.m - computes reconstruction error for each training observation using different bandwidths
Scr2_fFtBandwidthTuningInterpolants.m - models how the bandwidth should change as a function of age to maninatin acceptable reconstruction error for each model
Scr3_MakeGrowthCurves.m - illustrates how to to produce an age and sex appropriate model for a given age and sex.


