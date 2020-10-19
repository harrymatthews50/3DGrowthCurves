%% this script demonstrates how to make an age and sex-appropriate model fo normal variation


clear all; close all; restoredefaultpath;

% add path to supporting MATLAB Code
addpath(genpath('/usr/local/avalok/tmp/hmatth5/Chapter5/Submission1/3DGrowthCurves/ManuscriptSupportingCode/SupportingMATLABCode'));

%add path to patient assessment toolbox
addpath(genpath('/usr/local/avalok/tmp/hmatth5/Chapter5/Submission1/3DGrowthCurves/PatientAssessmenToolbox/'));

% add path to meshmonk toolbox - some classes for visualisation are needed
addpath(genpath('/usr/local/avalok/tmp/hmatth5/Projects/meshmonk'));

%% load data
% load training data
load('/usr/local/avalok/tmp/hmatth5/Chapter5/Submission1/3DGrowthCurves/ManuscriptSupportingCode/DATA/Simulated_Normative_Data/SimulatedNormativeData.mat');

%load MeshMonk Template
load('/usr/local/avalok/tmp/hmatth5/Projects/meshmonk/demo/Template.mat');


% load interpolants odf bandwidth from age
Interps = load('/usr/local/avalok/tmp/hmatth5/Chapter5/Submission1/BandwidthInterpolants/BandwidthInterpolants.mat');



%%
% set some properties
GC = growthCurve3D;
GC.TrainingLandmarks = shapes;
GC.Age = metadata.age;
GC.Sex = metadata.sexNumeric;
GC.ScaleShapes = true;
GC.RefScan = clone(Template);

GC.Speedy = false;

% set bandwidth adaptively given the evaluation age
GC.SetBandwidthFromInterpolant = true;
GC.SetBandwidthFromInterpolantInterpolant1 = Interps.maleInterp;
GC.SetBandwidthFromInterpolantInterpolant2 = Interps.femaleInterp;



%% Initialise
GC.initialize();

%% Create model for a given age

% 12 year old male
age = 12;
sex = 1;

GC.EvaluationSex = sex;
GC.EvaluationAge = age;

disp(['Model created with bandwidth = ',num2str(GC.BWSigma)]);
Model =GC.VariationModel;
% show the expected face
obj = clone(Model.Average);
obj.ViewMode = 'Solid';
obj.Material = 'Dull';
v = viewer(obj);
v.Tag = 'Expected Face';
v.SceneLightVisible = true;
v.SceneLightLinked = true;

%% plot color-map of standard deviations in the direction of the surface normals
obj = clone(Model.Average);
obj.ViewMode = 'Solid';
obj.Material = 'Dull';
v = viewer(obj);
v.Tag = 'Standard deviations in surface normal direction';
v.SceneLightVisible = true;
v.SceneLightLinked = true;
v.BackgroundColor = [1,1,1];
obj.VertexValue = Model.PointStandardDevs(:,4);
obj.ColorMode = 'Indexed';
colormap('hot');
colorbar();



%% plot the first PC

Model.animatePC(1,[-3,3])



