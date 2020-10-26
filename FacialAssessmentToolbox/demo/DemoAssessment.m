
%% Cell 1 SET PATHS
clear all; close all; restoredefaultpath


% add path to 'meshmonk' toolbox on your machine
path2meshmonk = '/Users/hmatth5/Documents/meshmonk';
addpath(genpath(path2meshmonk));

% add path to facial assessment toolbox on your machine
path2FAT = '/usr/local/avalok/tmp/hmatth5/Chapter5/Submission1/3DGrowthCurves/FacialAssessmentToolbox/';
addpath(genpath(path2FAT));


% save output path
saveOutputPath = '/usr/local/avalok/tmp/hmatth5/Chapter5/Submission1/3DGrowthCurves/FacialAssessmentToolbox/demo';


%% Cell 2 Detect Growth Curves


% determine for which ages models are available
malePath = [path2FAT,filesep,'GrowthCurves',filesep,'Male'];
femalePath = [path2FAT,filesep,'GrowthCurves',filesep,'Female'];

% determine the files that exist that contain a model
maleModelFiles = dir([malePath,filesep,'*.mat']);
femaleModelFiles = dir([femalePath,filesep,'*.mat']);

% remove unreal files beginning with '._' that sometimes appear on mac
maleModelFiles = removeInvisibleMacFiles(maleModelFiles);
femaleModelFiles = removeInvisibleMacFiles(femaleModelFiles);

assert(numel(maleModelFiles)>0,'Unable to find any models for males. Check paths are correct')
assert(numel(femaleModelFiles)>0,'Unable to find any models for females. Check paths are correct')


% detect ages from the file names
% pull floating point numbers from filenames

maleAges = arrayfun(@(x) regexp(x.name,'\d+\.?\d*','match'),maleModelFiles,'UniformOutput',false);
maleAges = arrayfun(@(x) str2num(x{1}{1}),maleAges);

femaleAges = arrayfun(@(x) regexp(x.name,'\d+\.?\d*','match'),femaleModelFiles,'UniformOutput',false);
femaleAges = arrayfun(@(x) str2num(x{1}{1}),femaleAges);
%% Cell 3 Edit Patient Info

% path to patient file
path2patient = '/usr/local/avalok/tmp/hmatth5/Chapter5/Submission1/3DGrowthCurves/FacialAssessmentToolbox/demo/demofaces/demoFace.obj';
% path to patient landmarks
path2patientLandmarks = '/usr/local/avalok/tmp/hmatth5/Chapter5/Submission1/3DGrowthCurves/FacialAssessmentToolbox/demo/demofaces/demoFace.xml';

% specify patient's age and sex
sex = 'M';
age = 35;



assert(exist(path2patient,'file')==2,'unable to find patient image file, verify the file exists and that the path is correct');
[~,fn,ext] =  fileparts(path2patient);

assert(strcmp(ext,'.obj'),'file should be in .obj format, you should be able to convert the format either in the software provided by your camera manufacturers or the free software MeshLab')

assert(exist(path2patientLandmarks,'file')==2,'unable to find patient landmark file, verify the file exists and the path is correct');

assert(ismember(sex,{'M','F'}), 'Specify sex of patient as either ''M'' or ''F'''); 
switch sex
    case 'M'
        modelAges = maleAges;
        modelFiles = maleModelFiles;
    case 'F'
        modelAges = femaleAges;
        modelFiles = femaleModelFiles;
end

% check that patient's age is within the range of the models 
if ~all([age<max(modelAges),age>min(modelAges)])
    warning('Patient''s age is outside the age range of the growth curves');
end
        



%% Cell 4 load template face and landmarks
Template = shape3D;
Template.importWavefront('Template.obj',[path2FAT,filesep,'demo',filesep,'demofaces']);
% read landmrks from xml file
templateLandmarks = readLandmarksMeVisLabXML([path2FAT,filesep,'demo',filesep,'demofaces',filesep,'Template.xml']);

%% Cell 5 load patient and landmarks

[path, name, ext] = fileparts(path2patient);
Patient= shape3D;
Patient.importWavefront([name,ext],path);

patientLandmarks = readLandmarksMeVisLabXML(path2patientLandmarks);


%%
%%%%%%%%%%%%%%%%%%%%%%%% Map patient with MeshMonk non-rigid registration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up shape mappers that will handle the registration

% The ShapeMappers handle the registration according to the given transformation
% type. We define one for the rigid step and one for the non-rigid step
% for this demo we use the default settings we have found optimal for faces
% where the meshes are regular
% see supplementary material of White et al (2019): https://doi.org/10.1038%2Fs41598-019-42533-y

%% Cell 6 Process Images with MeshMonk

% Set up Rigid ICP step using ShapeMapper with 'rigid' transformation type
% % see also the tutorial on the MeshMonk Github page for further info
RM= ShapeMapper;
RM.NumIterations = 30;

RM.TransformationType = 'rigid';
RM.UseScaling = true;

% settings for determining the correspondences
RM.CorrespondencesNumNeighbours = 3; % number of neighbours to use to estimate correspondences
RM.CorrespondencesFlagThreshold = 0.9; 
RM.CorrespondencesSymmetric = true; % if true correspondences are estimated from the template to target and target to template and combined - 
                                    % this can help with mapping structures
                                    % such as long bones and allows the
                                    % target to 'pull' the template
RM.CorrespondencesEqualizePushPull = false;

% settings that determine which points are 'outliers' not used to estimate the
% transformation. This is based on 1. whether the point is an outlier in the
% distribution of differences between floating and target and 2. Whether
% the points corresponds to a point that has been 'flagged'.
RM.InlierKappa = 3;
RM.InlierUseOrientation = true; % use surface normal direction to determine of point is inlier/outlier

% ignore points that correspond to the edges of the mesh - only applicable
% to 'open' surfaces like the face
RM.FlagFloatingBoundary = true; % ignore points that correspondences to the edge of Floating surface - only applicable if using SymmetricCorrespndences
RM.FlagTargetBoundary = true;% ignore points that correspondences to the edge of Target surface

RM.FlagTargetBadlySizedTriangles = true; % ignore points that match to regions of abnormally sized triangles
RM.TriangleSizeZscore = 6; % threshold to determine which triangles are abnormally sized

RM.UpSampleTarget = false; % will upsample the target mesh. if meshes are irregular it can help to set this to true
  

% Set up non rigid ICP with ShapeMapper with 'nonrigid' transformation type
NRM = ShapeMapper;
NRM.TransformationType = 'nonrigid';
NRM.NumIterations = 200; 
NRM.CorrespondencesSymmetric = true;
NRM.CorrespondencesNumNeighbours = 3;
NRM.CorrespondencesFlagThreshold = 0.9;
NRM.CorrespondencesUseOrientation = true;
NRM.CorrespondencesEqualizePushPull =false;
NRM.UseInlierWeights = false; 
NRM.InlierUseOrientation = true;
NRM.FlagFloatingBoundary = true;
NRM.FlagTargetBoundary = true;
NRM.FlagTargetBadlySizedTriangles = true;
NRM.TriangleSizeZscore = 6;
NRM.UpSampleTarget = true;




% parameters specific to non-rigid step. In general, if computing time is not a factor
% set the NumIterations and the TransformNumViscousIterationsStart and
% TransformNumElasticIterationsStart as high as possible. This will ensure
% a slow, gradual deformation of the template to target.
% please see the tutorial 'UnderstandingNonRigidMappingWithMeshmonk' for
% some more details on how this all works. 
% TransformNumViscousIterationsEnd, TransformNumElasticIterationsEnd should
% always be around 1. If they are higher than the floating shape may not
% actually match the shape of the target at the end of the algorithm. They
% can be 0 but this is only advisable with very high quality meshes


NRM.TransformSigma = 3;
NRM.TransformNumViscousIterationsStart = 200;
NRM.TransformNumViscousIterationsEnd = 1;
NRM.TransformNumElasticIterationsStart = 200;
NRM.TransformNumElasticIterationsEnd = 1;
NRM.TransformNumNeighbors = 80;

% step 1 initialise alignments
 %%%% step 1 Initialisation from landmarks

% calculate rigid transformation from templateLandmarks to
% paientLandmarks
T = computeTransform(templateLandmarks,patientLandmarks,true);

% apply rigid transform to Template
forTemplate = applyTransform(Template,T);


%check normals of floating and target consistently point inward or
%outward. If they are inconsistent the normals will be flipped. 
%THIS CHECK WORKS IN MOST CASES BUT IS NOT INFALLIBLE. PARTICULARLY
% IT MAY FAIL IF THE INITIALISATION IS BAD OR TEMPLATE AND PATIENT ARE
% VERY DIFFERENT SHAPES
% The only way to exactly check this is to plot
% each shape with their normals
% e.g v = viewer(Shape)
%    plotVectorField(Shape.Vertices,Shape.VertexNormals,v)

if ~normalsConsistent(Patient, forTemplate)
   Patient.FlipNormals = true;
end

% step 2 refine initialisation with rigid ICP
RM.FloatingShape = clone(forTemplate);
RM.TargetShape = Patient;
RM.map();

NRM.Display = true; % set to false if you don't want to see the mapping
NRM.TargetShape = Patient;
NRM.FloatingShape = clone(RM.FloatingShape);
NRM.map();

mappedPatient = NRM.FloatingShape;

% export mapped image to .obj
mappedPatient.exportWavefront([fn,'_mapped.obj'],saveOutputPath);

%% Cell 7 Inspect correspondences
% load the pre-recorded locations of some landmarks on the template face in
% barycentric co-ordinates
LMs = load([path2FAT,filesep,'matlab',filesep,'data',filesep,'FacialSparseLandmarks.mat']);
bary = LMs.Bary;
index = LMs.Index;
% reconstruct on unmodifed template
LMsOnTemplate = Template.Vertices(index(:,1),:).*bary(:,1) + Template.Vertices(index(:,2),:).*bary(:,2) + Template.Vertices(index(:,3),:).*bary(:,3);
% make shape3D of landmark points
corrLandmarksTemplate = shape3D;
corrLandmarksTemplate.Vertices = LMsOnTemplate;
corrLandmarksTemplate.SingleColor = [1,0,0];
corrLandmarksTemplate.VertexSize = 30;
Template.SingleColor = [0,0,1];
Template.Material = 'Dull';
Template.SingleColor = [.8,.8,.8];
Template.ViewMode = 'SolidWireframe';
close all
v = setUpViewer(Template);
viewer(corrLandmarksTemplate,v);

% reconstruct landmarks in mapping co-ordinate system
LMsOnresult = mappedPatient.Vertices(index(:,1),:).*bary(:,1) + mappedPatient.Vertices(index(:,2),:).*bary(:,2) + mappedPatient.Vertices(index(:,3),:).*bary(:,3);
corrLandmarksPatient = shape3D;
corrLandmarksPatient.Vertices = LMsOnresult;
corrLandmarksPatient.SingleColor = [1,0,0];
corrLandmarksPatient.VertexSize = 30;

% set view properties of obejcts
mappedPatient.ViewMode = 'wireframe';
mappedPatient.SingleColor = [0,0,1];
Patient.ViewMode = 'Solid';
Patient.ColorMode = 'Texture';
Patient.Material = 'Dull';
Patient.SingleColor = [.8,.8,.8];

v2 = setUpViewer(Patient);
viewer(mappedPatient,v2);
viewer(corrLandmarksPatient,v2);

% position second viewer right of the first
v2.Figure.Position(1) = v.Figure.Position(1)+v.Figure.Position(3);

%% Cell 8 Facial Signature
% select which facial signature to view 'X', 'Y', 'Z' or 'Normal';
direction = 'Normal';

% load model that is closest to the age of the patient
[i] = knnsearch(modelAges,age);
fileInfo = modelFiles(i);
M = load([fileInfo.folder,filesep,fileInfo.name]);
M = M.model;

%%%%%%%%%%%%%% facial signature %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% calculate facial signature
facialSignature = computeFacialSignature(M,mappedPatient);
% plot color-maps
% 
%fieldNames = fieldnames(facialSignature);
cLim = [-3,3]; % set limits of color-scale


Zvalues= facialSignature.(direction);
shp = clone(Template);
shp.ViewMode = 'Solid';
shp.Material = 'Dull';
shp.VertexValue = Zvalues;
shp.ColorMode = 'Indexed';
v = setUpViewer(shp);
caxis(v.Figure.CurrentAxes,cLim);
colormap(facialSignatureKULColormap2(1.5,cLim));
cb = colorbar();

% compute signature weight (summary measure of dysmorphism)
sigWeight = sqrt(sum(Zvalues.^2));

% save result
print([saveOutputPath,filesep,fn,'signature_','signatureWeight=',num2str(sigWeight),'.png'],'-dpng');




%% Cell 9 Normal Equivalent
%%%%%%%%%%%   calculate normal equivalent with 'NormalEquivalent' class;
close all
% set up object
NEObj = NormalEquivalent;
NEObj.MorphableModel = clone(M);
NEObj.TargetMesh = clone(mappedPatient);

% make normal equivalent
NEObj.fit();


% view original and normal equivalent
v1 = setUpViewer(NEObj.TargetMesh);
v1.Figure.Position(1) = 50;
NEObj.TargetMesh.ViewMode = 'Solid';
NEObj.TargetMesh.Material = 'Dull';
NEObj.TargetMesh.SingleColor= [.8,.8,.8];
v1.Tag = 'Mapped Original';
v2 = setUpViewer(NEObj.NormalEquivalentMesh);
NEObj.NormalEquivalentMesh.ViewMode = 'Solid';
NEObj.NormalEquivalentMesh.Material = 'Dull';
NEObj.NormalEquivalentMesh.SingleColor= [.8,.8,.8];
v2.Tag = 'Normal Equivalent';
v2.Figure.Position(1) = v1.Figure.Position(1)+v1.Figure.Position(3);

cMapObj = clone(mappedPatient);
cMapObj.ViewMode = 'Solid';
v3 = setUpViewer(cMapObj);
v3.Figure.Position(1) = v2.Figure.Position(1)+v2.Figure.Position(3);

v3.Tag = 'NormalEquivalentDifference (mm)';
diffs = NEObj.NormalEquivalentMesh.Vertices-mappedPatient.Vertices;
%project on surface normals
nProj = sum(diffs.*mappedPatient.VertexNormals,2);

cMapObj.VertexValue = nProj;
cMapObj.ColorMode = 'Indexed';
cMapObj.Material = 'dull';

% set appropriate color-scale limits
cLim = [-5,5]; % mm units
caxis(v3.Figure.CurrentAxes,cLim);

colormap(seismicColormap());
colorbar();


% export normal equivalent to obj
mappedPatient.exportWavefront([fn,'_normalEquivalent.obj'],saveOutputPath);
