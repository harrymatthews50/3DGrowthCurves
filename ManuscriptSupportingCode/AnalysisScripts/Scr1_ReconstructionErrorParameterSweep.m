%%% This script runs a parameter sweep testing different bandiwdths of the weighting function and computing the
% generalisation of the models in function of age.  This is used later in 'Scr_2' to adaptively determine the correct bandwidth to use for a given age
% see S1 Text section 2 of the manuscript for more information.

%This is demonstrated using simulated data, as such the results will not be

%%%%%%%%%%%%%%%%%% To create growth curves from your own data
%%% 1) map each scan with MeshMonk
%%% 2) load the resulting shapes into an n vertices x 3 dimensions x k
%%% scans array, called 'shapes'
%%% 3) create a table called 'metadata' with k rows, corresponding to the third dimension of
%%% 'shapes'. This should have 1 column called 'age' containing the age of
%%% the subject and another column called 'sex' containing te numerically
%%% coded sex of the subject (1=male; 2=female)




clear all; close all; restoredefaultpath;
% add path to supporting matlab code
addpath(genpath('/usr/local/avalok/tmp/hmatth5/Chapter5/Submission1/3DGrowthCurves/ManuscriptSupportingCode/SupportingMATLABCode'));

% add path to Patient Assessment Toolbox
addpath(genpath('/usr/local/avalok/tmp/hmatth5/Chapter5/Submission1/3DGrowthCurves/PatientAssessmenToolbox/matlab'));

% add path to MeshMonk
addpath(genpath('/usr/local/avalok/tmp/hmatth5/Projects/meshmonk'));

%% If running to test code will use only a subset of the data and only test a small number of bandwidths
test = false;


%% load data and build growth curve of normal data


GCTrain = load('/usr/local/avalok/tmp/hmatth5/Chapter5/Submission1/3DGrowthCurves/ManuscriptSupportingCode/DATA/Simulated_Normative_Data/SimulatedNormativeData.mat');

%load the template obj from MeshMonk toolbox
load('/usr/local/avalok/tmp/hmatth5/Projects/meshmonk/demo/Template.mat');
%%
% set some properties
GC = growthCurve3D;
GC.TrainingLandmarks = GCTrain.shapes;
GC.Age = GCTrain.metadata.age;
GC.Sex = GCTrain.metadata.sexNumeric;
GC.ScaleShapes = true; % make growth curves of shape not form
GC.RefScan = clone(Template);

GC.TrainingLandmarks = GCTrain.shapes;
originalShapes = GC.TrainingLandmarks;
GC.Age = GCTrain.metadata.age;
GC.Sex = GCTrain.metadata.sexNumeric;
GC.Speedy = false;
% end
%% Initialise
GC.initialize();
%%
GC.TrainingLandmarks = [];
testBandwidths = linspace(1,15,15);

minTrainValues = 10;


%% set locations for saving output
dstPath = '/uz/data/avalok/mic/tmp/hmatth5/Chapter5/Submission1/SimDataJobs';

if ~exist(dstPath,'dir')
    mkdir(dstPath)
end

recomputeExistingResults = false;% if results already exist for some bandwidths then skip

existingFiles = dir([dstPath,filesep,'*.mat']);
if recomputeExistingResults
    for f = 1:numel(existingFiles)
        delete([existingFiles(f).folder,filesep,existingFiles(f).name])
    end
else
    executedJobs = strrep({existingFiles.name},'.mat','');
    executedJobs = cellfun(@str2num,executedJobs);
    testBandwidths = setdiff(testBandwidths,executedJobs);
end

        
       

nBW = numel(testBandwidths);



cellResults = cell(1,nBW);
for i = 1:nBW
   forGC = clone(GC);
   forGC.BWSigma = testBandwidths(i);
   disp(['computing for bandwidth',num2str(forGC.BWSigma)]);
   
   % work out for which cases there is not enough training data around the
   % age to bother computing the generalization error (as it will be very
   % we start calculating when the number of training cases is withing +/-3
   % BW of the target age is greater than 10
   
   % number of training cases within range each observation
   nTrainValues =  nWithinRange(forGC.Sex,forGC.Age,forGC.Sex,forGC.Age,forGC.BWSigma*forGC.MaxZ);
   mask = nTrainValues>minTrainValues;
   forResults = cell(numel(mask),1);

   forResults(mask) = batchResidualiseNonlinearLOO(forGC,originalShapes,{'project2Model'},[],[],mask);

    



   XYZResiduals = getElementsFromResidualsCell(forResults,'project2Model','XYZDisplacements');
   save([dstPath,filesep,num2str(testBandwidths(i)),'.mat'],'XYZResiduals')

end
%%


function nIn = nWithinRange(testCat,testValues,trainCat,trainValues,range);
nIn = zeros(1,numel(testValues));
for i = 1:numel(testValues)
    a = testValues(i);
    catMask = trainCat==testCat(i);
    inRangeMask=(trainValues>=(a-range)).*(trainValues<=(a+range));
    nIn(i) = sum(inRangeMask.*catMask);
end
    
end


