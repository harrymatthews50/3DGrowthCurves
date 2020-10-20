%% This script takes the reconstruction errors calculated in 'Scr1' and uses these to estimate the optimal bandwidth for each age
% bandwidth as a function of age is stored separately for each sex in a griddedInterpolant object. 

% to run, this also requires a toolbox for scattered data interpolation
% from the MATLAB file exchange

%https://www.mathworks.com/matlabcentral/fileexchange/10056-scattered-data-interpolation-and-approximation-using-radial-base-functions


%% set MATLAB paths
clear all; close all; restoredefaultpath;
% add path to supporting matlab code
addpath(genpath('/usr/local/avalok/tmp/hmatth5/Chapter5/Submission1/3DGrowthCurves/ManuscriptSupportingCode/SupportingMATLABCode'));

% add path to Patient Assessment Toolbox
addpath(genpath('/usr/local/avalok/tmp/hmatth5/Chapter5/Submission1/3DGrowthCurves/PatientAssessmenToolbox/matlab'));

% add path to MeshMonk
addpath(genpath('/usr/local/avalok/tmp/hmatth5/Projects/meshmonk'));


% add path to 'rbfinterp' toolbox from the MATLAB file exchange
addpath(genpath('/usr/local/avalok/tmp/hmatth5/Chapter5/Submission1/rbfinterp_v1'));
%% load reconstruction errors estimated in the previous script
% get list of files

src = '/usr/local/avalok/tmp/hmatth5/Chapter5/Submission1/SimDataJobs';
files = dir([src,filesep,'*.mat']);

% remove files preceded by '._' that appear on mac sometimes
files = removeInvisibleMacFiles(files);
testBWs = strrep({files.name},'.mat','');
testBWs = cellfun(@str2num,testBWs);

% sort ascending
[testBWs,inds] = sort(testBWs);
files = files(inds);

% load results
results = arrayfun(@(x) load([x.folder,filesep,x.name]),files,'UniformOutput',false);


%% load growth curve metadata

load('/usr/local/avalok/tmp/hmatth5/Chapter5/Submission1/3DGrowthCurves/ManuscriptSupportingCode/DATA/Simulated_Normative_Data/SimulatedNormativeData.mat','metadata')

%%

% compute generalisation error for each subject for each bandwidth

genError = cellfun(@(x) squeeze(mean(sqrt(sum(x.XYZResiduals.^2,2)),1)),results,'UniformOutput',false);
genError = horzcat(genError{:});





% create matrices of metadata for each entry in genError by repeating
% metadata for each row of genError (corresponding to a particular bandwidth)
ages = repmat(metadata.age,[1,size(genError,2)]);
sex = repmat(metadata.sexNumeric,[1,size(genError,2)]);
BWs = repmat(testBWs,[size(genError,1),1]);

% ravel to vectors
genError = genError(:);
ages = ages(:);
sex = sex(:);
BWs = BWs(:);
%% 

% remove any nan (from where the reconstruction error was not calculated
% due to very small amounts of training data) % in paper this only ocurred
% for cases aged > 70 years

mask = ~isnan(genError);
genError = genError(mask);
ages = ages(mask);
sex = sex(mask);
BWs = BWs(mask);





%% Create interpolants for each sex modelling change in generalisation error as a function of age and bandwidth
% this can take a long time 

ErrorRBFs = cell(1,2);
for s = 1:2
    X = [ages,BWs,genError];
    X = X(sex==s,:);
    ErrorRBFs{s}= rbfcreate(X(:,1:2)',X(:,3)','RBFFunction','linear','RBFSmooth',10);
end


%%

f = figure;
threshold = 1.; % set acceptable average generalization error (mm), in the manuscript we set this to .6 but is set higher here because the simulated data don't give generalsation errors that low

for s = 1:2 % for each sex
    
    %%%%%%%%%%%%%% plot scattered data and surface fitted to scattered data
    ax = subplot(2,2,s);
    
    X = [ages,BWs,genError];
    X = X(sex==s,:);
    
    
    % plot scattered data
    scatter3(X(:,1),X(:,2),X(:,3),1,'filled','MarkerEdgeColor',[1,0,0]);

    hold on
    rbf = ErrorRBFs{s};
    
    % work out range over which to visulaise the rbf interpolants
    interpX = linspace(min(X(:,1)),max(X(:,1)),60);
    interpY = linspace(min(X(:,2)),max(X(:,2)),60);
    [XI,YI] = meshgrid(interpX,interpY);
    
    % interpolate generalisation error at the selected points
    ZI = rbfinterp([XI(:)';YI(:)'],rbf);
    ZI = reshape(ZI,size(XI));
    
    % plot fitted surface
    surf(XI,YI,ZI);
    
    %colorbar and color-scale
    cb = colorbar();
    caxis([0,2]);
    t = get(cb,'Title');
    t.Interpreter = 'latex';
    t.String = 'Error (mm)';
    t.FontSize = 8;
    
    %axis labels 
    xlabel('Age');
    ylabel('Bandwidth');
    zlabel('Error (mm)','Interpreter','latex');
    set(ax,'FontSize',6);
    if s==1
        label = 'Male';
    elseif s==2
        label = 'Female';
    end
    title(label)
    hold off
    
    
  %%% plot isoline at level of acceptable generalisation error 

    ax = subplot(2,2,s+2);

    hold on
    
    % plot elevation map
    [~,contf] = contourf(XI,YI,ZI,100,'EdgeColor','none');
    xlabel('Age');
    ylabel('Bandwidth');
    colormap('default');
    caxis([0,2]);

    
    % find isoline tracing bandwidth as a function of age at the level of
    % the acceptable generalization error
    
    c = contour(XI,YI,ZI,[threshold,threshold]);%,'EdgeColor','None');
    
    % extract only the x(age)  and y(bandwidth) co-ordinates from the resulting contour
    % matrix
    xy = trimContourMatrix2Vertices(c);
    
    
    
    
    % sort by ascending age
    [~,I] = sort(xy(1,:),'ascend');
    xy = xy(:,I);
   % keyboard
    % create interpolant predicting bandwidth from age
    interA = griddedInterpolant(xy(1,:),xy(2,:));
    interA.ExtrapolationMethod = 'nearest';
    
    
    
    
    
    % use this to interpolate any missing sections of the contour
    evalAges = linspace(min(ages(sex==s)),max(ages(sex==s)),100);
    interpBW = inter(evalAges);
    
    % create final interpolant including missing contour sections
    inter = griddedInterpolant(evalAges,interpBW);
    
    % plot bandwidth in function of age as modelled in the interpolant

    pl = plot(evalAges,interpBW,'r','LineWidth',3);
    legend(pl,num2str(threshold),'Location','NorthWest')
    hold off
    
    % assign interpolant to appropriate variable name
    if s==1
        maleInterp = inter;
    elseif s==2
        femaleInterp = inter;
    end

end
% save interpolants for use in next script

save('/usr/local/avalok/tmp/hmatth5/Chapter5/Submission1/BandwidthInterpolants/BandwidthInterpolants.mat','maleInterp','femaleInterp');



%%

function out = trimContourMatrix2Vertices(in)
    inds = [];
    i = 1;
    while true
       nPointsinContour = in(2,i);
      
       inds = [inds,(i+1):(i+nPointsinContour)];
       i=(i+1+nPointsinContour);
       if i>size(in,2)
           break
       end
    end
    out = in(:,inds);
           
        
        
end







