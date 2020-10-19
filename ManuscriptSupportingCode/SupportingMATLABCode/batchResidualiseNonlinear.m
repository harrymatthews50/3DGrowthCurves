function [residuals,success, msg] = batchResidualiseNonlinear(shapes,ages,sexes,growthCurve3DObj,methods,normalEquivalent,RRT,useGlobalModel)
%%% for each observation in shapes this residualises the shape against an age and sex-specific model of normal vraiation as as evaluated by the 3DGrowthCurve object
% INPUTS:
%       shapes -an n (vertices) x 3 (x,y,z) co-ordinates x k shapes x l repeat measurements array
%       of shapes
%       ages - a numeric vector of the ages corresponding to each shape
%       sexes - a numeric vector of the ages corresponding to each shape
%       growthCurve3DObj- a growthCurve3D object that will evaluate the
%           models of normal variation, based on its training data and
%           settings
%       methods - the shape can be residualised agins the model in
%                different ways. This is a cell array of strings, each
%                indicating a method (default = {'normalEquivalent','dist2Expected','facialSignature'})
%                options are: 
%                   'normalEquivalent' - residuals are the distances from
%                               the normal equivalent to the subject
%                   'facialSignature' - residuals are z-scores for each point on the face for the subject 
%                   'dist2Expected'  - residuals are the distance between
%                               the subject and the expected face computed
%                   'project2Model' - residuals are the distances between
%                           the subject and their reconstruction from the modes of
%                           variation of the model
%      normalEquivalent - a normalEquivalent object (the settings of this object will determine the settings used to compute normal equivalents)
%      RRT - a 'RobustRigidTransform' object that will control the
%           settings for aligning the shape to the model for method 'facialSignature' and 'dist2Expected' 
%      UseGlobalModel -if true all settings of the growthCurve3D object will
%           be ignored and the model will be based on all of the training data -
%           this is intended to be used only for comparison
%
%       OUTPUTS:
%           residuals - the residuals for each subject in a cell array
%           success - a boolean array indicating (for each subject) if they
%               were processed without error
%           msg - contains the error messages for subjects for whom success(i)==false.           

% set defaults if arguments not specified
if nargin<8 || isempty(useGlobalModel)
    useGlobalModel = false;
end
if nargin<7 || isempty(RRT)
    RRT = RobustRigidTransform;
end
if nargin<6 || isempty(normalEquivalent)
    normalEquivalent = NormalEquivalent; % use defaults
end
if nargin<5 || isempty(methods)
    methods = {'normalEquivalent','dist2Expected','facialSignature'};
end
nShapes = size(shapes,3);
nReps = size(shapes,4);
assert(numel(ages)==nShapes,'Wrong number of ages');
assert(numel(sexes)==nShapes,'Wrong number of sexes');

growthCurve3DObj = clone(growthCurve3DObj);
%travel light into parfor loop
growthCurve3DObj.clearToBasicSettings();
growthCurve3DObj.TrainingLandmarks = [];
growthCurve3DObj.GlobalAverageShape = [];


success = true(nShapes,1);
msg = cell(nShapes,1);
residuals =cell(nShapes,1);
if useGlobalModel
    growthCurve3DObj.computeGlobalVariationModel();
end
for i = 1:nShapes
   % i
    try
        forGrowthCurve3D = clone(growthCurve3DObj);
       if ~useGlobalModel % instead compute age-specific model (default)
            forGrowthCurve3D.EvaluationAge = ages(i);
            forGrowthCurve3D.EvaluationSex = sexes(i);
       end
        outRow = cell(1,nReps);
        for r = 1:nReps
            outElement = struct;
            testShape = shapes(:,:,i,r);
            shapeObj=clone(growthCurve3DObj.RefScan);
            shapeObj.Vertices = testShape;

            for m = 1:numel(methods)
               methodName = methods{m};
               switch methodName
                   case 'normalEquivalent'
                       outElement.(methodName) = residualizeNormalEquivalent(clone(forGrowthCurve3D.VariationModel),normalEquivalent,testShape);     
                   case 'dist2Expected'
                       outElement.(methodName) = residualizeDist2Expected(forGrowthCurve3D.ExpectedShapeObj,shapeObj,RRT,forGrowthCurve3D.VariationModel.PointStandardDevs(:,5));
                   case 'facialSignature'
                       outElement.(methodName) = residualizeFacialSignature(shapeObj,forGrowthCurve3D.VariationModel,RRT);
                   case 'project2Model'
                       outElement.(methodName) = residualizeProject2Model(clone(forGrowthCurve3D.VariationModel),testShape);
                   otherwise
                       error('Invalid method name');
               end
            end
            outRow{r} = outElement;
        end
        residuals{i} = [outRow{:}];
     


    catch me
        success(i) = false;
        msg{i} = me;
        disp(me);
    end

    
    
end
%out = vertcat(out{:});


end


    



function out=residualizeNormalEquivalent(model,NEObj,shape)
    NEObj.Weights = [];
    NEObj.MorphableModel = clone(model);
    
    Target = clone(model.Average);
    Target.Vertices = shape;
    
    NEObj.TargetMesh = Target;
    NEObj.fit();
    normalEquivalent = NEObj.NormalEquivalentMesh.Vertices;
    weights = NEObj.Weights;
    
    diffs = shape-normalEquivalent;
    
    residuals = struct('XYZDisplacements',[],'X',[], 'Y',[],'Z',[],'Normal',[],'Total',[]);
    residuals.XYZDisplacements = diffs;
    residuals.X = diffs(:,1);
    residuals.Y = diffs(:,2);
    residuals.Z = diffs(:,3);
    
    normals = NEObj.NormalEquivalentMesh.VertexNormals;
    normProj = sum(normals.*diffs,2);
    residuals.Normal = normProj;
    residuals.Total = sqrt(sum(diffs.^2,2));
    out = struct;
    out.Residuals = residuals;
    out.AdditionalInfo = struct;
    out.AdditionalInfo.NormalEquivalent = normalEquivalent;
    out.AdditionalInfo.Weights = weights;
 



end

function out = residualizeProject2Model(model,shape)
    [p,MD,sc,diffs] = model.computeOutOfSamplePValue(shape);
    residuals = struct('XYZDisplacements',[],'X',[], 'Y',[],'Z',[],'Normal',[],'Total',[]);
    residuals.XYZDisplacements = diffs;
    residuals.X = diffs(:,1);
    residuals.Y = diffs(:,2);
    residuals.Z = diffs(:,3);
    
    residuals.Total = sqrt(sum(diffs.^2,2));
    out = struct;
    out.Residuals = residuals;
    out.AdditionalInfo = struct;
    out.AdditionalInfo.p = p;
    out.AdditionalInfo.MD = MD;
    out.AdditionalInfo.Scores=sc;
    


end


function out = residualizeFacialSignature(shape,model,RRT)
    RRT.Weights = [];
   % RRT.ModelRMSDevs = model.PointStandardDevs(:,5);
    RRT.FloatingMesh = clone(shape);
    RRT.TargetMesh = clone(model.Average);
    RRT.fit();
    diffs = RRT.FloatingMesh.Vertices-RRT.TargetMesh.Vertices;
    residuals.XYZDisplacements = diffs;
    
    Zscore = diffs./model.PointStandardDevs(:,1:3);
    residuals.X = Zscore(:,1);
    residuals.Y = Zscore(:,2);
    residuals.Z = Zscore(:,3);
    normals = RRT.TargetMesh.VertexNormals;
    normProj = sum(normals.*diffs,2);
    residuals.Normal = normProj./model.PointStandardDevs(:,4);
    residuals.Total = sqrt(sum(diffs.^2,2))./model.PointStandardDevs(:,5);
    out.Residuals = residuals;
    out.AdditionalInfo = struct;

    out.AdditionalInfo.Weights = RRT.Weights;


    
    

end


function out = residualizeDist2Expected(expected,shape,RRT,expectedModelDevs)
   RRT.Weights = [];
  % RRT.UseModelExpectedDevs = true;
   RRT.ModelRMSDevs = expectedModelDevs;
   RRT.FloatingMesh = clone(shape);
   RRT.TargetMesh = clone(expected);
   RRT.fit();
   diffs = RRT.FloatingMesh.Vertices-RRT.TargetMesh.Vertices;
    
    residuals = struct('XYZDisplacements',[],'X',[], 'Y',[],'Z',[],'Normal',[],'Total',[]);
    residuals.XYZDisplacements = diffs;
    residuals.X = diffs(:,1);
    residuals.Y = diffs(:,2);
    residuals.Z = diffs(:,3);
    
    normals = RRT.TargetMesh.VertexNormals;
    normProj = sum(normals.*diffs,2);
    residuals.Normal = normProj;
    residuals.Total = sqrt(sum(diffs.^2,2));
    out = struct;
    out.Residuals = residuals;
    out.AdditionalInfo = struct;

    out.AdditionalInfo.Weights = RRT.Weights;
    out.AdditionalInfo.ExpectedShape = expected.Vertices;
   


end

