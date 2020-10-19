function [out] = batchResidualiseNonlinearLOO(growthCurve3DObj,landmarks,methods,normalEquivalent,RRT,returnCasesMask)
    % in batch residualises the shapes in landmarks against an age and sex appropriate model, as produced by the growth curve 3d
    % this assumes that the shapes in landmarks are the same individuals as
    % as for the growthCurve3D training data.  Each case is left out of the
    % training data used to compute the model against which it is
    % residualised.
    % INPUTS:
    %       growthCurve3DObj - a growthCurve3D object trained using the sam
    %                           subjects as in 'landmarks'
    %       landmarks - an n (vertices) x 3 (x,y,z) co-ordinates x k shapes x l repeat measurements array
    %               of shapes
    %       methods - is passed to 'batchResidualiseNonlinear' (see documentation of this function)
    %       normalEquivalent - is passed to 'batchResidualiseNonlinear' (see documentation of this function)
    %       RRT - is passed to 'batchResidualiseNonlinear' (see documentation of this function)
    %       returnCasesMask - a true/false array indicating which cases to
    %       return results for (true) or skip (false)
    % OUTPUTS:
    %       out - residuals for all cases for all repeat measurements.
    % copyright Harold Matthews harry.matthews@kuleuven.be (2020)
    if nargin<2 || isempty(landmarks)
        landmarks = growthCurve3DObj.TrainingLandmarks;
    end
    % these empty arguments will e later passed into
    % batchResidualiseNonlinear which will set the defaults
    if nargin<3 || isempty(methods)
        methods = [];
    end
    if nargin<4 || isempty(normalEquivalent)
        normalEquivalent = [];
    end
    if nargin <5 || isempty(RRT)
        RRT = [];
    end
    if nargin<6 || isempty(returnCasesMask)
       returnCasesMask =  true(numel(growthCurve3DObj.Age),1);
       
    end
    assert(size(landmarks,3)==numel(growthCurve3DObj.Age))

    growthCurve3DObj = clone(growthCurve3DObj);
    growthCurve3DObj.TrainingLandmarks = []; % travel light in parfor loop
    growthCurve3DObj.clearToBasicSettings();
    growthCurve3DObj.GlobalAverageShape = [];
    landmarks = landmarks(:,:,returnCasesMask,:);
    nAllCases = numel(growthCurve3DObj.Age);
    nCasesToReturn = size(landmarks,3);
    tmpResults = cell(1,nCasesToReturn);
  
    returnCasesInds = find(returnCasesMask);
    parfor n = 1:nCasesToReturn
        
        mask = true(nAllCases,1);
        mask(returnCasesInds(n)) = false;
        forGrowthCurve3DObj = clone(growthCurve3DObj);
        forGrowthCurve3DObj.UserExclusionMask = mask;
        shp = landmarks(:,:,n,:);
        age = forGrowthCurve3DObj.Age(returnCasesInds(n));
        sex = forGrowthCurve3DObj.Sex(returnCasesInds(n));
        tmpResults{n} = batchResidualiseNonlinear(shp,age,sex,forGrowthCurve3DObj,methods,normalEquivalent,RRT,false);

    end
    out = vertcat(tmpResults{:});
end

