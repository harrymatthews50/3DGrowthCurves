classdef morphableShape3D < superHandleClass
%   morphableShape3D - fits and represents a statistical shape model
%               morphableShape3D also supports weighted and regularized
%               projection onto the modes of variation
%       
%   morphableShape3D Properties
%           %% User settings
%           ProjectionRegularisationType - controls the type of
%               regularisation used when projecting a new face onto the
%               model when calling 'compute_morphable_transformation' or 
%               'projectIntoSpace'. See below 'Description of
%               regularization types' for options and explanations.
%           ProjectionRegularisationValue - controls the amount of
%               regularisation used when projecting a new face onto the
%               model when calling 'compute_morphable_transformation' or 
%               'projectIntoSpace'. See below 'Description of
%               regularization types' for options and explanations.  
%            ModelTrimmingType - 'PercentageVariance' controls how PCs are
%               removed from the model when 'stripModel' is called. Options
%               are 'None' or 'PercentageVariance' (default). If 'None' then
%               no PCs are removed if 'PercentageVariance' then PCs are
%               removed that explain up to a set amount of variation
%            ModelTrimmingArgs - A cell array of arguments that control how many PCs are
%               retained:
%               if ModelTrimmingType = 'PercentageVariance' set this to contain
%               the percentage of variance you want retained (default={98}) 
%          
%           %% Model Properties
%           Average - a shape3D object representing the model center shape
%               (by default this will be the average but it can also be 
%               specified by the user when calling fitMorphableModel)
%           EigVec - The modes of variation of the morphable model, each
%               column of this matrix is a mode.
%           EigVal - the variance of the projections/scores along each mode
%                of variation.
%           Projections - The projections/scores of observations used to
%                train the model. Rows correspond to observations, 
%                columns to modes of variation. These are removed when 
%                models are publicly released.
%           PointStandardDevs - an n points x 5 matrix of the root mean
%               squared displacements of the training shapes from obj.Average
%               shape. These are calculated along the 'X' (first column),
%               'Y' (second column) and 'Z' (third column) axes and also
%               along the direction of the corresponding surface normal of
%               obj.Average (fourth column) and also the total root mean
%               square displacement (fifth column)
%         PointStandardDevsDescription - a text deescription of the columns
%            of PointStandardDevs. Is always {'X','Y','Z','Normal','TotalRMS'}
%         Dim - the number of dimensions (modes) of the model
%         EigStd - the standard deviation of projections/scores along each
%                  mode. Is sqrt(EigVal).
%        PercVarianceExplained - the percentage of variance explained by
%                   each mode
%        OK2Share - if true then obj.Projections is empty and the model can
%        be widely shared.
%        MDPerObservation - the Mahalanobis distance of each training
%               observation from the origin.
%        PValuesPerObservation - a p-value for each training observation
%               according to the multinormal pdf modelled by the modes and EigVal 
%
%    Description of regularisation types
%       'none' - no regularisation is applied
%       'neta' - using this the projection minimises a combined cost of 
%               both reconstruction error and mahalanobis distance of the 
%               result from the origin (Eq 25 of Blanz et al (2004);doi: 10.1109/TDPVT.2004.1335212 )
%               Sensible corresponding values of obj.ProjectionRegularisationValue​​ are between 0-1. This 
%               value is multiplied by the first eigenvalue of the model to 
%               determine a coefficient of the cost function that varies 
%               the penalty associated with the mahalanobis distance from 
%               the origin. High values equals more regularisation.
% 
%       'nopcs' - Performs the projection on only the first N pcs. 
%               Sensible values of ProjectionRegularisationValue are integers 
%               between 1:obj.Dim.
%           
%       'mahalanobis' if the projection exceeds a specified critical 
%               mahalanobis distance the vector of scores is uniformly 
%               scaled so that its mahalanobis distance is equal to the critical value.
%               sensible values ProjectionRegularisationValuedepend on the dimensionality of the model, 
%               are roughly 0-sqrt(chi2inv(.99,obj.Dim));
%               For a nice explanation Mahalanobis distance and alos its 
%               relationship to principal components analysis and the 
%               Chi-squared statistic see Brereton (2015a; doi: 10.1002/cem.2692, 2015b; doi:
%               10.1002/cem.2680)
%           
%        'pvalue' the same idea as 'mahalanobis'. In this case the critical
%               mahalanobis distance is set as a p-value by the user and then 
%               converted into the critical mahalanobis distance: 
%               sqrt(chi2inv(1-p,obj.Dim)).
%               sensible values are between 0-1.
%
%        'clamp' the score/projection (p) onto a given mode (mode number l) 
%               falls outside of a specified number of standard deviations (s) the projection is 
%               clamped to sign(p)*obj.EigStd(l)*s. Sensible values are (obviously) 0-3.
%
% Written by Peter Claes (peter.claes@kuleuven.be) and Harold Matthews (harry.matthews@kuleuven.be)
        
    properties % user settings
        ProjectionRegularisationType = 'None';
        ProjectionRegularisationValue = [];
        ModelTrimmingType = 'PercentageVariance';
        ModelTrimmingArgs = {98.};
    end


    properties  % model properties
        Average;
        EigVec;
        EigVal;
        Projections;
        PointStandardDevs;
        
    end
    properties (Constant=true)
        PointStandardDevsDescription = {'X','Y','Z','Normal','TotalRMS'};
    end
    properties (Dependent = true) % model dependent properties
        Dim;
        EigStd;
        PercVarianceExplained;
        OK2Share;
        MDPerObservation;
        PValuesPerObservation;
    end
    properties (Hidden = true)
        AlignmentTransformation;
        TotalInitialVariance; %  the total variance present in the model when fitted
        
    end
    properties (Hidden = true, Dependent = true)
        AverageVEC;
    end
    methods %CONSTRUCTOR
        function obj = morphableShape3D(varargin)
            obj = obj@superHandleClass(varargin{:});
        end
    end 
    methods %GETTING
        function out = get.AverageVEC(obj)
            if isempty(obj.Average), out = []; return; end
            out = obj.Average.Vertices';out = out(:);
        end
        function out = get.Average(obj)
            out = obj.Average;
            if ~superHandleClass.isH(out), out = []; end
        end
        function out = get.Dim(obj)
            if isempty(obj.EigVal), out = []; return;end
            out = length(obj.EigVal);
        end
        function out = get.EigStd(obj)
            out = sqrt(obj.EigVal);
        end
        function out = get.PValuesPerObservation(obj)
             if isempty(obj.MDPerObservation), out = []; return;end
             out = 1-chi2cdf(obj.MDPerObservation.^2,obj.Dim);
        end
        function out = get.MDPerObservation(obj)
             if isempty(obj.Projections), out = []; return;end
                out = sqrt(sum((obj.Projections./obj.EigStd').^2,2));  
        end
        function out = get.PercVarianceExplained(obj)
                out = (cumsum(obj.EigVal)./obj.TotalInitialVariance)*100;
        end
        function out = get.OK2Share(obj)
            if isempty(obj.Projections)
                out = true;
            else
                out = false;
            end
        end
        
    end
    
    methods %SETTING
        function obj = set.ProjectionRegularisationType(obj,regtype)
            
            validtypes = {'neta','clamp','nopcs','pvalue','mahalanobis', 'none'};
            % check if a valid value

            if ~sum(strcmp(lower(regtype), validtypes))==1
                error(strcat('Invalid regularisation type. Must be empty or one of Neta, Clamp, NoPCs,Malanobis, Pvalue or None'));
            end
            
            obj.ProjectionRegularisationType = regtype;
            
        end
        
        function obj = set.ModelTrimmingArgs(obj,varargin)
            if iscell(varargin{:})
                obj.ModelTrimmingArgs = varargin{:};
            else
                obj.ModelTrimmingArgs = {varargin{:}};
            end
        end
        
    end
        methods %FOR FITTING
        function fitMorphableModel(obj,vertices,faces,center,avgShape,weights,normaliseVariationFactor)
%            fitMorphableModel - fits the morphable mode to the training
%               Performs a (weighted) Principal Components Analysis of the
%               shapes.
%            INPUTS: any input argument (except vertices) can be
%            substituted with [] to skip or revert to default values.
%                 vertices - either a k (vertices) x 3 x n (shapes) array
%                       of co-ordinates or an n (shapes) x 3k matrix or
%                       co-ordinates to which to fit the model. If the
%                       latter, it is assumed that each row is ordered:
%                       [x vertex 1, y vertex 1, z vertex 1 ... x vertex k, y vertex k, z vertex k]
%                 faces - an l (faces) x 3 (vertices) array specifiying 
%                       the connectivity of the points. This is not essential
%                       for fitting the model but visualisations will not work without it.
%                 center - if true  the verticies will be centered either on
%                       the (weighted if weights are specified) mean shape
%                       or (if a user-specified avgShape is specified) on
%                       the 'avgShape' specified.
%                avgShape - optional user-specified shape on
%                       which to center the model, should be either an
%                       instance of the shape3D class or a k(vertices) x 3
%                       array of co-ordinates. if center=false
%                weights - an n(shapes) x 1 array of weights. elements
%                       correspond to shapes in 'vertices'. Shapes with high
%                       weights have more influenece on the model than shapes with
%                       low weights. defaults to all ones, having no
%                       effect.
%                normaliseVariationFactor - the variances 
%                       will be normalised by dividing by this factor. This is useful
%                       to specify if weights have been applied to 'vertices' outside of this function. 
%                       If not specified this defaults to the sum of the weights. 'weights'
%                       and 'normaliseVariationFactor' are not intended to be
%                       specified together and one or the other should always
%                       be left to default.
% 
            if numel(size(vertices))==3 %reshape to 2D array
                vertices = arrayStructure2Vector(vertices);
            end
            
          
            % check for incompatible arguments
            if nargin>6 
                if ~isempty(weights) && ~isempty(normaliseVariationFactor)
                    error('Observation weights and variation normalisation are not intended to be applied together')
                end
            end
           
             

            if nargin<6 || isempty(weights)
                % dummy values
                weights = ones(size(vertices,1),1);
            end
            if nargin<7 || isempty(normaliseVariationFactor)
                normaliseVariationFactor= sum(weights);
            end
            
            %%%%%Centering
            if nargin<4 || isempty(center)
                center = true;
            end
 
            if nargin<5 || isempty(avgShape)
             % create average obj
             avgShape = shape3D;
             avgShape.Vertices = arrayVector2Structure(weightedColMean(vertices,weights));
            else
                if isobject(avgShape)
                    avgShape = clone(avgShape);
                else
                    tmp = avgShape;
                    avgShape = shape3D;
                    avgShape.Vertices = tmp;
                end
            end
            if nargin>2 && ~isempty(faces)
                avgShape.Faces = faces;
            else
                warning('No faces specified for the meshes so averages and reconstructed shape instances will not be returned correctly')
            end   
                       
            obj.Average = avgShape;
            
            %center and do svd
            if center
                vertices = bsxfun(@minus, vertices, obj.AverageVEC');
            end
            
            % do svd on weighted vertices
            [~,S,V] = svd(bsxfun(@times,vertices,weights),'econ');
            obj.EigVec = V;
            % project unweighted vertices onto eigvect
            obj.Projections = vertices*obj.EigVec;
            
            % calculate the weighted variance of the projections along each
            % PC
            
            %compute weighted SS (if no weights have been specified this is the num of squares)
            wSS = sum((obj.Projections.^2.*weights),1);
            eV = wSS./normaliseVariationFactor;

            
            
            obj.EigVal = eV'; 
            obj.TotalInitialVariance = sum(obj.EigVal); % this will be used to deal with multiple calls of stripPercVar
            
            % reshape and claculate pointwise standard deviations
            vertsStruc = arrayVector2Structure(vertices);
            
            pointStandardDevs = nan(size(vertsStruc,1),5);
            % normal projections
            if ~ isempty(avgShape.Faces)
                normals = avgShape.VertexNormals;
                normProjs = squeeze(sum(normals.*vertsStruc,2));
            else
                normProjs = [];
            end
            
            displacements = {squeeze(vertsStruc(:,1,:)),squeeze(vertsStruc(:,2,:)), squeeze(vertsStruc(:,3,:)),normProjs, squeeze(sqrt(sum(vertsStruc.^2,2)))};
            for d = 1:numel(displacements)
               vals = displacements{d}; 
               if ~isempty(vals)
                    % calculate weighted SS for each point
                     wSS = sum(vals.^2.*weights',2);
                     pointStandardDevs(:,d) = sqrt(wSS./normaliseVariationFactor);

               end
            end
            obj.PointStandardDevs = pointStandardDevs;

               
             
            
                        
            
            obj.stripModel();
        end
        function obj = stripModel(obj)
%               stripModel - removes modes of variation from the model
%                   the bavhiour here is controlled by obj.ModelTrimmingType
%                   and obj.ModelTrimmingArgs.
            switch obj.ModelTrimmingType
                case 'PercentageVariance'
                    obj.stripPercVar(obj.ModelTrimmingArgs{:});
                case 'None'
                    % do nothing
                otherwise
                    error(strcat('ModelTrimmingType ', obj.ModelTrimmingType,'not implemented or not valid'))
            end
         end 
      
        function obj = stripPercVar(obj,perc)
%           stripPercVar - removes modes after a certain threshold of variance is explained
%           INPUTS:
%                  perc - must be betweeb 0-100. This the amount of
%                  variation the retained modes should explain
            if perc<0. || perc > 100.
                error('Percentage entered is outside of range 1-100')
            end
            
            
            %index of eignevector for which variance explained exceeds perc
            nPCs =find(obj.PercVarianceExplained>perc,1,'first');
            if isempty(nPCs)
                nPCs = obj.Dim;
            end
            obj.trimNPCs(nPCs);
        end
        
        function obj =trimNPCs(obj,nPCs)
%          trimNPCs - remove modes of variation beyond a specified number           
%           INPUTS:
%               nPCs - number of modes to retain, must be an integer
%                   1-obj.Dim inclusive .
            %strip model
            assert(nPCs<=obj.Dim,'number of PCs to retain is greter than the number available')
           
            obj.EigVec = obj.EigVec(:,1:nPCs);
            obj.EigVal = obj.EigVal(1:nPCs);
            if ~isempty(obj.Projections)
                obj.Projections = obj.Projections(:,1:nPCs);
            end
        end
        
        function obj = prepare2Share(obj)
%           prepare2Share - deletes obj.Projections in preparation for
%                           wide sharing of models
            obj.Projections = [];
        end
       

    end
    methods %INTERFACING
        function out = alignIntoSpace(obj,in,w)
%           alignIntoSpace - estimate a weighted scaled rigid transformation of a shape onto the model average            
%           INPUTS: 
%                  in - an instance of a shape3D class to transform, or an
%                  k(vertices) x 3 matrix of vertex co-ordinates
%                  w - a vector of weights each corresponding to a vertex
%                     of in.
%           OUTPUTS:
%                  out - a transformed copy of in of the same type as in

            if nargin<3, w = ones(1,in.nVertices);end
            obj.AlignmentTransformation = morphableShape3D.getTransformation(obj.Average,in,true,w); 
            out = morphableShape3D.evalTransformation(obj.AlignmentTransformation,in);  
        end
        function out = projectIntoSpace(obj,in,w)
%           projectIntoSpace - performs a weighted regularized projection of a shape onto the model
%                   The behaviour of this is controlled by the values of  obj.ProjectionRegularisationType
%                   obj.ProjectionRegularisationValue        
%           INPUTS:   
%                   in - an instance of a shape3D class to transform, or an
%                   k(vertices) x 3 matrix of vertex co-ordinates
%                   w - a vector of weights each corresponding to a vertex
%                     of in.
%           OUTPUTS: 
%                   out - a vector of scores/projections onto each mode of
%                   variation.

            regtype = obj.ProjectionRegularisationType;
            regval = obj.ProjectionRegularisationValue;
            if isobject(in)
                inVEC = in.Vertices';
            else
                inVEC = in';
            end

            if nargin<3, w = ones(1,size(inVEC,2));end
            % converting input to Vector representation
            assert(size(w,2)==size(in,1),'Weights are the wrong shape') 
            inVEC = double(inVEC(:));
            inVEC = inVEC-obj.AverageVEC;
            wVEC = repmat(w,3,1);
            wVEC = double(wVEC(:));
            
            
            

            
            %Setting parameters for regularisation, 
            %irrelevant parameters for a given regularistion type are set so as to have no effect on the regularisation
            
            if strcmp(lower(regtype),'neta')
                sigma = regval*obj.EigVal(1);
            else
                sigma = 0; % will have no effect
            end     
            if strcmp(lower(regtype),'nopcs')
                EigVecs = obj.EigVec(:,1:regval);
                EigVals = obj.EigVal(1:regval);
                numPCs = regval; 
            else
                EigVecs = obj.EigVec;
                EigVals = obj.EigVal;
                numPCs = obj.Dim;
            end
            
            
            if strcmp(lower(regtype),'pvalue')
               critD =  sqrt(chi2inv(1.-regval,obj.Dim));
            elseif strcmp(lower(regtype),'mahalanobis')
               critD = regval;
            else
               critD = [];
            end    
            if strcmp(lower(regtype),'clamp')
                clampval = regval;
            else
                clampval = [];
            end
            
            
            
            % getting weight matrix and perform projection
            A = spdiags(wVEC,0,speye(length(inVEC),length(inVEC)));
            AQ = A*EigVecs*diag(EigVals);
            [U,W,V] = svd(AQ,'econ');
            W = diag(W);
            W = W./((W.^2)+ones(size(W))*sigma);
            coefs = diag(EigVals)*V*diag(W)*U'*A*inVEC;
            
            out = zeros(obj.Dim,1);
            % only has effect if numPCs < obj.Dim (i.e. if regularisation is based on the number of PCs)
            out(1:numPCs) = coefs;
            
            %scale or clamp vec if relevant, otherwise this will have no
            %effect
            out = scaleVec(obj,out,critD);
            out = clampVec(obj,out,clampval);
            
            
        end

        function out = getShapeInstance(obj,scores)
%       getShapeInstance - given a vector of scores on each mode reconstructs the shape
%             INPUTS:
%                   scores - a vector of scores shaped obj.Dim x 1;
%             OUTPUTS:
%                   out - an instance of the shape3D class representing the
%                        shape
           if isempty(obj.Average), out = []; return; end
           VEC = obj.AverageVEC + obj.EigVec*scores;
           out = clone(obj.Average);
           out.Vertices = reshape(VEC,3,obj.Average.nVertices)';
        end
        function [out,scores] = compute_morphable_transformation(obj,in,w)
%           compute_morphable_transformation - fits input shape to model.
%               aligns a shape to the model average, fits it to the model
%               with a weighted regularized projection and returns the
%               fitted shape in its original spatial location
%           INPUTS: 
%                INPUTS:   
%                   in - an instance of a shape3D class to transform, or an
%                   k(vertices) x 3 matrix of vertex co-ordinates
%                   w - a vector of weights each corresponding to a vertex
%                     of in
%                OUPUTS: an shape3D regresenting the fitted shape

            if nargin<3 
                if isobject(in)
                    nW = in.nVertices;
                else
                    nW = size(in,1);
                end
                w = ones(1,nW);
            end
            % step 1, align with average
            out = alignIntoSpace(obj,in,w);
            % step 2, project into PCA space applying regularisation
            out = projectIntoSpace(obj,out,w);
            scores = out;
           
            
            % step 4, recreate a shape
            out = getShapeInstance(obj,out);
            % step 5, transform back to original shape coordinate system
            out = morphableShape3D.evalInverseTransformation(obj.AlignmentTransformation,out);
            out.Vertices = double(out.Vertices);
        end
        
%         
        function p = batchComputeOutOfSamplePvalues(obj,shapes)
            if iscell(shapes)
                nShapes = numel(shapes);
            else
                nShapes = size(shapes,3);
            end
            p = nan(nShapes,1);
            for i = 1:nShapes
               if iscell(shapes)
                   shp = shapes{i};
               else
                   shp = shapes(:,:,i);
               end
               p(i) = obj.computeOutOfSamplePValue(shp);
                
            end
            
            
        end
        
        function [p,MD,pr,residuals] = computeOutOfSamplePValue(obj,shape)
           if iscell(shape)
               shape = shape{1};
           end
           if ismatrix(shape)
              shapeObj = shape3D;
              shapeObj.Vertices = shape;
              
           else
               shapeObj = shape;
           end
           shapeObj = obj.alignIntoSpace(shapeObj);
           v0 = arrayStructure2Vector(shapeObj.Vertices-obj.Average.Vertices);
           pr = v0*obj.EigVec;
           prZ = pr./obj.EigStd;
           MDsqu = sum(prZ.^2);
           p = chi2cdf(MDsqu,obj.Dim,'upper');
           MD = sqrt(MDsqu);
           
           % get target face residuals
           recObj = obj.getShapeInstance(pr');
           recObj = obj.evalInverseTransformation(obj.AlignmentTransformation,recObj);
           shapeObj = obj.evalInverseTransformation(obj.AlignmentTransformation,shapeObj);
           residuals = shapeObj.Vertices-recObj.Vertices;

        end
        
        end
    
methods % for visualisation
            function movie = runPC(obj,pc,range,scan,f)
             movie = struct('cdata',[],'colormap',[]);
             countFrame = 1;
             for i=range(1):0.5:range(2)
                 coeff = zeros(1,obj.Dim);
                 coeff(pc) = i*obj.EigStd(pc);
                 updateShowPC(obj,scan,coeff);
                 drawnow;pause(0.01);
                 movie(countFrame) = getframe(f.Figure);
                 countFrame = countFrame+1;
             end
             for i=range(2):-0.5:range(1)
                 coeff = zeros(1,obj.Dim);
                 coeff(pc) = i*obj.EigStd(pc);
                 updateShowPC(obj,scan,coeff);
                 drawnow;pause(0.01);
                 movie(countFrame) = getframe(f.Figure);
                 countFrame = countFrame+1;
             end
            end
          
        
        
        function movie = animatePC(obj,pc,range, viewMode)
            % a
            % animate principal component, obj, pc = principal component
            % range [-x x], 
            if nargin == 1
               pc = 1;range = [-3 3]; viewMode = 'SolidWireframe';
            elseif nargin == 2
               range = [-3 3];; viewMode = 'SolidWireframe';
                      
            end
            coeff = zeros(1,obj.Dim);
            coeff(pc) = range(1)*obj.EigStd(pc);
            [f,scan] = initializeShowPC(obj,coeff',viewMode);
            in = input('Ready? y/n:','s');
            while strcmp(in,'y')
                  movie = runPC(obj,pc,range,scan,f);
                  in = input('Again? y/n:','s');
            end


            try delete(f);catch, end %#ok<CTCH>
        end
         function [f,scan] = initializeShowPC(obj,Scores,ViewMode)
                 f = viewer3D;
                 scan = obj.getShapeInstance(Scores);

                 scan.RenderAxes = f.RenderAxes;scan.Visible = true;scan.Selected = true;
                 scan.SingleColor = [0.8 0.8 0.8];scan.ColorMode = 'Single';
                 scan.ViewMode = ViewMode;
                 scan.Material = 'Dull';
                 %disp('initializing');
                 f.SceneLightVisible = true;f.SceneLightLinked = true;
                 % transfer any landmarks from average onto image
%                  if isfield(obj.Average.UserData,'LandmarkSelection') && ~isempty(obj.Average.UserData.LandmarkSelection)
%                    
%                      [bar,index] = cart2baryKNN(obj.Average.Vertices,obj.Average.UserData.LandmarkSelection.Vertices);
%                      LMs = clone(obj.Average.UserData.LandmarkSelection);
%                      newVertices = bary2cartKNN(scan.Vertices, index, bar);
%                      LMs.Vertices = newVertices;
%                      LMs.Visible = true;
%                  else
%                      LMs = [];
%                  end
%                  if ~isempty(LMs)
%                     viewer(LMs,f)
%                  end
%                  
                 set(f.Toolbar.light_toggle,'State','on');
                 set(f.Toolbar.link_toggle,'State','on');
         end
      function updateShowPC(obj,scan, Scores)
%                     if ~isempty(LMs)
%                         [bar,index] = cart2baryKNN(scan.Vertices,LMs.Vertices);
%                     end
                     shape = obj.getShapeInstance(Scores');
                    scan.Vertices = shape.Vertices; 
%                     if ~isempty(LMs)
%                         newVertices = bary2cartKNN(scan.Vertices, index, bar);
%                         LMs.Vertices = newVertices;
%                     end
       end            
        
end
  
   methods %OTHER
%         function obj = checkRegularisationParams(obj)
%             %Checks consistency of specified regularisation type and value
%              switch obj.RegularisationType
%                 case 'None'
%                     if ~isempty(obj.RegularisationValue)
%                         warning('RegularisationType is None so RegularisationValue Should be empty')
%                     end
%                  case 'Neta'
%                      % No checks, in theory Neta can be any value
%                      %TODO Warning if it exceeds usual values?
%                  case 'NoPCs'
%                      % check if integer between 1 and the number of PCs in
%                      % the model
%                      if floor(obj.RegularisationValue)~=obj.RegularisationValue || obj.RegularisationValue < 1 || obj.RegularisationValue > obj.Dim
%                          warning('For regularisation type NoPCs Regularisation value should be an integer between 1 and the number of dimensions of the model')
%                      end
%                  case 'Clamp'
%                      if obj.RegularisationValue <0.
%                          warning('For regularisation type clamp regularisation value should be an absolute value (i.e. greater that zero)');
%                      end
%                  case 'Pvalue'
%                      if obj.RegularisationValue > 1. || obj.RegularisationValue <0.
%                         warning('For regularisation pvalue regularisation value should be 0.<=value<=1.')
%                      
%                      end
%                      
%                      
%             end
%  
%             
%         end
        
        function out = mahaDist(obj,a,b)
            if nargin<3
                b = zeros(size(a));
            
            end
                
            z=(a-b)./obj.EigStd;
            out = sqrt(sum(z.^2));
        end
        
        function out = clampVec(obj,vec,clampval)

              out = vec;
              if ~isempty(clampval)
                  % get sign and z-socre of each element in vec
                  s = sign(vec);
                  z = abs(vec./obj.EigStd);

                  %mask of values that need clamping
                  mask = clampval<z;

                  % vector of values if clamping applied
                  clampedVec = (clampval.*obj.EigStd);

                  % insert clamped values where appropriate and keep correct
                  % sign

                  out(mask) = clampedVec(mask).*s(mask);
              end
              
        end
        
        function out = scaleVec(obj,vec,critD)
            
            out = vec;
            if ~isempty(critD)
                D = mahaDist(obj,zeros(size(vec)),vec);
                if D>critD
                    ratio = critD/D;
                    out = out*ratio;
                end
            end
        end
                
            
            
            
    end
    
    methods (Static = true)
        function T = getTransformation(Target,Floating,scale,w)
%           getTransformation - computes a weighted Procrustes transformation of a floating shape onto a target shape
%           INPUTS:
%                    Target - a shape3D object or kx3 array of vertices of the target shape
%                    Floating - a shape3D object or kx3 array of vertices of the floating shape
%                    scale - if true scaling will  be estimated, if false it will not
%                    w - a vectir of k weights each corresponding to a
%                    vertex of the target and floating shapes. Those points
%                    with high weightes will have more influence on the
%                    resulting transformation
%           OUTPUTS:
%               T - a struct containing the parameters of the transformation:
%                   'Scale' a scalar representing the scaling from floating
%                   shape to target shape; 'Rotation' - a 3x3 rotation
%                   matrix; and 'Translation' a 3 element translation
%                   vector
%               out - a transformed copy of in, of the same type
%                  
            if isobject(Target)
                q = Target.Vertices';
             else
                q = Target';
             end
             
             if isobject(Floating)
                 p = Floating.Vertices';
             else
                p = Floating';
             end
             nbpts = size(p,2);
             index = find(w);% find points not having weights == 0;
             if length(index)<nbpts
                 p = p(:,index);
                 q = q(:,index);
                 w = w(:,index);
             end
             totalw = sum(w,2);
             rw = repmat(w,3,1);
             nbpts = size(p,2);
             % align centers
             centerq = sum(rw.*q,2)./repmat(totalw,3,1);
             centerp = sum(rw.*p,2)./repmat(totalw,3,1);
             q = q-repmat(centerq,1,nbpts);
             p = p-repmat(centerp,1,nbpts);
             % Scale both by making mean length 1
             if scale
                 lengths = sqrt(sum((rw.*p).^2,1));
                 meanLength1 = sum(lengths,2)/totalw;
                 p = p./meanLength1;
                 lengths = sqrt(sum((rw.*q).^2,1));
                 meanLength2 = sum(lengths,2)/totalw;
                 q = q./meanLength2;
                 scalefactor = meanLength2/meanLength1;
             else
                 scalefactor = 1;% keep scale fixed
             end
             % Rotate
             [U,S,V] = svd(rw.*q*p');
             H = V*sign(S)*U';
             H = H';
             % putting it all together
             T.Scale=scalefactor;
             T.Rotation = H;
             transform = T.Scale*H;
             Ta=eye(4);
             Ta(1:3,4)=centerq;
             Tb=eye(4);
             Tb(1:3,4)=-centerp;
             R=eye(4);
             R(1:3,1:3)=transform;
             Tout=Ta*R*Tb;
             T.Translation = Tout(1:3,4)/T.Scale;
        end
        function out = evalTransformation(T,in)
%           evalTransform - applies a scaled rigid transformation to a shape
%               INPUTS:
%                 T - a struct containing the parameters of the transformation:
%                   'Scale' a scalar representing the scaling from floating
%                   shape to target shape; 'Rotation' - a 3x3 rotation
%                   matrix; and 'Translation' a 3 element translation
%                   vector.
%                 in - the shape to transform, either a shape3D or a kx3
%                   matrix of co-ordinates.
%               OUTPUTS:
%                 out - a transformed copy of in, of the same type
%                   

             if isobject(in)
                 v = in.Vertices;
             else
                 v = in;
             end
             outv = (T.Scale*(T.Rotation*(v')+repmat(T.Translation,1,size(v,1))))';
             
             if isobject(in)
                out = clone(in);
                out.Vertices = outv;
             else
                 out = outv;
             end
        end
        function out = evalInverseTransformation(T,in)
%       evalInverseTransform - applies the inverse of the specified scaled rigid transformation to a shape%               INPUTS:
%                 T - a struct containing the parameters of the transformation:
%                   'Scale' a scalar representing the scaling from floating
%                   shape to target shape; 'Rotation' - a 3x3 rotation
%                   matrix; and 'Translation' a 3 element translation
%                   vector.
%                 in - the shape to transform, either a shape3D or a kx3
%                   matrix of co-ordinates.
           out = clone(in);
             invRotation = T.Rotation^-1;
             out.Vertices = (invRotation*(out.Vertices'/T.Scale-repmat(T.Translation,1,out.nVertices)))';
        end
        
  
        
    end
end
