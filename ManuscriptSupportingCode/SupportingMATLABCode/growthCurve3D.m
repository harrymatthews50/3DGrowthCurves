classdef growthCurve3D <superHandleClass
    % this class implements calculating age and sex specific models of
    % shape variation
    % copyright Harold Matthews 2020
    properties % changes to any of these properties trigger the model to re-evaluate, unless one of them is not set
         EvaluationAge; % the age for which to produce the model of variation
         EvaluationSex; % the sex for which to produce the model of variation
         BWSigma; % the bandwidth; is either user-defined or determined dynamically if  SetBandwidthFromInterpolant = true
    end
    properties % training data always needs to be set by the user
        RefScan; % a shape3D object specifying the topology of all the training images (e.g. the template image used for registration)
        GlobalAverageShape = []; % the average shape3D of all the training shapes
        TrainingLandmarks; % shapes used for training
        Age; % a vector containing the age for each training case
        Sex; % a vector containing the sex for each training case (1=male; 2=female)
    end
    properties % optional user settings
        UserExclusionMask; % a user-defined mask to manually exclude training observations must contain true or false for each training observation (or be empty). If false, the observation will always be excluded
        MaxZ = 3; % will exclude training observations from calculating the model 
        Speedy = false; % if true will run slightly faster - for testing code only
        ScaleShapes = true; % if true shapes will be scaled to unit size, if false no scaling will occur
    end 
    properties % settings for age-adaptive setting of bandwidth
        SetBandwidthFromInterpolant = false; % determines if BWSigma will be set dynamically from the interpolant objects in 
        SetBandwidthFromInterpolantInterpolant1 = []; % the GriddedInterpolant object to interpolate the bandwidth to use given age for males
        SetBandwidthFromInterpolantInterpolant2 = [];% the GriddedInterpolant object to interpolate the bandwidth to use given age for females
    end
    properties % output
        VariationModel; % a morphable shape3D of the variatioon model for the current EvaluationAge and EvaluationSex 
        ExpectedShapeObj = []; % a shape3D opf the expected shape
        ShapeBeta; % the beta coeffiecients of the local regression model arranged per vertex (n landmarks x 3 (x,y,z) co-ordinates
    end
          
         % 

      
    properties (Dependent=true)
       ObservationMask;
       ObservationWeights;
       ZScores;
       ExpectedShape;
     
       
    end
        
    properties (Hidden=true)
        EigVecs;
        Colormap;
        Mean;
        CorrectedShapes;
        PCScores;
        ExpectedPC;
        CorrectedPC;
        PCBeta;
        ViewerObj;
        ApplyWeights = true; % for development only
        rootWeights = true; % for development only
        applyAgeCorrection = true; % for development only
        applyWeightsToModel = false; % for development only
        method = 'residualise'; % for development only
        noExtrapolation = true; % for development only
    end
    
    properties (Hidden=true, Dependent=true)
       X;
       Y;
       WeightedXMean;
       WeightedYMean;
       SexMask;
    end
    

    methods % setters
        function set.EvaluationAge(obj,in)
            obj.EvaluationAge = in;
            obj.update()
        end
        function set.EvaluationSex(obj,in)
            obj.EvaluationSex = in;
            obj.update();
        end
        function set.BWSigma(obj,in)
            obj.BWSigma = in;
            obj.update()
        end
        
        
    end
    methods % getters
        function out = get.ZScores(obj)
           % convert age into z-scores mu = obj.EvaluationAge, sigma = obj.BWSigma 
           if any([isempty(obj.Age),isempty(obj.EvaluationAge),isempty(obj.BWSigma)]); out = []; return; end
           out = (obj.Age-obj.EvaluationAge)./obj.BWSigma;
        end
        function out = get.SexMask(obj)
            if any([isempty(obj.Sex),isempty(obj.EvaluationSex)]); out = []; return; end
            out = obj.Sex==obj.EvaluationSex;
        end
        function out = get.ObservationMask(obj)
            if any([isempty(obj.SexMask),isempty(obj.ZScores),isempty(obj.MaxZ)]); out = []; return; end
            zMask = abs(obj.ZScores)<=obj.MaxZ;
            if isempty(obj.UserExclusionMask)
                out = logical(obj.SexMask.*zMask);
            else
                out = logical(obj.SexMask.*zMask.*obj.UserExclusionMask);
            end
        end
        function out = get.ObservationWeights(obj)
            if obj.ApplyWeights
                if any([isempty(obj.ZScores),isempty(obj.ObservationMask)]); out = []; return; end 
                z= obj.ZScores(obj.ObservationMask);
                out = (1/sqrt(2*pi)).*exp(-z.^2./2);
            else 
                if isempty(obj.ObservationMask); out = []; return; end 
                out = ones(sum(obj.ObservationMask),1);
            end
        end
%         function out = get.ScalarGrowthRate(obj)
%             if isempty(obj.obj.PCBeta); out = []; return; end
%             out = norm(obj.PCBeta,2);
%         end
%         
        function out = get.X(obj)
           if any([isempty(obj.Age),isempty(obj.ObservationMask)]); out = []; return; end
           out = obj.Age(obj.ObservationMask);
        end
        function out = get.Y(obj)
            if any([isempty(obj.PCScores),isempty(obj.ObservationWeights)]); out = []; return; end
            out = obj.PCScores(obj.ObservationMask,:);
        end
        function out = get.WeightedYMean(obj)
            if any([isempty(obj.Y),isempty(obj.ObservationWeights)]); out = []; return; end
            out = weightedColMean(obj.Y,obj.ObservationWeights);
        end
        function out = get.WeightedXMean(obj)
            if any([isempty(obj.X),isempty(obj.ObservationWeights)]); out = []; return; end
            out = weightedColMean(obj.X,obj.ObservationWeights);
        end
        function out = get.ExpectedShape(obj)
            if isempty(obj.ExpectedPC); out=[]; return; end
            out = obj.reconstructShapeFromPC(obj.ExpectedPC);
        end
            
        function out = get.BWSigma(obj)
            if ~obj.SetBandwidthFromInterpolant
                out = obj.BWSigma;
            else
                out = obj.getBandwidthFromInterpolant(obj.EvaluationAge,obj.EvaluationSex);

            end
          
         end
                    
                    
                    
                
            
            
            
        
                    
        
        
        
        
    end
         
    
    methods % constructor
        function obj = growthCurve3D(varargin)
            obj@superHandleClass(varargin{:});               
        end
        
    end
    

    
    methods %PCA conversion
        function fitPCA(obj,landmarks) % fit PCA (keping all components)
            vecs = arrayStructure2Vector(landmarks);
            [obj.EigVecs,obj.PCScores,~,~,~,obj.Mean] = pca(vecs);
        
        end
        function out = reconstructShapeFromPC(obj, Scores)
             out = arrayVector2Structure(Scores*obj.EigVecs'+obj.Mean);
        end
        function out = reconstructVectorFromPCs(obj, Scores)
            out = arrayVector2Structure(Scores*obj.EigVecs');
        end
        function out = convertShapeToPCs(obj,shape)
            out = (arrayStructure2Vector(shape)-obj.Mean)*obj.EigVecs;
        end
        function out = convertVectorToPCs(obj,vector);
            if numel(size(vector))>2
                vector = arrayStructure2Vector(vector);
            end
            out = vector*obj.EigVecs;
        end
    end
    methods % pipeline steps
        function initialize(obj)
           [obj.TrainingLandmarks,obj.ExpectedShapeObj] = GeneralizedProcrustesAnalysis(obj.TrainingLandmarks,obj.RefScan,3,obj.ScaleShapes);
           obj.GlobalAverageShape = clone(obj.ExpectedShapeObj);
           obj.fitPCA(obj.TrainingLandmarks);
           
        end
        function update(obj,src,evt)
            if any([isempty(obj.EvaluationAge),isempty(obj.EvaluationSex),isempty(obj.BWSigma)]); return; end
            obj.computeExpectedShape();
            obj.ageCorrectShapes();
            obj.computeVariationModel();
            %obj.ExpectedShapeObj.VertexValue = obj.ExpectedShapeVertexValues;
        end
        
        function ageCorrectShapes(obj)
            shapes = obj.reconstructShapeFromPC(obj.Y);
            if ~obj.Speedy
                % center all shapes on the expected shape
                recenteredShapes = GeneralizedProcrustesAnalysis(shapes,obj.ExpectedShape,1,obj.ScaleShapes,false);
            else
                recenteredShapes = shapes;
            end
            [X0,Y0,beta,YResiduals] = obj.localRegression(obj.X,arrayStructure2Vector(recenteredShapes),obj.EvaluationAge,arrayStructure2Vector(obj.ExpectedShape),obj.ObservationWeights);

            if obj.applyAgeCorrection
                % produce age corrected y values
                predValues = X0*beta;
                switch obj.method
                    case 'ageCorrect'
                        obj.CorrectedShapes = arrayVector2Structure((Y0-predValues)+arrayStructure2Vector(obj.ExpectedShape));
                    case 'residualise'
                        obj.CorrectedShapes = arrayVector2Structure(YResiduals);
                    otherwise
                        error('Invalid method');
                end
                        else


                obj.CorrectedShapes = recenteredShapes;
            end
            obj.ShapeBeta = arrayVector2Structure(beta);
            obj.PCBeta = obj.convertVectorToPCs(beta);
            obj.CorrectedPC = obj.convertShapeToPCs(obj.CorrectedShapes);
                
        
        end

        function computeExpectedShape(obj)
            [~,~,beta] = obj.localRegression(obj.X,obj.Y,obj.WeightedXMean,obj.WeightedYMean,obj.ObservationWeights);
            obj.ExpectedPC = (obj.EvaluationAge-obj.WeightedXMean)*beta+obj.WeightedYMean;
            if isempty(obj.ExpectedShapeObj) || ~isvalid(obj.ExpectedShapeObj)
                obj.ExpectedShapeObj = clone(obj.RefScan);
            end
            obj.ExpectedShapeObj.Vertices = obj.ExpectedShape;
%             
         end

        function computeVariationModel(obj)

            VM = morphableShape3D();
            VM.ModelTrimmingType = 'PercentageVariance';
            VM.ModelTrimmingArgs = {98.};
            if obj.applyWeightsToModel
                if strcmp(obj.method,'residualise')
                    error('Cannot apply weights while defining model if method=residualise - observations are already weighted')
                else
                    VM.fitMorphableModel(obj.CorrectedShapes,obj.RefScan.Faces,true,clone(obj.ExpectedShapeObj),obj.ObservationWeights);
                end
            else
                switch obj.method
                    case 'ageCorrect'
                        VM.fitMorphableModel(obj.CorrectedShapes,obj.RefScan.Faces,true,clone(obj.ExpectedShapeObj));
                    case 'residualise'
                        VM.fitMorphableModel(obj.CorrectedShapes,obj.RefScan.Faces,false,clone(obj.ExpectedShapeObj),[],sum(obj.ObservationWeights));
                    otherwise
                        error('invalid method')
                end
            end
            obj.VariationModel =  VM;
        end
        function computeGlobalVariationModel(obj)
           % ignore all settings and compute a variation model from all the training data
           % escept those excluded by the user's obj.UserExclusionMask
           mask = obj.UserExclusionMask;
           if isempty(mask)
               mask = true(size(obj.Age));
           end
           shp = obj.reconstructShapeFromPC(obj.PCScores(mask,:));
           VM = morphableShape3D();
           VM.ModelTrimmingType = 'PercentageVariance';
           VM.ModelTrimmingArgs = {98.};
           VM.fitMorphableModel(shp,obj.RefScan.Faces)
           obj.VariationModel = VM;
           obj.ExpectedPC = mean(obj.PCScores(mask,:));
           if isempty(obj.ExpectedShapeObj) || ~isvalid(obj.ExpectedShapeObj)
                obj.ExpectedShapeObj = clone(obj.RefScan);
            end
            obj.ExpectedShapeObj.Vertices = obj.ExpectedShape;
    
            
            
        end
        function clearToBasicSettings(obj)
            toClear = {'VariationModel','ShapeBeta','CorrectedShapes','ExpectedShapeObj','EvaluationAge','EvaluationSex'};
            for i = 1:numel(toClear)
               obj.(toClear{i}) = []; 
            end
            
        end
        
        function [X0,Y0,beta,Yresiduals] = localRegression(obj,X,Y,centerX,centerY,weights)
            X0 = X-centerX;
            Y0 = Y-centerY;

            % apply weights
            if obj.rootWeights
                w = sqrt(weights);
            else
                w = weights;
            end
            wX0 = bsxfun(@times,X0,w);
            wY0 = bsxfun(@times,Y0,w);

            % do regression
            [~,~,~,~,beta,~,~,stats] = plsregress_noCentering(wX0,wY0);
            Yresiduals = stats.Yresiduals;
        end

        
        
    end
    
    methods
         function out = getBandwidthFromInterpolant(obj,age,sex)
                if ~isempty(sex)
                   if obj.EvaluationSex==1
                       interp = obj.SetBandwidthFromInterpolantInterpolant1;
                   elseif obj.EvaluationSex==2
                       interp = obj.SetBandwidthFromInterpolantInterpolant2;
                   end
                else
                    out = [];
                    return
                end
                if isempty(interp)
                    out = [];
                    return
                end
                if ~isempty(age)
                    if obj.noExtrapolation
                       minx = min(interp.GridVectors{:});
                       maxx = max(interp.GridVectors{:});
                      
                       if age<minx
                           error('Evaluation Age less than range modelled with interpolant')
                       elseif age>maxx
                           error('Evaluation Age greater than range modelled with interpolant')
                       else
                           out = interp(obj.EvaluationAge);%;obj.SetBandwidthFromRBFPrecision],rbf); 
                       end
                    else
                         out = interp(obj.EvaluationAge);
                    end
                else
                    out = [];
                end
        
        end
    end
  
methods (Static=true)
%     function out = computePValues(shapes,model)
%         shapes = arrayStructure2Vector(shapes);
%         shapes = shapes-model.AverageVEC';
%         scores = shapes*model.EigVec;
%         normScores = scores./model.EigStd';
%         MD = sqrt(sum(normScores.^2,2));
%         out = 1-chi2cdf(MD.^2,model.Dim);
%         
%     end
% end

    
    
%     methods % for visualisation
%         function plotExpectedShape(obj)
%             if isempty(obj.ViewerObj) || ~isvalid(obj.ViewerObj)
%                 obj.ViewerObj = setUpViewer(obj.ExpectedShapeObj);
%             end
%             obj.updateColorMap()
%         end
%         function updateColorMap(obj)
%             if any([isempty(obj.ColorModeClim),~obj.FreezeClim]) && ~isempty(obj.ExpectedShapeObj.VertexValue)
%                     [directed,obj.ColorModeClim] = getAdaptiveColorScaleLimits(obj.ExpectedShapeObj.VertexValue);
%                     if directed
%                         obj.Colormap = seismicColormap();
%                     else
%                         obj.Colormap = yellowsColorMap();
%                     end
%                     
%                 
%             end
%            if ~isempty(obj.ViewerObj) && isvalid(obj.ViewerObj)
%             caxis(obj.ViewerObj.Figure.CurrentAxes,obj.ColorModeClim);
%             colormap(obj.ViewerObj.Figure.CurrentAxes,obj.Colormap);
%            end
%             
%         
%     end
    end
end


