classdef NormalEquivalent < superHandleClass
% NormalEquivalent - This class creates a 'Normal Equivalent' of a given shape via 
%               a robust fitting to a morphable shape model. The algorithm 
%               iteratively 1. calculates weights, that indicate which 
%               regions of the target face are poorly approximated by the 
%               fitting, performs a weighted and regularized fitting 
%               (projection) onto the model. The algorithm repeats until convergence or set number of iterations is
%               reached.
% 
%  
%REFERENCES:  Claes P, Daniels K, Walters M, Clement J, Vandermeulen D, Suetens P.
%               Dysmorphometrics: the modelling of morphological abnormalities. 
%               Theor Biol Med Model. 2012;9:5. Published 2012 Feb 6. doi:10.1186/1742-4682-9-5
% 
%   NormalEquivalent Properties
%   %% User Defined Settings
%      MorphableModel - an instance of the morphableShape3D, the normal 
%               equivalent is constructed by projection onto this model       
%      RegularizationAmount - the amount of regularization to apply to the
%               projection - see help for morphableShape3D.ProjectionRegularizationAmount 
%               for more information (default=.05; this default is only 
%               sensible in conjunction with 'pvalue' RegularizationMethod)            
%      RegularizationMethod - the method of regularization to apply to the projection see help for 
%               morphableShape3D.ProjectionRegularizationType for more information
%               (default = 'pvalue')
%      TargetMesh - An instance of the shape3D class (from the 'meshmonk' 
%               toolbox). This represents the shape to be fitted to the model
%      MaxIterations - An integer that sets the maximum number of iterations of the
%                      algorithm (default = 40)
%      FlagNormal - for development only
%      FlagAbnormal - for development only
%      Display - should be either true or false. If Display=true a display
%                opens in which you can view the progression of the algorithms over multiple iterations 
%                default = true
%      Kappa - sets a 'soft' cut-off (in standard deviation units) within the distribution of distances
%               between the normal equivalent and the target mesh used to determine
%               the weights. Points lower than kappa are given relatively high
%               weights, points greater than kappa are given relatively low weights.
%               default = 2.
%   %% Algorithm outputs
%        NormalEquivalentMesh - an instance of the shape3D class
%               representing the normal equivalent shape
%        Weights - vector of weights corresponding to each point on the
%               target mesh
%        Iteration - the current iteration of the algorithm, after the
%               algorithm has executed this indicates the total number
%               of iterations used
% Written by Harold Matthews harry.matthews@kuleuven.be Copyright 2020

    properties % user-defined settings
        MorphableModel=[];
        RegularizationAmount = 0.05;
        RegularizationMethod = 'pvalue';
        TargetMesh = [];
        MaxIterations = 40;
        FlagNormal = [];
        FlagAbnormal = [];
        Display = false;
        Kappa = 2;
    end
    properties  % algorithm outputs
        NormalEquivalentMesh = []; 
        Weights = [];
        Iteration = [];

    end
    
    
    properties (Dependent = true)
        MahalanobisDistance;
        MahalanobisPvalue;
        nVertices;
    end
    
    
    properties (Hidden=true)
        
        RigidTransformation = [];
        PCTransformation = [];
        WeightViewMesh = [];
        WeightsLastIteration = [];
        MeshAdjacencyMatrix = [];
    end
    
    methods
        function obj = NormalEquivalent(obj, varargin)
            % constructor
            obj@superHandleClass(varargin);
        end
    end
    
    
    methods %getters
        function out = get.nVertices(obj)
            if ~isempty(obj.TargetMesh)
                out = obj.TargetMesh.nVertices;
            else
                out = [];
            end
        end
        
        function out = get.MahalanobisDistance(obj)
            if ~isempty(obj.PCTransformation)
               out = sqrt(sum((obj.PCTransformation./obj.MorphableModel.EigStd).^2)); 
            else
               out = [];
            end
            
        end
        
        function out = get.MahalanobisPvalue(obj)
            if ~isempty(obj.MahalanobisDistance)
                out = 1-chi2cdf(obj.MahalanobisDistance^2,obj.MorphableModel.Dim);
            else
                out = [];
            end
        end
        
        
    end
    methods 
        function fit(obj)
            % fit - execute Normal Equivalent Algorithm
           obj.Iteration = 0;
           obj.initialize();
           
           while true % repeat until convergence
               obj.Iteration = obj.Iteration+1;
               obj.updateWeights();
               obj.updateRigidTransform();
               obj.updatePCTransform();
               obj.updateNormalEquivalentMesh();
               if obj.Display
                   pause(1);
               end
               
               %%%%%% Assess convergence or if numer of iterations exceeded
               if ~isempty(obj.MaxIterations)
                  if obj.Iteration==obj.MaxIterations
                      break
                  end
               end
               if ~isempty(obj.WeightsLastIteration)
                   maxdiff = max(abs(obj.Weights-obj.WeightsLastIteration));
                   if maxdiff<.0001
                       break
                   end
               end
               obj.WeightsLastIteration=obj.Weights;
           end 
        end
        
        function initialize(obj)
            % initialize - initializes the algorithm
            % initializes weights, normal equivalent and rigid transformation of the target onto the model
            %           
           
            
            
            

            % initialise weights with all ones
            obj.Weights = ones(obj.nVertices,1);
            % adjust weights for user-specified  inliers/outliers - for
            % development only
           if ~isempty(obj.FlagNormal)
               obj.Weights(logical(obj.FlagNormal))=1;
               
            end
            if ~isempty(obj.FlagAbnormal)
               obj.Weights(logical(obj.FlagAbnormal))=0;
               
            end
            % initialise normal equivalent with the average shape of the
            % morphable model
           obj.NormalEquivalentMesh = clone(obj.MorphableModel.Average);
           
           % Initialise rigid transformation from the target onto the
           % average of the model
           obj.RigidTransformation = computeTransform(obj.TargetMesh,obj.MorphableModel.Average,true,obj.Weights');
           
           %bring normal equivalent into co-ordinate system of target -
           % this is necessaru to calculate distances between them, and
           % weights for the next iteration
           obj.NormalEquivalentMesh.Vertices = applyInverseTransform(obj.NormalEquivalentMesh.Vertices,obj.RigidTransformation);
           if obj.Display % set up viewers for watching the algorithm iterations
                shapeViewer = viewer(obj.TargetMesh);
                viewer(obj.NormalEquivalentMesh,shapeViewer);
                shapeViewer.SceneLightVisible = true;
                shapeViewer.SceneLightLinked = true;
                %create a copy of the target to view as a weight map
                obj.WeightViewMesh = clone(obj.TargetMesh);
                obj.WeightViewMesh.ViewMode = 'Solid';
                
                obj.WeightViewMesh.VertexValue = obj.Weights;
                obj.WeightViewMesh.ColorMode = 'indexed';
                weightsViewer = viewer(obj.WeightViewMesh);
                weightsViewer.SceneLightVisible = true;
                weightsViewer.SceneLightLinked = true;
                while true
                    out = input('Are you ready y/n','s');
                    switch lower(out)
                        case 'y'
                            break
                        case 'n'
                            return

                    end
                end
           end
           
                   
           
        end
        function updateWeights(obj)
         % updateWeights - computes weights to use in the projection of the target onto the model
                          
            
            % calculate distances between corresponding points on the
            % normal equivalent and target mesh
            diffs = obj.NormalEquivalentMesh.Vertices-obj.TargetMesh.Vertices;
            dists = sqrt(sum(diffs.^2,2));
            
            % calculate weighted root mean square distance
            weightedStd = sqrt(sum((dists.^2).*obj.Weights)/sum(obj.Weights));
            
            % convert to z-score
            z = dists./weightedStd;
            
            % convert z-score into weights, weights are set to decrease
            % sharply as z-scores exceed kappa
            lambda = exp((-obj.Kappa^2)*0.5)/sqrt(2*pi);
            inlierprob = exp((-z.^2).*0.5)/sqrt(2*pi);
            obj.Weights = inlierprob./(inlierprob+lambda);
            
            % adjust weights for user-specified  inliers/outliers - for
            % development only
            if ~isempty(obj.FlagNormal)
               obj.Weights(logical(obj.FlagNormal))=1;
               
            end
             if ~isempty(obj.FlagAbnormal)
               obj.Weights(logical(obj.FlagAbnormal))=0;
               
             end            
            
            % update visualisation of weights 
            if obj.Display
                obj.WeightViewMesh.VertexValue = obj.Weights;
                drawnow();
            end
        end
        
        function updateRigidTransform(obj)
            % updateRigidTransform - calculates weighted scaled rigid Procrustes alignment of the target mesh onto the model average
            obj.RigidTransformation = computeTransform(obj.TargetMesh,obj.MorphableModel.Average, true,obj.Weights');
        end
        
        function updatePCTransform(obj)
            % updatePCTransform - calculates scores/projections of the target face onto the model 
            
            
            % align vertices to space
            inVertices = applyTransform(obj.TargetMesh.Vertices,obj.RigidTransformation);
            
            % regularization is handled in the method of morphableShape3D
            % 'projectIntoSpace' for this we need to set some attributes of
            % the morphable model
            obj.MorphableModel.ProjectionRegularisationType = obj.RegularizationMethod;
            obj.MorphableModel.ProjectionRegularisationValue = obj.RegularizationAmount;
            
            % calculate regularized, weighted projections onto the modes of
            % the morphable model
            obj.PCTransformation = obj.MorphableModel.projectIntoSpace(inVertices,obj.Weights);
        end
        
        function updateNormalEquivalentMesh(obj)
            % updateNormalEquivalentMesh - update the normal equivalent mesh by reconstructing it from the current scores
            
            % get shape3D reconstructed from the scores - representing the
            % shape of normal equivalent shape
            out = obj.MorphableModel.getShapeInstance(obj.PCTransformation);
            
            % go from co-ordinate system of model average to co-ordinate
            % system of target mesh
            out = applyInverseTransform(out,obj.RigidTransformation);
            
            % update vertices of the stored normal equivalent
            obj.NormalEquivalentMesh.Vertices = out.Vertices;
            
            % update display
            if obj.Display
                drawnow();
            end
            
    end
    end
end


