classdef RobustRigidTransform < superHandleClass
% RobustRigidTransform - computes a robust Procrustes transformation of a floating shape onto a target shape
%               The algorithm interatively estimates a weighted Procrustes 
%               transformation and weights that signify if a given point 
%               of the shape is an inlier/outlier. the process is repeated until   
%               convergence or a maximum number of iterations is reached.
%               The Floating shape object is transformed in situ to align
%               it with the target.
% REFERENCES: P. Claes, K. Daniels, M. Walters, J. Clement, D. Vandermuelen, 
%            and P. Suetens, ?Dysmorphometrics: The modelling of morphological 
%             abnormality,? Theor. Biol. Med. Model., vol. 9, 2012, doi: 10.1186/1742-4682-9-5
% 
% RobustRigidTransform Properties
%    %% User-defined settings   
%       FloatingMesh - an instance of the shape3D class (from the
%                'meshmonk' toolbox), or an nx3 matrix of point co-ordinates that represents the floating shape        
%       TargetMesh - an instance of the shape3D class or an nx3 matrix of 
%                   point co-ordinates that represents the target shape
%        MaxIterations - the maximum number of iterations
%        FlagNormal - for development only
%        FlagAbnormal - for development only
%        Display - either true or false (default = true) if true a viewer
%               will open where the progression of the algorithm is
%               displayed. this only works if TargetMesh and FloatingMesh
%               are both shape3D objects.
%        Robust - either true or false (default=true). If true the robust
%               Procrustes transformation will be estimated iteratively, if false
%               the regular Procrustes transformation will be estimated in one
%               iteration
%        Kappa - sets a 'soft' cut-off (in standard deviation units) within the distribution of distances
%               between the floating shape and the target shape used to determine
%               the weights. Points lower than kappa are given relatively high
%               weights, points greater than kappa are given relatively low weights.
%               default = 2.
%        Scale - either true or false - if true scaling from the floating
%               shape will be scaled to the size of the target, if false the
%               floating shape will only be rotated and translated
%             
%        ModelRMSDevs - for development only
%        UseModelExpectedDevs - for development only
%   %% algorithm output
%       % RigidTransformation - a struct with fields:  Scale-the scaling from
%               the floating shape to the target shape; Rotation - a 3x3 
%               rotation matrix of the rotation from floating shape to target
%               shape; Translation - a 3 element vector of the translation 
%               from the floating shape to the target.
%       % Weights - Outlier weights
%                                                   
% Written by Harold Matthews harry.matthews@kuleuven.be Copyright 2020

properties % user settings
        FloatingMesh=[];
        TargetMesh = [];
        MaxIterations = [];
        FlagNormal = [];
        FlagAbnormal = [];
        Display = false;
        Robust = true;
        Kappa = 2;
        Scale = true;
        ModelRMSDevs;
        UseModelExpectedDevs = false;
end
properties  % algorithm output
    RigidTransformation = []
    Weights = []; 
end
    
    
    properties (Dependent = true)
        
        nVertices;
    end
    
    
    properties (Hidden=true)
        Iteration = [];
        InitialFloatingMesh = [];
        WeightViewMesh = [];
        WeightsLastIteration = [];

    end
    
     methods
        function obj = RobustRigidTransform(varargin)
             obj@superHandleClass(varargin);
        end
    end
    
    
    methods %getters
        function out = get.nVertices(obj)
            if ~isempty(obj.TargetMesh)
                if isobject(obj.TargetMesh)
                    out = obj.TargetMesh.nVertices;
                else
                    out = size(obj.TargetMesh,1);
                end
            else
                out = [];
            end
        end
        
         
    end
    methods
        function fit(obj)
         % fit - execute algorithm
           obj.Iteration = 0;
           obj.initialize(); % initialize algorithm performing non-robust Procrustes
          
           if obj.Robust
               while true % repeat until convergence or maximum iterations exceeded
                   obj.Iteration = obj.Iteration+1;
                   obj.updateWeights(); % re-estimate weights
                   obj.updateRigidTransform(); % re-estimate transform 
                   obj.updateFloatingMesh(); % apply transform to floating mesh
                   if obj.Display 
                       pause(1);
                   end
                   
                   % check if at maximum number of iterations
                   if ~isempty(obj.MaxIterations)
                      if obj.Iteration==obj.NumIterations
                          break
                      end
                   else % check convergence
                       if ~isempty(obj.WeightsLastIteration)
                           maxdiff = max(abs(obj.Weights-obj.WeightsLastIteration));
                           if maxdiff<.0001
                               break
                           end
                       end

                   end

                   obj.WeightsLastIteration=obj.Weights;
               end 
           end
        end
        
        function initialize(obj)
           % initialize - initialize algorithm 
           
           % initialize weights with ones
            
           obj.Weights = ones(obj.nVertices,1);
            
           % adjust weights for user-specified  inliers/outliers - for
           % development only
           if ~isempty(obj.FlagNormal)
               obj.Weights(logical(obj.FlagNormal))=1;
               
           end
           if ~isempty(obj.FlagAbnormal)
               obj.Weights(logical(obj.FlagAbnormal))=0;
               
           end
           
           % record initial location of the floating mesh prior to any
           % transformation
           if isobject(obj.InitialFloatingMesh)
                obj.InitialFloatingMesh = clone(obj.FloatingMesh);
           else
                obj.InitialFloatingMesh = obj.FloatingMesh;
           end
           % Initialise rigid transform
           obj.RigidTransformation = computeTransform(obj.InitialFloatingMesh,obj.TargetMesh,true,obj.Weights');
           
           % initialise display
           if obj.Display
                shapeViewer = viewer(obj.TargetMesh);
                viewer(obj.FloatingMesh, shapeViewer);
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
            % update floating mesh object
            obj.updateFloatingMesh()
             
           
           
        end
        function updateWeights(obj)
%           updateWeights - calculates outlier weights based on the distribution of distances between the template and the target            
            

            % get matrices of vertices of the floating and target meshes
            if isobject(obj.FloatingMesh)
                p = obj.FloatingMesh.Vertices;
            else
                p = obj.FloatingMesh;
            end
            if isobject(obj.TargetMesh)
                q = obj.TargetMesh.Vertices;
            else
                q = obj.TargetMesh;
            end
            
            
            % calculate distances
            diffs = p-q;
            dists = sqrt(sum(diffs.^2,2));
            
            % calculate z-scores
            if obj.UseModelExpectedDevs % alternative method for calculating z-scores (for development only)
                z=dists./obj.ModelRMSDevs;
            else
                % weighted root mean square distances
                weightedStd = sqrt(sum((dists.^2).*obj.Weights)/sum(obj.Weights));
                % convert to z-scores
                z = dists./weightedStd;
            end
            % calculate weights so as to decline steeply as z exceeds kappa
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
                          
            % update display 
            if obj.Display
                obj.WeightViewMesh.VertexValue = obj.Weights;
                drawnow();
            end
        end
        
        function updateRigidTransform(obj)            
            obj.RigidTransformation = computeTransform(obj.InitialFloatingMesh,obj.TargetMesh, obj.Scale,obj.Weights');
        end
        

        function updateFloatingMesh(obj)
           if isobject(obj.FloatingMesh)
                obj.FloatingMesh.Vertices = applyTransform(obj.InitialFloatingMesh.Vertices,obj.RigidTransformation);
           else
                obj.FloatingMesh = applyTransform(obj.InitialFloatingMesh,obj.RigidTransformation);
           end
           % update display     
           if obj.Display
                drawnow();
           end
            
        end
    end
end


