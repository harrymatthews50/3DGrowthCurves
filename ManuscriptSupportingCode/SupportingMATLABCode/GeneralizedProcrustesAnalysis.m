function [AlignedLandmarks,AvgLandmarks,CentroidSizes,v] = GeneralizedProcrustesAnalysis(Landmarks,TemplateLandmarks,iter,scale,reflect,display)
  %%% implements a generalized Procrsutes analysis to co-align all shapes to the smaple mean shape
  % INPUTS: Landmarks - an n (vertices) x 3 (x,y,z) co-ordinates a x k
  %             shapes array of landmark co-ordinates
  %         TemplateLandmarks - the landmarks with which to initialise the
  %             algorithm - must either be an n x 3 matrix or a shape3D 
  %         iter - number of iterations to perform (default = 3)
  %         scale - if true (default) all faces will be scaled so that the
  %             root-mean-square distance of all points from the center is equal
  %             to one
  %         reflect - if true Landmark configurations will be allowed to reflect, if
  %              false (default) they won't
  %         display - if true the progresssion of the algorithm will be
  %         visible in a figure window (default = false)
  % OUTPUTS:
  %       AlignedLandmarks - an n (vertices) x 3 (x,y,z) co-ordinates a x k
  %             shapes array of landmark co-ordinates aligned (rotated, translated and optionally scaled) to the sample
  %             mean.
  %       AvgLandmarks - the sample mean configuration returned as either a
  %         matrix or shape3D (in the same format as TemplateLandmarks)
  %         
  %%%copyright Peter Claes and Harold Matthews (harry.matthews@kuleuven.be; peter.claes@kuleuven.be) (2020)
         if nargin<6, display = false;end
         if nargin<5, reflect = false;end
         if nargin<4, scale = true;end
         if nargin<3, iter = 3; end
               N = size(Landmarks,3);
            nLM = size(Landmarks,1);

         if ~isempty(TemplateLandmarks)
             if isobject(TemplateLandmarks)
                AvgLandmarks = TemplateLandmarks.Vertices;
             else
                 AvgLandmarks = TemplateLandmarks;
             end
         else
            AvgLandmarks = Landmarks(:,:,randi(size(Landmarks,3)));
         end
         CentroidSizes = zeros(1,N);
         AlignedLandmarks = Landmarks;
         for i=1:iter
             if display, disp([num2str(i) ' out of ' num2str(iter)]);end
             parfor n=1:1:N
                 tmp = squeeze(AlignedLandmarks(:,:,n));
                 if i==1
                    avg = mean(tmp,1);
                    Differences = repmat(avg,nLM,1)-tmp;
                    Distances = sqrt(sum(Differences.^2,2));
                    CentroidSizes(n) = sqrt(sum(Distances.^2)/nLM);      
                 end
                 [~,~,transform] = procrustes(AvgLandmarks,tmp,'Scaling',scale,'Reflection',reflect);
                 AlignedLandmarks(:,:,n) = transform.b*tmp*transform.T + repmat(transform.c(1,:),nLM,1);
             end
             AvgLandmarks = mean(AlignedLandmarks,3);
             if scale
                tmp = shape3D;tmp.Vertices = AvgLandmarks;
                AvgLandmarks = AvgLandmarks./tmp.CentroidSize;
             end
         end
         v = [];
         if display
            shape = shape3D;
            shape.Vertices = AvgLandmarks;
            shape.VertexSize = 20;
            shape.SingleColor = [0 1 0];
            v = viewer(shape);
            for n=1:1:N
               shape = shape3D;
               shape.Vertices = squeeze(AlignedLandmarks(:,:,n));
               shape.VertexSize = 10;
               viewer(shape,v);
            end
            drawnow;
         end
         if ~isempty(TemplateLandmarks) && isobject(TemplateLandmarks) % return as shape 3D
             tmp = AvgLandmarks;
             AvgLandmarks = clone(TemplateLandmarks);
             AvgLandmarks.Vertices = tmp;
         end
end