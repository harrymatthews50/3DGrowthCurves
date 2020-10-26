function [Z] = computeFacialSignature(model,shp)
%   computeFacialSignature - calculates facial signatures against a given
%           model of normal variation
%    INPUTS: model - a morphableShape3D object representing normal
%               variation
%           shp - a 'shape3D' of the individual for whom to calculate the facial signature
%   OUTPUTS:
%           Z - the facial signatures; this is a struct. Each field
%               contains z-scores for each point in 'shp'. fields correspond to
%               facial signatures in different dimensions ('X','Y','Z' and 'Normal')


% align patient to expected face
RRT = RobustRigidTransform;
RRT.FloatingMesh = clone(shp);
RRT.TargetMesh = clone(model.Average);
RRT.fit();
shp = RRT.FloatingMesh;
diffs = shp.Vertices-model.Average.Vertices;

% compute projections on normals
normals = model.Average.VertexNormals;
nProj = sum(diffs.*normals,2);

% calculate z scores
Z = struct;
Z.X = diffs(:,1)./model.PointStandardDevs(:,1);
Z.Y = diffs(:,2)./model.PointStandardDevs(:,2);
Z.Z = diffs(:,3)./model.PointStandardDevs(:,3);
Z.Normal = nProj./model.PointStandardDevs(:,4);







end

