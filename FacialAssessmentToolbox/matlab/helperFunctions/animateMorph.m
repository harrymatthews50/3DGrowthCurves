function [frames] = animateMorph(shp1,shp2,caricature,reps,viewerTag)
% animates morphing from shape 1 to shape 2
% INPUTS: shp1 - a shape3D of the first shape;
%         shp2 - a shape3D of the shape to morph into;
%         caricature - any value  (default=1) If caricature ~= 1 shape 2 is either 
%             exaggerated (caricature > 1) or de-exaggerated (caricature<1) by 
%             multiplying the difference between shp1 and shp 2 by the value of caricature
%             negative values will morph to an anti-face of the subject.
%         reps - number of times to repeat the morphing between the two shape
%         viewerTag - a string to identify the created viewer (and will be displayed at the to of the viewer window)
% OUTPUTS:frames -  a series of screenshots of the animation that can be saved as a movie file

if nargin<3
    caricature=1;
end

if nargin<4
    reps = 1;
end
if nargin<5
    viewerTag = ['Morph_shape1_to_shape2x',num2str(caricature)];
end

% align shape 2 to shape 1
RRT = RobustRigidTransform;
RRT.Target = shp1;
RRT.Floating = shp2;
RRT.fit();


% caricature shape 2 vertices

diff = RRT.Floating.Vertices-RRT.Target.Vertices;

shp2.Vertices = RRT.Target+diff*caricature;

obj = clone(shp1);
obj.ColorMode = 'Single';
obj.Material = 'Dull';
obj.SingleColor = [.8,.8,.8];

v = setUpViewer(obj);
v.Figure.Units = 'Pixels';
figurePos = v.Figure.Position;
v.Tag = viewerTag;
while true
    % get user input
    in = input('Are you ready? Y/N','s');
    switch in
        case 'y'
            break
        case 'n'
            return
    end
end

% determe scaling factors to iterate through
scalars = [linspace(0,1,15),linspace(1,0,15)];
scalars = repmat(scalars,[1,reps]);

diff = shp2.Vertices-shp1.Vertices;

frames = cell(1,numel(scalars));

for i = 1:numel(scalars)
    obj.Vertices=shp1.Vertices+diff*scalars(i);
    % fix pixel size
    v.Figure.Position = figurePos;
    frames{i} = getframes(v.Figure);
    pause(0.5);
end

frames = [frames{:}];


end

