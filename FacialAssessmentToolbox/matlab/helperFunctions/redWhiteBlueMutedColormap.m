function [cmap] = redWhiteBlueMutedColormap(cyanLoc,yellowLoc,transitionSpan)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
cmap = zeros(255,3);

blue = [0,0,1];
cyan = [0.55,0.83,1]; % muted version of cyan
white= [.85,.85,.85];
yellow = [1,0.74,0.25]; % muted orangey yellow
red = [1,0,0];

% get cyan location
cyanInd = floor(cyanLoc*255);
yellowInd = floor(yellowLoc*255);

transitionSpan = floor(transitionSpan*255);

cmap(1:cyanInd,:) = linInterp(blue,cyan,cyanInd);
span = (cyanInd+1):(cyanInd+transitionSpan);
cmap(span,:) = linInterp(cyan,white,numel(span));
span = (cyanInd+transitionSpan+1):(yellowInd-transitionSpan);
cmap(span,:) = repmat(white,[numel(span),1]);
span = (yellowInd-transitionSpan+1):yellowInd;
cmap(span,:) = linInterp(white, yellow,numel(span));
span = (yellowInd+1):255;
cmap(span,:) = linInterp(yellow,red,numel(span));





end

function out = linInterp(p1,p2,numPoints)
    t = linspace(0,1,numPoints)';
    out = (1-t)*p1+t*p2;
end