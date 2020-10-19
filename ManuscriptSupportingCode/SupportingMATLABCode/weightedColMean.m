function [out] = weightedColMean(X,w)
% calculate the weighted column means of X
if numel(w) ~= size(X,1)
    error('the number of weights should equal the number of rows in X');
end
if size(w,1)==1
    w = w';
end

% weight x
wX = bsxfun(@times,X,w);

% sum along columns
sumX = sum(wX,1);

out = sumX./sum(w); 



end

