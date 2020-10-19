function [out] = getElementsFromResidualsCell(inCell,method,value,field)
if nargin<4
    field = 'Residuals';
end

sz =size(inCell);
if numel(sz)>2
    error('bahaviour not checked for nDims >2')
end
inCell = inCell(:);
emptyInds = cellfun(@isempty,inCell);
values = getElementsFromResidualsStruct([inCell{~emptyInds}]',method,value,field);

% set up output array
szValues = size(values);
dim1dim2 = szValues(1:2);
nObs = sz(1);
nReps = sz(2);
out = nan([dim1dim2,nObs,nReps]);
notEmptyInds= find(~emptyInds);
[obsNum,repNum] = ind2sub(sz,notEmptyInds);

for i = 1:numel(notEmptyInds)
   out(:,:,obsNum(i),repNum(i)) = values(:,:,i);
end
    
end



