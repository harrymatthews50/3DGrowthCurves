function [out] = getElementsFromResidualsStruct(inStruct,method,value,field)
if nargin<4
    field = 'Residuals';
end
cellVals = arrayfun(@(x) x.(method).(field).(value),inStruct,'UniformOutput',false);
nObjects = size(cellVals,1);

%concatenate replicates along 4th dimension
repCat = arrayfun(@(x) cat(4,cellVals{x,:}),1:nObjects,'UniformOutput',false);
out = cat(3, repCat{:});

end

