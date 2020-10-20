function [out] = removeInvisibleMacFiles(in)

    out = in(~contains({in.name},'._'));


end

