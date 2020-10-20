function [out] = facialSignatureKULColormap(kappa, clim)
    if nargin<2
        clim = [-3,3];
    end
    if nargin<1
        kappa = 2;
    end
    % put clim to range 0-max
    clim = clim-clim(1);
    climRange = clim(2)-clim(1) ;
    climMid = climRange/2;
    yellowLoc = (climMid+kappa)/climRange;
    cyanLoc = (climMid-kappa)/climRange;
    out = redWhiteBlueColormap(cyanLoc,yellowLoc,0.1);


end

