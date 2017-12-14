%--------------------------------------------------------------------------
%
% NutMatrix: Transformation from mean to true equator and equinox
%
% Input:
%   Mjd_TT    Modified Julian Date (Terrestrial Time)
%
% Output:
%   NutMat    Nutation matrix
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
function NutMat = NutMatrix(Mjd_TT)

% Mean obliquity of the ecliptic
ep = MeanObliquity(Mjd_TT);

% Nutation in longitude and obliquity
[dpsi deps] = NutAngles(Mjd_TT);

% Transformation from mean to true equator and equinox
NutMat = R_x(-ep-deps)*R_z(-dpsi)*R_x(ep);

