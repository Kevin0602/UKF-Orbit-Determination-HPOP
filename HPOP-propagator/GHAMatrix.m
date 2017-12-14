%--------------------------------------------------------------------------
%
% GHAMatrix: Transformation from true equator and equinox to Earth equator
%            and Greenwich meridian system 
%
% Input:
%   Mjd_UT1   Modified Julian Date UT1
% 
% Output:
%   GHAmat    Greenwich Hour Angle matrix
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
function GHAmat = GHAMatrix(Mjd_UT1)

GHAmat = R_z(gast(Mjd_UT1));

