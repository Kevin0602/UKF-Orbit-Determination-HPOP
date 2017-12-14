%--------------------------------------------------------------------------
%
% PoleMatrix: Transformation from pseudo Earth-fixed to Earth-fixed
%             coordinates for a given date
%
% Input:
%   Pole coordinte(xp,yp)
%
% Output:
%   PoleMat   Pole matrix
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
function PoleMat =  PoleMatrix(xp,yp)

PoleMat = R_y(-xp) * R_x(-yp);

