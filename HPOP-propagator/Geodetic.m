%--------------------------------------------------------------------------
%
% Geodetic: geodetic coordinates (Longitude [rad], latitude [rad],
%           altitude [m]) from given position vector (r [m])
%
% Last modified:   2015/08/12   M. Mahooti
% 
%--------------------------------------------------------------------------
function [lon, lat, h] = Geodetic(r)

SAT_Const

R_equ = R_Earth;
f     = f_Earth;

epsRequ = eps*R_equ;      % Convergence criterion
e2      = f*(2-f);        % Square of eccentricity

X = r(1);                 % Cartesian coordinates
Y = r(2);
Z = r(3);
rho2 = X*X + Y*Y;         % Square of distance from z-axis

% Check validity of input data
if (norm(r)==0)
    disp ( ' invalid input in Geodetic constructor\n' );
    lon = 0;
    lat = 0;
    h   = -R_Earth;
    return
end

% Iteration 
dZ = e2*Z;

while(1)
    ZdZ    = Z + dZ;
    Nh     = sqrt ( rho2 + ZdZ*ZdZ ); 
    SinPhi = ZdZ / Nh;    % Sine of geodetic latitude
    N      = R_equ / sqrt(1-e2*SinPhi*SinPhi);
    dZ_new = N*e2*SinPhi;
    if (abs(dZ-dZ_new) < epsRequ)
        break
    end
    dZ = dZ_new;
end

% Longitude, latitude, altitude
lon = atan2(Y, X);
lat = atan2(ZdZ, sqrt(rho2));
h   = Nh - N;

