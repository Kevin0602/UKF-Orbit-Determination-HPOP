%--------------------------------------------------------------------------
%
% gmst: Greenwich Mean Sidereal Time
%
% Input:
%  Mjd_UT1    Modified Julian Date UT1
%
% Output:
%  gmstime	   GMST in [rad]
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
function gmstime = gmst(Mjd_UT1)

Secs = 86400;                       % Seconds per day
MJD_J2000 = 51544.5;

Mjd_0 = floor(Mjd_UT1);
UT1   = Secs*(Mjd_UT1-Mjd_0);       % [s]
T_0   = (Mjd_0  -MJD_J2000)/36525;
T     = (Mjd_UT1-MJD_J2000)/36525;

gmst  = 24110.54841 + 8640184.812866*T_0 + 1.002737909350795*UT1...
        + (0.093104-6.2e-6*T)*T*T;  % [s]

gmstime = 2*pi*Frac(gmst/Secs);     % [rad], 0..2pi

