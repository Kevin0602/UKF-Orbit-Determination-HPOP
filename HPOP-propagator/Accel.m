%--------------------------------------------------------------------------
%
% Accel: Computes the acceleration of an Earth orbiting satellite due to 
%    	 - Earth's harmonic gravity field (including Solid Earth Tides and
%      	   Ocean Tides), 
%    	 - gravitational perturbations of the Sun, Moon and planets
%    	 - solar radiation pressure
%    	 - atmospheric drag and
%	 	 - relativity
%
% Inputs:
%   Mjd_UTC     Modified Julian Date (UTC)
%   Y           Satellite state vector in the ICRF/EME2000 system
%   Area        Cross-section 
%   mass        Spacecraft mass
%   Cr          Radiation pressure coefficient
%   Cd          Drag coefficient
%
% Output:
%   dY		    Acceleration (a=d^2r/dt^2) in the ICRF/EME2000 system
%
% Last modified:   2015/08/12   M. Mahooti
% 
%--------------------------------------------------------------------------
function dY = Accel(x, Y)

global AuxParam eopdata

SAT_Const

[UT1_UTC,TAI_UTC,x_pole,y_pole,ddpsi,ddeps] = IERS(eopdata,AuxParam.Mjd_UTC+x/86400,'l');
[UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
Mjd_TT = AuxParam.Mjd_UTC+x/86400+TT_UTC/86400;
Mjd_UT1 = AuxParam.Mjd_UTC+x/86400+UT1_UTC/86400;

P = PrecMatrix(MJD_J2000,Mjd_TT);
N = NutMatrix(Mjd_TT);
T = N * P;
E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

[r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, ...
 r_Neptune,r_Pluto,r_Moon,r_Sun,r_SunSSB] = JPL_Eph_DE405(AuxParam.Mjd_UTC+x/86400);

% Acceleration due to harmonic gravity field
a = AccelHarmonic(AuxParam.Mjd_UTC+x/86400,r_Sun,r_Moon,Y(1:3),E);

% Luni-solar perturbations
if (AuxParam.sun)
    a = a + AccelPointMass(Y(1:3),r_Sun,GM_Sun);
end

if (AuxParam.moon)
    a = a + AccelPointMass(Y(1:3),r_Moon,GM_Moon);
end

% Planetary perturbations
if (AuxParam.planets)
    a = a + AccelPointMass(Y(1:3),r_Mercury,GM_Mercury);
    a = a + AccelPointMass(Y(1:3),r_Venus,GM_Venus);
    a = a + AccelPointMass(Y(1:3),r_Mars,GM_Mars);
    a = a + AccelPointMass(Y(1:3),r_Jupiter,GM_Jupiter);
    a = a + AccelPointMass(Y(1:3),r_Saturn,GM_Saturn);
    a = a + AccelPointMass(Y(1:3),r_Uranus,GM_Uranus);    
    a = a + AccelPointMass(Y(1:3),r_Neptune,GM_Neptune);
    a = a + AccelPointMass(Y(1:3),r_Pluto,GM_Pluto);
end

% Solar radiation pressure
if (AuxParam.sRad)
    a = a + AccelSolrad(Y(1:3),r_Earth,r_Moon,r_Sun,r_SunSSB, ...
        AuxParam.area_solar,AuxParam.mass,AuxParam.Cr,P_Sol,AU);
end

% Atmospheric drag
if (AuxParam.drag)
    % Atmospheric density
	Omega = 7292115.8553e-11+4.3e-15*( (AuxParam.Mjd_UTC+x/86400-MJD_J2000)/36525 );    
    [~,dens] = JB2008(AuxParam.Mjd_UTC+x/86400,r_Sun,T*Y(1:3));
    a = a + AccelDrag(dens,Y(1:3),Y(4:6),T,AuxParam.area_drag,AuxParam.mass,AuxParam.Cd,Omega);
end

% Relativistic Effects
if (AuxParam.Relativity)
    a = a + Relativity(Y(1:3),Y(4:6));
end

dY = [Y(4:6);a];

