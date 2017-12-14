%--------------------------------------------------------------------------
%
% AccelHarmonic: Computes the acceleration due to the harmonic gravity
%                field of the central body (Solid Earth Tides and Ocean
%                Tides effects are considered)
%
% Inputs:
%   r_Sun       Geocentric equatorial position (in [m]) referred to the
%               mean equator and equinox of J2000 (EME2000, ICRF)
%   r_Moon      Geocentric equatorial position (in [m]) referred to the
%               mean equator and equinox of J2000 (EME2000, ICRF)
%   Mjd_UTC     Modified Julian Date of UTC
%   r           Satellite position vector in the inertial system
%   E           Transformation matrix to body-fixed system
%   Cnm,Snm     Spherical harmonics coefficients (normalized)
%   n_max       Maximum degree
%   m_max       Maximum order (m_max<=n_max; m_max=0 for zonals, only)
%
% Output:
%   a           Acceleration (a=d^2r/dt^2)
%
% Last modified:   2015/08/12   M. Mahooti
% 
%--------------------------------------------------------------------------
function [a] = AccelHarmonic(Mjd_UTC, r_Sun, r_Moon, r, E)

global Cnm Snm eopdata AuxParam

SAT_Const

r_ref = 6378.1363e3;   % Earth's radius [m]; EGM96
gm    = 398600.4415e9; % [m^3/s^2]; EGM96

C = Cnm;
S = Snm;

r_Moon = E*r_Moon;
[lM, phiM, rM] = CalcPolarAngles(r_Moon);
r_Sun = E*r_Sun;
[lS, phiS, rS] = CalcPolarAngles(r_Sun);

[UT1_UTC,TAI_UTC,x_pole,y_pole,ddpsi,ddeps] = IERS(eopdata,Mjd_UTC,'l');
[UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC] = timediff(UT1_UTC, TAI_UTC);
Mjd_TT = Mjd_UTC + TT_UTC/86400;
Mjd_UT1 = Mjd_UTC + UT1_UTC/86400;

T  = (Mjd_TT-MJD_J2000)/36525;
T2 = T*T;
T3 = T2*T;
rev = 360*3600;  % arcsec/revolution

if (AuxParam.SolidEarthTides)
    % Effect of Solid Earth Tides    
    coeff0 = [...
    %l  l' F  D Om   dkf(R)     Amp(ip)     dkf(I)      AMP(op)
     1, 0, 2, 0, 1, -0.00044,    -0.1,     -0.00045,    -0.1;
     1, 0, 2, 0, 2, -0.00044,    -0.7,     -0.00046,    -0.7;   % Q1
    -1, 0, 2, 2, 2, -0.00047,    -0.1,     -0.00049,    -0.1;   % rho1
     0, 0, 2, 0, 1, -0.00081,    -1.2,     -0.00082,    -1.3;     
     0, 0, 2, 0, 2, -0.00081,    -6.6,     -0.00082,    -6.7;   % O1
     1, 0, 2,-2, 2, -0.00167,     0.1,     -0.00168,     0.1;   % Ntau1
    -1, 0, 2, 0, 2, -0.00193,     0.4,     -0.00193,     0.4;   % LK1
     1, 0, 0, 0, 0, -0.00196,     1.3,     -0.00197,     1.3;   % NO1
     1, 0, 0, 0, 1, -0.00197,     0.2,     -0.00198,     0.3;
    -1, 0, 0, 2, 0, -0.00231,     0.3,     -0.00231,     0.3;   % X1
     0, 1, 2,-2, 2, -0.00834,    -1.9,     -0.00832,    -1.9;   % pi1
     0, 0, 2,-2, 1, -0.01114,     0.5,     -0.01111,     0.5;
     0, 0, 2,-2, 2, -0.01135,   -43.3,     -0.01132,   -43.2;   % P1
     0, 1, 0, 0, 0, -0.01650,     2.0,     -0.01642,     2.0;   % S1
     0, 0, 0, 0,-1, -0.03854,    -8.8,     -0.03846,    -8.8;
     0, 0, 0, 0, 0, -0.04093,   472.0,     -0.04085,   471.0;   % K1
     0, 0, 0, 0, 1, -0.04365,    68.3,     -0.04357,    68.2;
     0, 0, 0, 0, 2, -0.04678,    -1.6,     -0.04670,    -1.6;
     0,-1, 0, 0, 0,  0.23083,   -20.8,      0.22609,   -20.4;   % sai1
     0, 0,-2, 2,-2,  0.03051,    -5.0,      0.03027,    -5.0;   % phi1
     1, 0, 0,-2, 0,  0.00374,    -0.5,      0.00371,    -0.5;   % theta1
    -1, 0, 0, 0, 0,  0.00329,    -2.1,      0.00325,    -2.1;   % J1
    -1, 0, 0, 0, 1,  0.00327,    -0.4,      0.00324,    -0.4;
     0, 0, 0,-2, 0,  0.00198,    -0.2,      0.00195,    -0.2;   % SO1
     0, 0,-2, 0,-2,  0.00187,    -0.7,      0.00184,    -0.6;   % OO1
     0, 0,-2, 0,-1,  0.00187,    -0.4,      0.00184,    -0.4;
    ];
    
    coeff1 = [...
    %l  l' F  D Om   dkf(el)   Amp(elas)   dkf(anel)    AMP(anel)
     0, 0, 0, 0, 1,  0.01347,    16.6,     -0.00541,    -6.7;
     0, 0, 0, 0, 2,  0.01124,    -0.1,     -0.00488,     0.1;
     0,-1, 0, 0, 0,  0.00547,    -1.2,     -0.00349,     0.8;   % Sa
     0, 0,-2, 2,-2,  0.00403,    -5.5,     -0.00315,     4.3;   % Ssa
     0, 0,-2, 2,-1,  0.00398,     0.1,     -0.00313,    -0.1;
     0,-1,-2, 2,-2,  0.00326,    -0.3,     -0.00296,     0.2;
     1, 0, 0,-2, 0,  0.00101,    -0.3,     -0.00242,     0.7;   % Msm
    -1, 0, 0, 0,-1,  0.00080,     0.1,     -0.00237,    -0.2;
    -1, 0, 0, 0, 0,  0.00080,    -1.2,     -0.00237,     3.7;   % Mm
    -1, 0, 0, 0, 1,  0.00079,     0.1,     -0.00237,    -0.2;
     1, 0,-2, 0,-2,  0.00077,     0.1,     -0.00236,    -0.2;
     0, 0, 0,-2, 0, -0.00009,     0.0,     -0.00216,     0.6;   % Msf
    -2, 0, 0, 0, 0, -0.00018,     0.0,     -0.00213,     0.3;
     0, 0,-2, 0,-2, -0.00019,     0.6,     -0.00213,     6.3;   % Mf
     0, 0,-2, 0,-1, -0.00019,     0.2,     -0.00213,     2.6;
     0, 0,-2, 0, 0, -0.00019,     0.0,     -0.00213,     0.2;
     1, 0,-2,-2,-2, -0.00065,     0.1,     -0.00202,     0.2;   % Mstm
    -1, 0,-2, 0,-2, -0.00071,     0.4,     -0.00201,     1.1;   % Mtm
    -1, 0,-2, 0,-1, -0.00071,     0.2,     -0.00201,     0.5;
     0, 0,-2,-2,-2, -0.00102,     0.1,     -0.00193,     0.2;   % Msqm
    -2, 0,-2, 0,-2, -0.00106,     0.1,     -0.00192,     0.1;   % Mqm
    ];
    
    coeff2 = [...
    %l  l' F  D Om   dkf(R)      Amp
     1, 0, 2, 0, 2,  0.00006,   -0.3;   % N2
     0, 0, 2, 0, 2,  0.00004,   -1.2;   % M2
    ];
    
    % Mean arguments of luni-solar motion
    %
    %   l   mean anomaly of the Moon
    %   l'  mean anomaly of the Sun
    %   F   mean argument of latitude
    %   D   mean longitude elongation of the Moon from the Sun
    %   Om  mean longitude of the ascending node of the Moon
    l  = mod (  485866.733 + (1325.0*rev +  715922.633)*T...
                              + 31.310*T2 + 0.064*T3, rev );
    lp = mod ( 1287099.804 + (  99.0*rev + 1292581.224)*T...
                              -  0.577*T2 - 0.012*T3, rev );
    F  = mod (  335778.877 + (1342.0*rev +  295263.137)*T...
                              - 13.257*T2 + 0.011*T3, rev );
    D  = mod ( 1072261.307 + (1236.0*rev + 1105601.328)*T...
                              -  6.891*T2 + 0.019*T3, rev );
    Om = mod (  450160.280 - (   5.0*rev +  482890.539)*T...
                              +  7.455*T2 + 0.008*T3, rev );
    
    % STEP1 CORRECTIONS
    [lgM, dlgM] = Legendre(4,4,phiM);
    [lgS, dlgS] = Legendre(4,4,phiS);
    dCnm20 = (0.29525/5)*( (GM_Moon/gm)*((r_ref/rM)^3)*lgM(3,1)...
           + (GM_Sun/gm)*((r_ref/rS)^3)*lgS(3,1) );
    dCnm21 = (0.29470/5)*( (GM_Moon/gm)*((r_ref/rM)^3)*lgM(3,2)*cos(lM)...
           + (GM_Sun/gm)*((r_ref/rS)^3)*lgS(3,2)*cos(lS) );
    dSnm21 = (0.29470/5)*( (GM_Moon/gm)*((r_ref/rM)^3)*lgM(3,2)*(sin(lM))...
           + (GM_Sun/gm)*((r_ref/rS)^3)*lgS(3,2)*(sin(lS)) );
    dCnm22 = (0.29801/5)*( (GM_Moon/gm)*((r_ref/rM)^3)*lgM(3,3)*cos(lM)...
           + (GM_Sun/gm)*((r_ref/rS)^3)*lgS(3,3)*cos(lS) );
    dSnm22 = (0.29801/5)*( (GM_Moon/gm)*((r_ref/rM)^3)*lgM(3,3)*(sin(2*lM))...
           + (GM_Sun/gm)*((r_ref/rS)^3)*lgS(3,3)*(sin(2*lS)) );
    
    dCnm20 = dCnm20+(0.30190/5)*( (GM_Moon/gm)*((r_ref/rM)^3)*lgM(3,1)...
           + (GM_Sun/gm)*((r_ref/rS)^3)*lgS(3,1) );
    dCnm21 = dCnm21+(0.29830/5)*( (GM_Moon/gm)*((r_ref/rM)^3)*lgM(3,2)*cos(lM)...
           + (GM_Sun/gm)*((r_ref/rS)^3)*lgS(3,2)*cos(lS) );
    dSnm21 = dSnm21+(-0.00144/5)*( (GM_Moon/gm)*((r_ref/rM)^3)*lgM(3,2)*(sin(lM))...
           + (GM_Sun/gm)*((r_ref/rS)^3)*lgS(3,2)*(sin(lS)) );
    dCnm22 = dCnm22+(0.30102/5)*( (GM_Moon/gm)*((r_ref/rM)^3)*lgM(3,3)*cos(lM)...
           + (GM_Sun/gm)*((r_ref/rS)^3)*lgS(3,3)*cos(lS) );
    dSnm22 = dSnm22+(-0.00130/5)*( (GM_Moon/gm)*((r_ref/rM)^3)*lgM(3,3)*(sin(2*lM))...
           + (GM_Sun/gm)*((r_ref/rS)^3)*lgS(3,3)*(sin(2*lS)) );
    
    dCnm30 = (0.093/7.0)*( (GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,1)...
           + (GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,1) );
    dCnm31 = (0.093/7.0)*( (GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,2)*cos(lM)...
           + (GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,2)*cos(lS) );
    dSnm31 = (0.093/7.0)*( (GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,2)*sin(lM)...
           + (GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,2)*sin(lS) );
    dCnm32 = (0.093/7.0)*( (GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,3)*cos(2*lM)...
           + (GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,3)*cos(2*lS) );
    dSnm32 = (0.093/7.0)*( (GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,3)*sin(2*lM)...
           + (GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,3)*sin(2*lS) );
    dCnm33 = (0.094/7.0)*( (GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,4)*cos(3*lM)...
           + (GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,4)*cos(3*lS) );
    dSnm33 = (0.094/7.0)*( (GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,4)*sin(3*lM)...
           + (GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,4)*sin(3*lS) );
    
    dCnm40 = (-0.00087/5.0)*( (GM_Moon/gm)*((r_ref/rM)^3.0)*lgM(5,1)...
           + (GM_Sun/gm)*((r_ref/rS)^3.0)*lgS(5,1) );
    dCnm41 = (-0.00079/5.0)*( (GM_Moon/gm)*((r_ref/rM)^3.0)*lgM(5,2)*cos(lM)...
           + (GM_Sun/gm)*((r_ref/rS)^3.0)*lgS(5,2)*cos(lS) );
    dSnm41 = (-0.00079/5.0)*( (GM_Moon/gm)*((r_ref/rM)^3.0)*lgM(5,2)*sin(lM)...
           + (GM_Sun/gm)*((r_ref/rS)^3.0)*lgS(5,2)*sin(lS) );
    dCnm42 = (-0.00057/5.0)*( (GM_Moon/gm)*((r_ref/rM)^3.0)*lgM(5,3)*cos(2*lM)...
           + (GM_Sun/gm)*((r_ref/rS)^3.0)*lgS(5,3)*cos(2*lS) );
    dSnm42 = (-0.00057/5.0)*( (GM_Moon/gm)*((r_ref/rM)^3.0)*lgM(5,3)*sin(2*lM)...
           + (GM_Sun/gm)*((r_ref/rS)^3.0)*lgS(5,3)*sin(2*lS) );
    
    dCnm40 = dCnm40+(-0.00089/5.0)*( (GM_Moon/gm)*((r_ref/rM)^3.0)*lgM(5,1)...
           + (GM_Sun/gm)*((r_ref/rS)^3.0)*lgS(5,1) );
    dCnm41 = dCnm41+(-0.00080/5.0)*( (GM_Moon/gm)*((r_ref/rM)^3.0)*lgM(5,2)*cos(lM)...
           + (GM_Sun/gm)*((r_ref/rS)^3.0)*lgS(5,2)*cos(lS) );
    dSnm41 = dSnm41+(-0.00080/5.0)*( (GM_Moon/gm)*((r_ref/rM)^3.0)*lgM(5,2)*sin(lM)...
           + (GM_Sun/gm)*((r_ref/rS)^3.0)*lgS(5,2)*sin(lS) );
    dCnm42 = dCnm42+(-0.00057/5.0)*( (GM_Moon/gm)*((r_ref/rM)^3.0)*lgM(5,3)*cos(2*lM)...
           + (GM_Sun/gm)*((r_ref/rS)^3.0)*lgS(5,3)*cos(2*lS) );
    dSnm42 = dSnm42+(-0.00057/5.0)*( (GM_Moon/gm)*((r_ref/rM)^3.0)*lgM(5,3)*sin(2*lM)...
           + (GM_Sun/gm)*((r_ref/rS)^3.0)*lgS(5,3)*sin(2*lS) );
    
    % STEP2 CORRECTIONS
    A0 = 1/(r_ref*sqrt(4*pi));
    dC20 = 0;
    
    for i=1:21
        theta_f = - coeff1(i,1:5)*[l lp F D Om]';
        dC20 = dC20 + A0*1e-12*(coeff1(i,7)*coeff1(i,6)*cos(theta_f)...
                               -coeff1(i,9)*coeff1(i,8)*sin(theta_f));
        theta_f = 0;
    end
    dCnm20 = dCnm20 + dC20;
    
    theta_g = gmst(Mjd_UT1);
    A1 = -1/(r_ref*sqrt(8*pi));
    dC21 = 0;
    dS21 = 0;
     
    for i=1:26
        theta_f = theta_f + (theta_g+pi)-coeff0(i,1:5)*[l lp F D Om]';
        dC21 = dC21 + A1*1e-12*coeff0(i,7)*(-coeff0(i,6)*sin(theta_f));
        dS21 = dS21 + A1*1e-12*coeff0(i,7)*(coeff0(i,6)*cos(theta_f));
        theta_f = 0;
    end
    dCnm21 = dCnm21 + dC21;
    dSnm21 = dSnm21 + dS21;
    
    dC21 = 0;
    dS21 = 0;
    
    for i=1:26
        theta_f = theta_f + (theta_g+pi)-coeff0(i,1:5)*[l lp F D Om]';
        dC21 = dC21 + A1*1e-12*coeff0(i,9)*(-coeff0(i,8)*sin(theta_f));
        dS21 = dS21 + A1*1e-12*coeff0(i,9)*(coeff0(i,8)*cos(theta_f));
        theta_f = 0;
    end
    dCnm21 = dCnm21 + dC21;
    dSnm21 = dSnm21 + dS21;
    
    A2 = 1/(r_ref*sqrt(8*pi));
    dC22 = 0;
    dS22 = 0;
    
    for i=1:2
        theta_f = theta_f + 2*(theta_g+pi)-coeff2(i,1:5)*[l lp F D Om]';
        dC22 = dC22 + A2*1e-12*coeff2(i,7)*(coeff2(i,6)*cos(theta_f));
        dS22 = dS22 + A2*1e-12*coeff2(i,7)*(coeff2(i,6)*sin(theta_f));
        theta_f = 0.0;
    end
    dCnm22 = dCnm22 + dC22;
    dSnm22 = dSnm22 + dS22;
    
    % Treatment of the Permanent Tide
    dC20 = 1e-12*4.4228e-8*(-0.31460)*0.29525; % elastic Earth
    dCnm20 = dCnm20 - dC20;
    dC20 = 1e-12*4.4228e-8*(-0.31460)*0.30190; % anelastic Earth
    dCnm20 = dCnm20 - dC20;
    
    % Effect of Solid Earth Pole Tide (elastic Earth)
    dC21 = -1.290e-9*(x_pole);
    dS21 = 1.290e-9*(y_pole);
    dCnm21 = dCnm21 + dC21;
    dSnm21 = dSnm21 + dS21;
    
    % Effect of Solid Earth Pole Tide (anelastic Earth)
    dC21 = -1.348e-9*(x_pole+0.0112*y_pole);
    dS21 = 1.348e-9*(y_pole-0.0112*x_pole);
    dCnm21 = dCnm21 + dC21;
    dSnm21 = dSnm21 + dS21;
    
    C(3,1) = C(3,1) + dCnm20;
    C(3,2) = C(3,2) + dCnm21;
    C(3,3) = C(3,3) + dCnm22;
    S(3,2) = S(3,2) + dSnm21;
    S(3,3) = S(3,3) + dSnm22;
    
    C(4,1) = C(4,1) + dCnm30;
    C(4,2) = C(4,2) + dCnm31;
    C(4,3) = C(4,3) + dCnm32;
    C(4,4) = C(4,4) + dCnm33;
    S(4,2) = S(4,2) + dSnm31;
    S(4,3) = S(4,3) + dSnm32;
    S(4,4) = S(4,4) + dSnm33;
    
    C(5,1) = C(5,1) + dCnm40;
    C(5,2) = C(5,2) + dCnm41;
    C(5,3) = C(5,3) + dCnm42;
    S(5,2) = S(5,2) + dSnm41;
    S(5,3) = S(5,3) + dSnm42;    
end

if (AuxParam.OceanTides)
    % Ocean Tides
    [lgM, dlgM] = Legendre(6,6,phiM);
    [lgS, dlgS] = Legendre(6,6,phiS);
    
    dCnm20 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.3075)/5*( (GM_Moon/gm)*((r_ref/rM)^3)*lgM(3,1)...
           + (GM_Sun/gm)*((r_ref/rS)^3)*lgS(3,1) );
    dCnm21 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.3075)/5*( (GM_Moon/gm)*((r_ref/rM)^3)*lgM(3,2)*cos(lM)...
           + (GM_Sun/gm)*((r_ref/rS)^3)*lgS(3,2)*cos(lS) );
    dSnm21 = -0.3075/5*( (GM_Moon/gm)*((r_ref/rM)^3)*lgM(3,2)*sin(lM)...
           + (GM_Sun/gm)*((r_ref/rS)^3)*lgS(3,2)*sin(lS) );
    dCnm22 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.3075)/5*( (GM_Moon/gm)*((r_ref/rM)^3)*lgM(3,3)*cos(2*lM)...
           + (GM_Sun/gm)*((r_ref/rS)^3)*lgS(3,3)*cos(2*lS) );
    dSnm22 = -0.3075/5*( (GM_Moon/gm)*((r_ref/rM)^3)*lgM(3,3)*sin(2*lM)...
           + (GM_Sun/gm)*((r_ref/rS)^3)*lgS(3,3)*sin(2*lS) );
    dCnm30 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.195)/7*( (GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,1)...
           + (GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,1) );
    dCnm31 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.195)/7*( (GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,2)*cos(lM)...
           + (GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,2)*cos(lS) );
    dSnm31 = -0.195/7*( (GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,2)*sin(lM)...
           + (GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,2)*sin(lS) );
    dCnm32 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.195)/7*( (GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,3)*cos(2*lM)...
           + (GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,3)*cos(2*lS) );
    dSnm32 = -0.195/7*( (GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,3)*sin(2*lM)...
           + (GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,3)*sin(2*lS) );
    dCnm33 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.195)/7*( (GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,4)*cos(3*lM)...
           + (GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,4)*cos(3*lS) );
    dSnm33 = -0.195/7*( (GM_Moon/gm)*((r_ref/rM)^4)*lgM(4,4)*sin(3*lM)...
           + (GM_Sun/gm)*((r_ref/rS)^4)*lgS(4,4)*sin(3*lS) );
    dCnm40 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.132)/9*( (GM_Moon/gm)*((r_ref/rM)^5)*lgM(5,1)...
           + (GM_Sun/gm)*((r_ref/rS)^5)*lgS(5,1) );
    dCnm41 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.132)/9*( (GM_Moon/gm)*((r_ref/rM)^5)*lgM(5,2)*cos(lM)...
           + (GM_Sun/gm)*((r_ref/rS)^5)*lgS(5,2)*cos(lS) );
    dSnm41 = -0.132/9*( (GM_Moon/gm)*((r_ref/rM)^5)*lgM(5,2)*sin(lM)...
           + (GM_Sun/gm)*((r_ref/rS)^5)*lgS(5,2)*sin(lS) );
    dCnm42 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.132)/9*( (GM_Moon/gm)*((r_ref/rM)^5)*lgM(5,3)*cos(2*lM)...
           + (GM_Sun/gm)*((r_ref/rS)^5)*lgS(5,3)*cos(2*lS) );
    dSnm42 = -0.132/9*( (GM_Moon/gm)*((r_ref/rM)^5)*lgM(5,3)*sin(2*lM)...
           + (GM_Sun/gm)*((r_ref/rS)^5)*lgS(5,3)*sin(2*lS) );
    dCnm43 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.132)/9*( (GM_Moon/gm)*((r_ref/rM)^5)*lgM(5,4)*cos(3*lM)...
           + (GM_Sun/gm)*((r_ref/rS)^5)*lgS(5,4)*cos(3*lS) );
    dSnm43 = -0.132/9*( (GM_Moon/gm)*((r_ref/rM)^5)*lgM(5,4)*sin(3*lM)...
           + (GM_Sun/gm)*((r_ref/rS)^5)*lgS(5,4)*sin(3*lS) );
    dCnm44 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.132)/9*( (GM_Moon/gm)*((r_ref/rM)^5)*lgM(5,5)*cos(4*lM)...
           + (GM_Sun/gm)*((r_ref/rS)^5)*lgS(5,5)*cos(4*lS) );
    dSnm44 = -0.132/9*( (GM_Moon/gm)*((r_ref/rM)^5)*lgM(5,5)*sin(4*lM)...
           + (GM_Sun/gm)*((r_ref/rS)^5)*lgS(5,5)*sin(4*lS) );
    dCnm50 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.1032)/11*( (GM_Moon/gm)*((r_ref/rM)^6)*lgM(6,1)...
           + (GM_Sun/gm)*((r_ref/rS)^6)*lgS(6,1) );
    dCnm51 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.1032)/11*( (GM_Moon/gm)*((r_ref/rM)^6)*lgM(6,2)*cos(lM)...
           + (GM_Sun/gm)*((r_ref/rS)^6)*lgS(6,2)*cos(lS) );
    dSnm51 = -0.1032/9*( (GM_Moon/gm)*((r_ref/rM)^6)*lgM(6,2)*sin(lM)...
           + (GM_Sun/gm)*((r_ref/rS)^6)*lgS(6,2)*sin(lS) );
    dCnm52 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.1032)/11*( (GM_Moon/gm)*((r_ref/rM)^6)*lgM(6,3)*cos(2*lM)...
           + (GM_Sun/gm)*((r_ref/rS)^6)*lgS(6,3)*cos(2*lS) );
    dSnm52 = -0.1032/9*( (GM_Moon/gm)*((r_ref/rM)^6)*lgM(6,3)*sin(2*lM)...
           + (GM_Sun/gm)*((r_ref/rS)^6)*lgS(6,3)*sin(2*lS) );
    dCnm53 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.1032)/11*( (GM_Moon/gm)*((r_ref/rM)^6)*lgM(6,4)*cos(3*lM)...
           + (GM_Sun/gm)*((r_ref/rS)^6)*lgS(6,4)*cos(3*lS) );
    dSnm53 = -0.1032/9*( (GM_Moon/gm)*((r_ref/rM)^6)*lgM(6,4)*sin(3*lM)...
           + (GM_Sun/gm)*((r_ref/rS)^6)*lgS(6,4)*sin(3*lS) );
    dCnm54 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.1032)/11*( (GM_Moon/gm)*((r_ref/rM)^6)*lgM(6,5)*cos(4*lM)...
           + (GM_Sun/gm)*((r_ref/rS)^6)*lgS(6,5)*cos(4*lS) );
    dSnm54 = -0.1032/9*( (GM_Moon/gm)*((r_ref/rM)^6)*lgM(6,5)*sin(4*lM)...
           + (GM_Sun/gm)*((r_ref/rS)^6)*lgS(6,5)*sin(4*lS) );
    dCnm55 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.1032)/11*( (GM_Moon/gm)*((r_ref/rM)^6)*lgM(6,6)*cos(5*lM)...
           + (GM_Sun/gm)*((r_ref/rS)^6)*lgS(6,6)*cos(5*lS) );
    dSnm55 = -0.1032/9*( (GM_Moon/gm)*((r_ref/rM)^6)*lgM(6,6)*sin(5*lM)...
           + (GM_Sun/gm)*((r_ref/rS)^6)*lgS(6,6)*sin(5*lS) );
    dCnm60 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.0892)/13*( (GM_Moon/gm)*((r_ref/rM)^7)*lgM(7,1)...
           + (GM_Sun/gm)*((r_ref/rS)^7)*lgS(7,1) );
    dCnm61 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.0892)/13*( (GM_Moon/gm)*((r_ref/rM)^7)*lgM(7,2)*cos(lM)...
           + (GM_Sun/gm)*((r_ref/rS)^7)*lgS(7,2)*cos(lS) );
    dSnm61 = -0.0892/9*( (GM_Moon/gm)*((r_ref/rM)^7)*lgM(7,2)*sin(lM)...
           + (GM_Sun/gm)*((r_ref/rS)^7)*lgS(7,2)*sin(lS) );
    dCnm62 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.0892)/13*( (GM_Moon/gm)*((r_ref/rM)^7)*lgM(7,3)*cos(2*lM)...
           + (GM_Sun/gm)*((r_ref/rS)^7)*lgS(7,3)*cos(2*lS) );
    dSnm62 = -0.0892/9*( (GM_Moon/gm)*((r_ref/rM)^7)*lgM(7,3)*sin(2*lM)...
           + (GM_Sun/gm)*((r_ref/rS)^7)*lgS(7,3)*sin(2*lS) );
    dCnm63 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.0892)/13*( (GM_Moon/gm)*((r_ref/rM)^7)*lgM(7,4)*cos(3*lM)...
           + (GM_Sun/gm)*((r_ref/rS)^7)*lgS(7,4)*cos(3*lS) );
    dSnm63 = -0.0892/9*( (GM_Moon/gm)*((r_ref/rM)^7)*lgM(7,4)*sin(3*lM)...
           + (GM_Sun/gm)*((r_ref/rS)^7)*lgS(7,4)*sin(3*lS) );
    dCnm64 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.0892)/13*( (GM_Moon/gm)*((r_ref/rM)^7)*lgM(7,5)*cos(4*lM)...
           + (GM_Sun/gm)*((r_ref/rS)^7)*lgS(7,5)*cos(4*lS) );
    dSnm64 = -0.0892/9*( (GM_Moon/gm)*((r_ref/rM)^7)*lgM(7,5)*sin(4*lM)...
           + (GM_Sun/gm)*((r_ref/rS)^7)*lgS(7,5)*sin(4*lS) );
    dCnm65 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.0892)/13*( (GM_Moon/gm)*((r_ref/rM)^7)*lgM(7,6)*cos(5*lM)...
           + (GM_Sun/gm)*((r_ref/rS)^7)*lgS(7,6)*cos(5*lS) );
    dSnm65 = -0.0892/9*( (GM_Moon/gm)*((r_ref/rM)^7)*lgM(7,6)*sin(5*lM)...
           + (GM_Sun/gm)*((r_ref/rS)^7)*lgS(7,6)*sin(5*lS) );
    dCnm66 = 4*pi*r_ref^2*1025/(5.9722e24)*(1-0.0892)/13*( (GM_Moon/gm)*((r_ref/rM)^7)*lgM(7,7)*cos(6*lM)...
           + (GM_Sun/gm)*((r_ref/rS)^7)*lgS(7,7)*cos(6*lS) );
    dSnm66 = -0.0892/9*( (GM_Moon/gm)*((r_ref/rM)^7)*lgM(7,7)*sin(6*lM)...
           + (GM_Sun/gm)*((r_ref/rS)^7)*lgS(7,7)*sin(6*lS) );
    
    C(3,1) = C(3,1) + dCnm20;
    C(3,2) = C(3,2) + dCnm21;
    C(3,3) = C(3,3) + dCnm22;
    S(3,2) = S(3,2) + dSnm21;
    S(3,3) = S(3,3) + dSnm22;
    
    C(4,1) = C(4,1) + dCnm30;
    C(4,2) = C(4,2) + dCnm31;
    C(4,3) = C(4,3) + dCnm32;
    C(4,4) = C(4,4) + dCnm33;
    S(4,2) = S(4,2) + dSnm31;
    S(4,3) = S(4,3) + dSnm32;
    S(4,4) = S(4,4) + dSnm33;
    
    C(5,1) = C(5,1) + dCnm40;
    C(5,2) = C(5,2) + dCnm41;
    C(5,3) = C(5,3) + dCnm42;
    C(5,4) = C(5,4) + dCnm43;
    C(5,5) = C(5,5) + dCnm44;
    S(5,2) = S(5,2) + dSnm41;
    S(5,3) = S(5,3) + dSnm42;
    S(5,4) = S(5,4) + dSnm43;
    S(5,5) = S(5,5) + dSnm44;
    
    C(6,1) = C(6,1) + dCnm50;
    C(6,2) = C(6,2) + dCnm51;
    C(6,3) = C(6,3) + dCnm52;
    C(6,4) = C(6,4) + dCnm53;
    C(6,5) = C(6,5) + dCnm54;
    C(6,6) = C(6,6) + dCnm55;
    S(6,2) = S(6,2) + dSnm51;
    S(6,3) = S(6,3) + dSnm52;
    S(6,4) = S(6,4) + dSnm53;
    S(6,5) = S(6,5) + dSnm54;
    S(6,6) = S(6,6) + dSnm55;
    
    C(7,1) = C(7,1) + dCnm60;
    C(7,2) = C(7,2) + dCnm61;
    C(7,3) = C(7,3) + dCnm62;
    C(7,4) = C(7,4) + dCnm63;
    C(7,5) = C(7,5) + dCnm64;
    C(7,6) = C(7,6) + dCnm65;
    C(7,7) = C(7,7) + dCnm66;
    S(7,2) = S(7,2) + dSnm61;
    S(7,3) = S(7,3) + dSnm62;
    S(7,4) = S(7,4) + dSnm63;
    S(7,5) = S(7,5) + dSnm64;
    S(7,6) = S(7,6) + dSnm65;
    S(7,7) = S(7,7) + dSnm66;    
end

% Body-fixed position 
r_bf = E * r;

% Auxiliary quantities
d = norm(r_bf);                     % distance
latgc = asin(r_bf(3)/d);
lon = atan2(r_bf(2),r_bf(1));

[pnm, dpnm] = Legendre(AuxParam.n,AuxParam.m,latgc);

dUdr = 0;
dUdlatgc = 0;
dUdlon = 0;
q3 = 0; q2 = q3; q1 = q2;
for n=0:AuxParam.n
    b1 = (-gm/d^2)*(r_ref/d)^n*(n+1);
    b2 =  (gm/d)*(r_ref/d)^n;
    b3 =  (gm/d)*(r_ref/d)^n;
    for m=0:AuxParam.m
        q1 = q1 + pnm(n+1,m+1)*(C(n+1,m+1)*cos(m*lon)+S(n+1,m+1)*sin(m*lon));
        q2 = q2 + dpnm(n+1,m+1)*(C(n+1,m+1)*cos(m*lon)+S(n+1,m+1)*sin(m*lon));
        q3 = q3 + m*pnm(n+1,m+1)*(S(n+1,m+1)*cos(m*lon)-C(n+1,m+1)*sin(m*lon));
    end
    dUdr     = dUdr     + q1*b1;
    dUdlatgc = dUdlatgc + q2*b2;
    dUdlon   = dUdlon   + q3*b3;
    q3 = 0; q2 = q3; q1 = q2;
end

% Body-fixed acceleration
r2xy = r_bf(1)^2+r_bf(2)^2;

ax = (1/d*dUdr-r_bf(3)/(d^2*sqrt(r2xy))*dUdlatgc)*r_bf(1)-(1/r2xy*dUdlon)*r_bf(2);
ay = (1/d*dUdr-r_bf(3)/(d^2*sqrt(r2xy))*dUdlatgc)*r_bf(2)+(1/r2xy*dUdlon)*r_bf(1);
az =  1/d*dUdr*r_bf(3)+sqrt(r2xy)/d^2*dUdlatgc;

a_bf = [ax ay az]';

% Inertial acceleration 
a = E'*a_bf;

