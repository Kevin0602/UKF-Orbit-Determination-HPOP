%--------------------------------------------------------------------------
% 
%                   Orbital Perturbations
%
% References:
% Montenbruck O., Gill E.; Satellite Orbits: Models, Methods and 
% Applications; Springer Verlag, Heidelberg; Corrected 3rd Printing (2005).
%
% Montenbruck O., Pfleger T.; Astronomy on the Personal Computer; Springer 
% Verlag, Heidelberg; 4th edition (2000).
%
% Seeber G.; Satellite Geodesy; Walter de Gruyter, Berlin, New York; 2nd
% completely revised and extended edition (2003).
%
% Vallado D. A; Fundamentals of Astrodynamics and Applications; McGraw-Hill
% , New York; 3rd edition(2007).
%
% http://sol.spacenvironment.net/jb2008/
%
% http://ssd.jpl.nasa.gov/?ephemerides
%   
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------

format long g


global Cnm Snm AuxParam eopdata SOLdata DTCdata APdata n_eqn PC

SAT_Const

load DE405Coeff.mat
PC = DE405Coeff;

Cnm = zeros(71,71);
Snm = zeros(71,71);
fid = fopen('egm96','r');
for n=0:70
    for m=0:n
        temp = fscanf(fid,'%d %d %f %f %f %f',[6 1]);        
        Cnm(n+1,m+1) = temp(3);
        Snm(n+1,m+1) = temp(4);
    end
end

% read Earth orientation parameters
fid = fopen('eopc1962-now.txt','r');
%  ----------------------------------------------------------------------------------------------------
% |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
% |(0h UTC)           "         "          s          s          "        "          "         "     s 
%  ----------------------------------------------------------------------------------------------------
eopdata = fscanf(fid,'%i %d %d %i %f %f %f %f %f %f %f %f %i',[13 inf]);
eopdata = load('eopc1962-now.txt');
eopdata = eopdata';
fclose(fid);

% read space weather data
fid = fopen('SOLFSMY.txt','r');
%  ------------------------------------------------------------------------
% | YYYY DDD   JulianDay  F10   F81c  S10   S81c  M10   M81c  Y10   Y81c
%  ------------------------------------------------------------------------
SOLdata = fscanf(fid,'%d %d %f %f %f %f %f %f %f %f %f',[11 inf]);
fclose(fid);

% READ GEOMAGNETIC STORM DTC VALUE
fid = fopen('DTCFILE.txt','r');
%  ------------------------------------------------------------------------
% | DTC YYYY DDD   DTC1 to DTC24
%  ------------------------------------------------------------------------
DTCdata = fscanf(fid,'%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d',[26 inf]);
fclose(fid);

% read space weather data
fid = fopen('SOLRESAP.txt','r');
%  ------------------------------------------------------------------------
% | YYYY DDD  F10 F10B Ap1 to Ap8
%  ------------------------------------------------------------------------
APdata = fscanf(fid,'%d %d %f %f %f %f %f %f %f %f %f',[12 inf]);
fclose(fid);

% model parameters
AuxParam = struct('Mjd_UTC',0,'area_solar',0,'area_drag',0,'mass',0,'Cr',0,...
                  'Cd',0,'n',0,'m',0,'sun',0,'moon',0,'sRad',0,'drag',0,...
                  'planets',0,'SolidEarthTides',0,'OceanTides',0,'Relativity',0);

% epoch state (Envisat)
fid = fopen('InitialState.txt','r');

tline = fgetl(fid);
if ~ischar(tline)
    exit
end
year = str2num(tline(1:4));
mon = str2num(tline(6:7));
day = str2num(tline(9:10));
hour = str2num(tline(12:13));
min = str2num(tline(15:16));
sec = str2num(tline(18:23));
Y0(1) = str2num(tline(29:40));
Y0(2) = str2num(tline(42:53));
Y0(3) = str2num(tline(55:66));
Y0(4) = str2num(tline(68:79));
Y0(5) = str2num(tline(81:92));
Y0(6) = str2num(tline(94:105));

% epoch
Mjd_UTC = Mjday(year, mon, day, hour, min, sec);
Y0 = Y0';
Y0 = ECEF2ECI(Mjd_UTC, Y0');

tline = fgetl(fid);
AuxParam.area_solar = str2num(tline(49:53));
tline = fgetl(fid);
AuxParam.area_drag = str2num(tline(38:41));
tline = fgetl(fid);
AuxParam.mass = str2num(tline(19:24));
tline = fgetl(fid);
AuxParam.Cr = str2num(tline(5:7));
tline = fgetl(fid);
AuxParam.Cd = str2num(tline(5:7));

fclose(fid);

AuxParam.Mjd_UTC  = Mjd_UTC;
AuxParam.n       = 40;
AuxParam.m       = 40;
AuxParam.sun     = 1;
AuxParam.moon    = 1;
AuxParam.planets = 1;
AuxParam.sRad    = 1;
AuxParam.drag    = 1;
AuxParam.SolidEarthTides = 1;
AuxParam.OceanTides = 1;
AuxParam.Relativity = 1;

Mjd0  = Mjd_UTC;
n_eqn = 6;

% Step   = 60;   % [s]
% %N_Step = 1588; % 26.47hours
% N_Step = 20;
% 
% Eph = Ephemeris(Y0, N_Step, Step);

% fid = fopen('SatelliteStates.txt','w');
% 
% for i=1:N_Step+1
%     [year,mon,day,hr,min,sec] = invjday(Mjd0+Eph(i,1)/86400+2400000.5);
%     fprintf(fid,'  %4d/%2.2d/%2.2d  %2d:%2d:%6.3f',year,mon,day,hr,min,sec);
%     fprintf(fid,'  %14.3f%14.3f%14.3f%12.3f%12.3f%12.3f\n',...
%             Eph(i,2),Eph(i,3),Eph(i,4),Eph(i,5),Eph(i,6),Eph(i,7));
% end
% 
% fclose(fid);
% 
% [n, m] = size(Eph);
% Eph_ecef = zeros(n,m);
% 
% for i=1:n
%     Eph_ecef(i,1) = Eph(i,1);
%     Eph_ecef(i,2:7) = ECI2ECEF(Mjd0+Eph_ecef(i,1)/86400, Eph(i,2:7));    
% end
% 
% True_EnvisatStates
% dd = True_Eph(1:N_Step+1,:)-Eph_ecef(:,2:7);
% 
% % Plot orbit in ECI reference
% figure(1)
% plot3(Eph(:,2),Eph(:,3),Eph(:,4),'o-r')
% grid;
% title('Orbit ECI (inertial) (m)')
% 
% % Plot orbit in ECEF reference
% figure(2)
% plot3(Eph_ecef(:,2),Eph_ecef(:,3),Eph_ecef(:,4),'-')
% title('Orbit ECEF (m)')
% xlabel('X');ylabel('Y');zlabel('Z');
% grid
% 
% % Plot Discrepancy of Precise and Propagated orbits
% figure(3)
% subplot(3,1,1);
% plot(dd(:,1));
% title('Discrepancy of Precise and Propagated Envisat Positions for 26.47 hours');
% axis tight
% xlabel('Time')
% ylabel('dX[m]')
% hold on
% subplot(3,1,2);
% plot(dd(:,2));
% axis tight
% xlabel('Time')
% ylabel('dY[m]')
% subplot(3,1,3);
% plot(dd(:,3));
% axis tight
% xlabel('Time')
% ylabel('dZ[m]')



