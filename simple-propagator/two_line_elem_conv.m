function [OrbitData] = two_line_elem_conv(file_name,SV)
%% DESCRIPTION
%
%       Written by:           Tyler Reid
%       Lab:                  Stanford GPS Lab
%       Project Title:        Arctic Navigation / WAAS
%       Project Start Date:   March 28, 2011
%       Last updated:         January 14, 2013
%
% -------------------------------------------------------------------------
% FUNCTION DESCRIPTION
%
% Given a two line element (TLE) file from NORAD and the satellite vehicle
% name, output the Keplerian orbital element, and date of ephemeris.
%
% -------------------------------------------------------------------------
% INPUT:
%
%       file_name = string which contains the file name of the TLE data
%                   e.g. 'sbas.txt'
%       SV        = string which contains the name of satellite in question
%                   'all' returns all in the given text file.
%
% ------------------------------------------------------------------------- 
% OUTPUT:
%
%       OrbitData = data structure of the form:
%        
%           OrbitData.ID          = NORAD spacecraft ID
%           OrbitData.designation = NORAD spacecraft designation
%           OrbitData.i           = inclination [degrees]
%           OrbitData.RAAN        = RAAN [degrees]
%           OrbitData.omega       = argument of perigee [degrees]  
%           OrbitData.M           = Mean anomaly at epoch [degrees]
%           OrbitData.a           = semi-major axis [meters]
%           OrbitData.e           = eccentricity [unitless]
%           OrbitData.date        = string containing the UTC epoch which
%                                   the TLE is referenced to in the form:
%                                   "dd-mmm-yyyy HH:MM:SS.SSS"
%           OrbitData.BC          = NORAD SPG-4 Ballastic coefficient
%
% -------------------------------------------------------------------------
% NOTES:
%
%       (1) See page 115 of Vallado's 'Fundamentals of Astrodynamics and
%       Applications'
%
%% IMPLEMENTATION

% Earth Gravitational Parameter mu = G*M_earth.
mu =       3.986004418e14;      % [m^3/s^2]

% Open the specified.
fid      = fopen(file_name);
txt_data = textscan(fid,'%s %s %s %s %s %s %s %s %s');

% Find the satellite in question.
if strcmp(SV,'all') == 1
    % Find the total number of satellites in the file.
    numSV = length(txt_data{1})/3;
    
    % Initialize array.
    ind = zeros(1,numSV);
    
    count = 1;
    for i = 1:numSV
        % Take every 3rd line.
        ind(i) = count;
        
        OrbitData.ID(i) = txt_data{1}(count);
        OrbitData.designation(i) = txt_data{2}(count);
        OrbitData.PRN(i) = txt_data{3}(count);
        
        count  = count + 3;
    end
else
    % TODO - code to find specified satellites.
    % This was never done as it was not needed. 
end

% Find the two lines corresponding to the spacecraft in question.
for j = 1:length(ind)
    % Find the first line of the first satellite.
    index = ind(j);
    
    % Translate two line element data into obital elements.
    OrbitData.i(j)     = str2double(txt_data{1,3}{index+2});      % [deg]
    OrbitData.RAAN(j)  = str2double(txt_data{1,4}{index+2});      % [deg]
    OrbitData.omega(j) = str2double(txt_data{1,6}{index+2});      % [deg]
    OrbitData.M(j)     = str2double(txt_data{1,7}{index+2});      % [deg]
    n                  = str2double(txt_data{1,8}{index+2});      % [rev/day]
    n                  = n*2*pi/24/60/60;                         % [rad/s]
    OrbitData.a(j)     = ( mu / n^2 )^(1/3);                      % [m]
    OrbitData.e(j)     = str2double(txt_data{1,5}{index+2})*1e-7; % [-]
    
    % Compute the UTC date / time.
    temp2             = txt_data{1,4}{index+1};
    yy                = str2double(temp2(1:2));
    yyyy              = 2000 + yy;
    start             = datenum([yyyy 0 0 0 0 0]);
    secs              = str2double(temp2(3:length(temp2)))*24*3600-2*24*3600;
    date1             = datevec(addtodate(start,floor(secs),'second'));
    remainder         = [0 0 0 0 0 mod(secs,1)];
    OrbitData.date{j} = datestr(date1+remainder,'dd-mmm-yyyy HH:MM:SS.FFF');
    
    % Compute ballistic coefficient in SI units.
    temp3 = txt_data{1,7}{index+1};
    if length(temp3) == 7
        base  = str2double(temp3(1:5));
        expo  = str2double(temp3(6:7));
    elseif length(temp3) == 8
        base  = str2double(temp3(2:6));
        expo  = str2double(temp3(7:8));
    else
        fprintf('Error in ballistic coefficient calculation\n')
        error('end program')
    end
    
    Bstar = base*10^expo;
    OrbitData.BC(j)    = 1/12.741621/Bstar; % [kg/m^2]
        
end
