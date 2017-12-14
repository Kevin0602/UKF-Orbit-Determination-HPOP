%% PHYSICAL CONSTANTS
%
%       Written by:           Tyler Reid
%       Lab:                  Stanford GPS Lab
%       Project Title:        Arctic Navigation / WAAS
%       Project Start Date:   March 28, 2011
%       Last updated:         April 14, 2011
%
% -------------------------------------------------------------------------
% DESCRIPTION
%
% Loading this file enters all of the physical constants needed to perform
% all of the neccessary calculations into the global workspace.  All
% constants are given in SI units.
%
%% DEFINE CONSTANTS

global mu omega_e R_e J_2 Earth_E2 

% Earth's Mean Equitorial Radius.
R_e =      6.37813649e6;        % [m]

% Earth Gravitational Parameter mu = G*M_earth.
mu =       3.986004418e14;      % [m^3/s^2]

% Mean Angular Velocity of the Earth.
omega_e =  7.29211585530e-5;    % [rad/s]

% J_2 - second zonal harmonic.
J_2 =      1.0826300e-3;        % [-]

% Earth's shape - eccentricity^2.
Earth_E2 = 0.006694385000;      % [-]


