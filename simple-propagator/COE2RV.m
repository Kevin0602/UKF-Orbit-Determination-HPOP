function [X,V] = COE2RV(a,e,i,OMEGA,w,M)
%% DESCRIPTION:
%
%       Written by:           Tyler Reid (tyreid@stanford.edu)
%       Lab:                  Stanford GPS Lab
%       Project Title:        Arctic Navigation / WAAS
%       Project Start Date:   March 28, 2011
%       Last updated:         March 31, 2011
%
% -------------------------------------------------------------------------
% FUNCTION DESCRIPTION
%
% Based on Vallado (2007) Algorithm 10
%
% This algorithm will compute the Earth Centered Inertial (ECI) position
% and velocity vectors which are equivalent to the given classical orbital
% elements
% 
% -------------------------------------------------------------------------
% INPUT:
%   
%       a = semi-major axis                          [length]*
%       e = eccentrity                               [rad]    
%       i = inclination                              [rad]
%   OMEGA = right ascension of the ascending node    [rad]
%       w = argument of perigee                      [rad]
%       M = mean anomaly                             [rad]
%
% ------------------------------------------------------------------------- 
% OUTPUT:
%      
%       X = ECI position vector of the spacecraft    [length]*
%       V = ECI veloicty vector of the spacecraft    [length / time]*
%  
% -------------------------------------------------------------------------
% NOTES:
%
% * This quantity can be expressed in either m or km or etc as long
%   as the global value of mu (the Earth's gravitational parameter) is in
%   consitant units.
%
%% GLOBAL VARIABLES USED

global mu

%% MAIN ALGORITHM

% Handle special cases:

% Circular equitorial orbit.
if e == 0 && i == 0
    w     = 0;
    OMEGA = 0;
  
% Circular inclined orbit.
elseif e == 0
    w     = 0;
    
% Elliptical equitorial.
elseif i == 0 
    OMEGA = 0;       
end

% Compute the semi-latus rectum.
p = a*(1-e^2);

% Convert mean anomaly to true anomaly.
% First, compute the eccentric anomaly.
Ea = Keplers_Eqn(M,e);

% Compute the true anomaly f.
y = sin(Ea)*sqrt(1-e^2)/(1-e*cos(Ea));
x = (cos(Ea)-e)/(1-e*cos(Ea));

f = atan2(y,x);

% Define the position vector in perifocal PQW coordinates.
X(1,1) = p*cos(f)/(1+e*cos(f));
X(2,1) = p*sin(f)/(1+e*cos(f));
X(3,1) = 0;

% Define the velocity vector in perifocal PQW coordinates.
V(1,1) = -sqrt(mu/p)*sin(f);
V(2,1) =  sqrt(mu/p)*(e+cos(f));
V(3,1) =  0;

% Define Transformation Matrix To IJK.
Trans(1,1) =  cos(OMEGA)*cos(w)-sin(OMEGA)*sin(w)*cos(i);
Trans(1,2) = -cos(OMEGA)*sin(w)-sin(OMEGA)*cos(w)*cos(i);
Trans(1,3) =  sin(OMEGA)*sin(i);

Trans(2,1) =  sin(OMEGA)*cos(w)+cos(OMEGA)*sin(w)*cos(i);
Trans(2,2) = -sin(OMEGA)*sin(w)+cos(OMEGA)*cos(w)*cos(i);
Trans(2,3) = -cos(OMEGA)*sin(i);

Trans(3,1) =  sin(w)*sin(i);
Trans(3,2) =  cos(w)*sin(i);
Trans(3,3) =  cos(i);

% Transform to the ECI coordinate system.
X = Trans*X;
V = Trans*V;
