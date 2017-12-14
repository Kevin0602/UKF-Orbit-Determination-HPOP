%--------------------------------------------------------------------------
%
% AccelDrag: Computes the acceleration due to the atmospheric drag in
%            ICRF/EME2000 system
% 
% Inputs:
%   dens        Atmospheric Density [kg/m^3]
%   r           Satellite position vector in the inertial system [m]
%   v           Satellite velocity vector in the inertial system [m/s]
%   T           Transformation matrix to true-of-date inertial system
%   Area        Cross-section [m^2]
%   mass        Spacecraft mass [kg]
%   CD          Drag coefficient
%   Omega		Angular velocity of the Earth
%
% Output:
%   a           Acceleration (a=d^2r/dt^2) [m/s^2]
%
% Last modified:   2015/08/12   M. Mahooti
% 
%--------------------------------------------------------------------------
function a = AccelDrag(dens, r, v, T, Area, mass, CD, Omega)

% Earth angular velocity vector [rad/s]
omega = [0, 0, Omega]';

% Position and velocity in true-of-date system
r_tod = T * r;
v_tod = T * v;
  
% Velocity relative to the Earth's atmosphere
v_rel = v_tod - cross(omega, r_tod);
v_abs = norm(v_rel);

% Acceleration
a_tod = -0.5*CD*(Area/mass)*dens*v_abs*v_rel;

a = T' * a_tod;

