%--------------------------------------------------------------------------
%
% Relativisty: Computes the perturbational acceleration due to relativistic
%              effects
%
% Inputs:
%   r           Satellite position vector
%   v           Satellite velocity vector
% 
% Output:
%   a    		Acceleration (a=d^2r/dt^2)
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
function a = Relativity(r, v)

SAT_Const

% Relative position vector of satellite w.r.t. point mass 
r_Sat = norm(r);
v_Sat = norm(v);

% Acceleration 
a = GM_Earth/(c_light^2*r_Sat^3)*((4*GM_Earth/r_Sat-v_Sat^2)*r+4*dot(r,v)*v);

