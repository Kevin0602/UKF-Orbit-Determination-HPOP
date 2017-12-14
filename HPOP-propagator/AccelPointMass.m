%--------------------------------------------------------------------------
%
% AccelPointMass: Computes the perturbational acceleration due to a point
%				  mass
%
% Inputs:
%   r           Satellite position vector 
%   s           Point mass position vector
%   GM          Gravitational coefficient of point mass
%
% Output:
%   a    		Acceleration (a=d^2r/dt^2)
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
function a = AccelPointMass(r, s, GM)

% Relative position vector of satellite w.r.t. point mass 
d = r - s;

% Acceleration 
a = -GM * ( d/(norm(d)^3) + s/(norm(s)^3) );

