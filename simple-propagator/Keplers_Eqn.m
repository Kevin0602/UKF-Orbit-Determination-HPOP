function E_a = Keplers_Eqn(M,e)
%% DESCRIPTION:
%
%   AA 279 - SPACE MECHANICS 
%   
%   Problem Set  # 2
%   Problem      2.4
%   Written by:  Tyler Reid
%   Date:        April 14, 2011
%
% -------------------------------------------------------------------------
% FUNCTION DESCRIPTION
%
% Based on Vallado (2007) Algorithm 2
%
% Given the mean anomaly M, eccentricity e, compute the eccentric anomaly E
% 
% -------------------------------------------------------------------------
% INPUT:
%   
%       e = eccentrity                               [rad]    
%       M = mean anomaly                             [rad]
%
% ------------------------------------------------------------------------- 
%
% OUTPUT:
%      
%       E = eccentric anomaly                        [rad]
%
%% MAIN ALGORITHM

% Select initial guess.
if -pi<M<0 || M>pi
    E_a = M-e;
else
    E_a = M+e;
end

% Define tolerance.
tol  = 1e-12;

test = 999; % Dummy variable.

% Implement Newton's method.
while test > tol
    E_new = E_a + (M-E_a+e*sin(E_a))/(1-e*cos(E_a));
    test = abs(E_new - E_a);
    E_a = E_new;
end
    
