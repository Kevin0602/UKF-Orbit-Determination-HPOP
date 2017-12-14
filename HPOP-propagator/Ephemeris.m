%--------------------------------------------------------------------------
%
% Ephemeris computation using variable-order Radau IIA integrator with 
% step-size control
%
% Last modified:   2015/08/12   M. Mahooti
% 
%--------------------------------------------------------------------------
function [Eph] = Ephemeris(Y0, N_Step, Step)

global n_eqn

Eph = zeros(N_Step,n_eqn+1);
Eph(1,1) = 0;
Eph(1,2:7) = Y0;

options = rdpset('RelTol',1e-13,'AbsTol',1e-16);

for i=1:N_Step	
	[~,yout] = radau(@Accel,[(i-1)*Step i*Step],Y0,options);
	Y0 = yout(end,:)';
    Eph(i+1,1) = i*Step;
    Eph(i+1,2:7) = Y0;
    

end


