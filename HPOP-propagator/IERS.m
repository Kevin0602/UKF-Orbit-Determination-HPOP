%--------------------------------------------------------------------------
%
% IERS: Management of IERS time and polar motion data
%  
% Last modified:   2015/08/12   M. Mahooti
% 
%--------------------------------------------------------------------------
function [UT1_UTC,TAI_UTC,x_pole,y_pole,ddpsi,ddeps] = IERS(eop,Mjd_UTC,interp)

if (nargin == 2)
   interp = 'n';
end
Arcs = 3600*180/pi;  % Arcseconds per radian

if (interp =='l')
    % linear interpolation
    mj = (floor(Mjd_UTC));
    nop = length(eop);
    
    for i=1:nop
        if (mj==eop(4,i))
            preeop = eop(:,i);
            nexteop = eop(:,i+1);
            break
        end
    end
    
    mfme = 1440*(Mjd_UTC-floor(Mjd_UTC));
    fixf = mfme/1440;
    
    % Setting of IERS Earth rotation parameters
    % (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
    UT1_UTC = preeop(7)+(nexteop(7)-preeop(7))*fixf;
    % TAI_UTC = preeop(13)+(nexteop(13)-preeop(13))*fixf;
	TAI_UTC = preeop(13);
    x_pole  = preeop(5)+(nexteop(5)-preeop(5))*fixf;
    y_pole  = preeop(6)+(nexteop(6)-preeop(6))*fixf;
    ddpsi   = preeop(9)+(nexteop(9)-preeop(9))*fixf;
    ddeps   = preeop(10)+(nexteop(10)-preeop(10))*fixf;
    
    x_pole  = x_pole/Arcs;  % Pole coordinate [rad]
    y_pole  = y_pole/Arcs;  % Pole coordinate [rad]
    ddpsi   = ddpsi/Arcs;
    ddeps   = ddeps/Arcs;	
elseif (interp =='n')    
    mj = (floor(Mjd_UTC));
    nop = length(eop);
    
    for i=1:nop
        if (mj==eop(4,i))
            eop = eop(:,i);
            break;
        end
    end
    
    % Setting of IERS Earth rotation parameters
    % (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
    UT1_UTC = eop(7);      % UT1-UTC time difference [s]
    TAI_UTC = eop(13);     % TAI-UTC time difference [s]
    x_pole  = eop(5)/Arcs; % Pole coordinate [rad]
    y_pole  = eop(6)/Arcs; % Pole coordinate [rad]
    ddpsi   = eop(9)/Arcs;
    ddeps   = eop(10)/Arcs;
end

