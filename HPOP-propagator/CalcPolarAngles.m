%--------------------------------------------------------------------------
%
% Calculate polar components
%
% Last modified:   2015/08/12   M. Mahooti
%
%--------------------------------------------------------------------------
function [m_phi, m_theta, m_r] = CalcPolarAngles(m_Vec)

% Length of projection in x-y-plane:
rhoSqr = m_Vec(1) * m_Vec(1) + m_Vec(2) * m_Vec(2); 

% Norm of vector
m_r = sqrt(rhoSqr + m_Vec(3) * m_Vec(3));

% Azimuth of vector
if ( (m_Vec(1)==0) && (m_Vec(2)==0) )
    m_phi = 0;
else
    m_phi = atan2(m_Vec(2), m_Vec(1));
end
if ( m_phi < 0 )
    m_phi = m_phi + 2*pi;
end

% Altitude of vector
rho = sqrt( rhoSqr );
if ( (m_Vec(3)==0) && (rho==0) )
    m_theta = 0;
else
    m_theta = atan2(m_Vec(3), rho);
end

