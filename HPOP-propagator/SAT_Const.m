%--------------------------------------------------------------------------
%
% SAT_Const: Definition of astronomical and mathematical constants
%
% Last modified:   2015/08/12   M. Mahooti
% 
%--------------------------------------------------------------------------

% Mathematical constants
pi2       = 2*pi;                % 2pi
Rad       = pi/180;              % Radians per degree
Deg       = 180/pi;              % Degrees per radian
Arcs      = 3600*180/pi;         % Arcseconds per radian

% General
MJD_J2000 = 51544.5;             % Modified Julian Date of J2000
T_B1950   = -0.500002108;        % Epoch B1950
c_light   = 299792457.999999984; % Speed of light  [m/s]; DE405
AU        = 149597870691.000015; % Astronomical unit [m]; DE405

% Physical parameters of the Earth, Sun and Moon

% Equatorial radius and flattening
R_Earth   = 6378.137e3;          % Earth's radius [m]; WGS-84
f_Earth   = 1/298.257223563;     % Flattening; WGS-84   
R_Sun     = 696000e3;            % Sun's radius [m]; DE405
R_Moon    = 1738e3;              % Moon's radius [m]; DE405

% Earth rotation (derivative of GMST at J2000; differs from inertial period by precession)
omega_Earth = 15.04106717866910/3600*Rad;  % [rad/s]; WGS-84

% Gravitational coefficients
GM_Earth   = 398600.4418e9;                % [m^3/s^2]; WGS-84
GM_Sun     = 132712440017.9870e9;          % [m^3/s^2]; DE405
GM_Moon    = GM_Earth/81.3005600000000044; % [m^3/s^2]; DE405
GM_Mercury = 22032.08048641792e9;          % [m^3/s^2]; DE405
GM_Venus   = 324858.5988264596e9;          % [m^3/s^2]; DE405
GM_Mars    = 42828.31425806710e9;          % [m^3/s^2]; DE405
GM_Jupiter = 126712767.8577960e9;          % [m^3/s^2]; DE405
GM_Saturn  = 37940626.06113726e9;          % [m^3/s^2]; DE405
GM_Uranus  = 5794549.007071872e9;          % [m^3/s^2]; DE405
GM_Neptune = 6836534.063879259e9;          % [m^3/s^2]; DE405
GM_Pluto   = 981.6008877070042e9;          % [m^3/s^2]; DE405
              
% Solar radiation pressure at 1 AU 
P_Sol      = 1367/c_light;       % [N/m^2] (1367 W/m^2); IERS 96

