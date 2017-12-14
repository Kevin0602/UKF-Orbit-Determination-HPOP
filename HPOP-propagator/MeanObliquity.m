%--------------------------------------------------------------------------
%
% MeanObliquity: Computes the mean obliquity of the ecliptic
%
% Input:
%  Mjd_TT    Modified Julian Date (Terrestrial Time)
% 
% Output:
%  MOblq     Mean obliquity of the ecliptic
%
% Last modified:   2015/08/12   M. Mahooti
% 
%--------------------------------------------------------------------------
function MOblq = MeanObliquity(Mjd_TT)

SAT_Const

T = (Mjd_TT-MJD_J2000)/36525;

MOblq = Rad*(23.43929111-(46.8150+(0.00059-0.001813*T)*T)*T/3600);

