%-------------------------------------------------------------------
%------------------------------- sgp4 ------------------------------
%-------------------------------------------------------------------
function [pos, vel] = sgp4(tsince, satdata)

ae = 1.0;
tothrd = (2.0/3.0);
XJ3 = -2.53881e-6;
e6a = 1.0E-6;
xkmper = 6378.135;
ge = 398600.8; % Earth gravitational constant
CK2 = (1.0826158e-3 / 2.0);
CK4 = (-3.0 * -1.65597e-6 / 8.0);

% Constants
s = ae + 78 / xkmper;
qo = ae + 120 / xkmper;
xke = sqrt((3600.0 * ge) / (xkmper^3));
qoms2t = ((qo - s)^2)^2;
temp2 = xke / (satdata.xno);
a1 = temp2^tothrd;
cosio = cos (satdata.xincl);
theta2 = (cosio^2);
x3thm1 = 3.0 * theta2 - 1.0;
eosq = (satdata.eo^2);
betao2 = 1.0 - eosq;
betao = sqrt(betao2);
del1 = 1.5 * CK2 * x3thm1 / ((a1^2) * betao * betao2);
ao = a1 * ( 1.0 - del1*((1.0/3.0) + del1 * (1.0 + (134.0/81.0) * del1)));
delo = 1.5 * CK2 * x3thm1 / ((ao^2) * betao * betao2);
xnodp = (satdata.xno)/(1.0 + delo);
aodp = ao/(1.0 - delo);
% Initialization
% For perigee less than 220 kilometers, the isimp flag is set and
% the equations are truncated to linear variation in sqrt a and
% quadratic variation in mean anomaly.  Also, the c3 term, the
% delta omega term, and the delta m term are dropped.
isimp = 0;
if ((aodp * (1.0 - satdata.eo)/ ae) < (220.0/xkmper + ae))
    isimp = 1;
end
% For perigee below 156 km, the values of s and qoms2t are altered.
s4 = s;
qoms24 = qoms2t;
perige = (aodp * (1.0 - satdata.eo) - ae) * xkmper;
if (perige < 156)
    s4 = perige - 78.0;
    if (perige <= 98)
        s4 = 20.0;
    end
    qoms24 = (((120.0 - s4) * ae / xkmper)^4.0);
    s4 = s4 / xkmper + ae;
end
pinvsq = 1.0 / ( (aodp^2) * (betao2^2) );
tsi = 1.0 / (aodp - s4);
eta = aodp * (satdata.eo) * tsi;
etasq = (eta^2);
eeta = (satdata.eo) * eta;
psisq = abs( 1.0 - etasq);
coef = qoms24 * (tsi^4.0);
coef1 = coef / (psisq^3.5);
c2 = coef1 * xnodp * (aodp * (1.0 + 1.5 * etasq + eeta * (4.0 + etasq)) + 0.75 * CK2 * tsi / psisq * x3thm1 * (8.0 + 3.0 * etasq * (8.0 + etasq)));
c1 = (satdata.bstar) * c2;
sinio = sin(satdata.xincl);
a3ovk2 = -XJ3 / CK2 * (ae^3.0);
c3 = coef * tsi * a3ovk2 * xnodp * ae * sinio / (satdata.eo);
x1mth2 = 1.0 - theta2;
c4 = 2.0 * xnodp * coef1 * aodp * betao2 * ( eta * (2.0 + 0.5 * etasq) + (satdata.eo) * (0.5 + 2.0 * etasq) - 2.0 * CK2 * tsi / (aodp * psisq) * ( -3.0 * x3thm1 * ( 1.0 - 2.0 * eeta + etasq * (1.5 - 0.5*eeta)) + 0.75 * x1mth2 * (2.0 * etasq - eeta * (1.0 + etasq)) * cos(2.0 * (satdata.omegao))));
c5 = 2.0 * coef1 * aodp * betao2 * (1.0 + 2.75 * (etasq + eeta) + eeta * etasq);
theta4 = (theta2^2);
temp1 = 3.0 * CK2 * pinvsq * xnodp;
temp2 = temp1 * CK2 * pinvsq;
temp3 = 1.25 * CK4 * pinvsq * pinvsq * xnodp;
xmdot = xnodp + 0.5 * temp1 * betao * x3thm1 + 0.0625 * temp2 * betao * (13.0 - 78.0 * theta2 + 137.0 * theta4);
x1m5th = 1.0 - 5.0 * theta2;
omgdot = -0.5 * temp1 * x1m5th + 0.0625 * temp2 * (7.0 - 114.0 * theta2 + 395.0 * theta4) + temp3 * (3.0 - 36.0 * theta2 + 49.0 * theta4);
xhdot1 = -temp1 * cosio;
xnodot = xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * theta2) + 2.0 * temp3 * (3.0 - 7.0 * theta2)) * cosio;
omgcof = (satdata.bstar) * c3 * cos(satdata.omegao);
xmcof = -(2.0/3.0) * coef * (satdata.bstar) * ae / eeta;
xnodcf = 3.5 * betao2 * xhdot1 * c1;
t2cof = 1.5 * c1;
xlcof = 0.125 * a3ovk2 * sinio * (3.0 + 5.0 * cosio) / (1.0 + cosio);
aycof = 0.25 * a3ovk2 * sinio;
delmo = ((1.0 + eta * cos(satdata.xmo))^3);
sinmo = sin(satdata.xmo);
x7thm1 = 7.0 * theta2 - 1.0;
if (isimp==0)	
    c1sq = (c1^2);
    d2 = 4.0 * aodp * tsi * c1sq;
    temp = d2 * tsi * c1 / 3.0;
    d3 = (17.0 * aodp + s4)*temp;
    d4 = 0.5 * temp * aodp * tsi * (221.0 * aodp + 31.0 * s4) * c1;
    t3cof = d2 + 2.0*c1sq;
    t4cof = 0.25 * (3.0 * d3 + c1 * (12.0 * d2 + 10.0 * c1sq));
    t5cof = 0.2 * (3.0 * d4 + 12.0 * c1 * d3 + 6.0 * d2 * d2 + 15.0 * c1sq * (2.0 * d2 + c1sq));
end
% Update for secular gravity and atmospheric drag.
xmdf = satdata.xmo + xmdot * tsince;
omgadf = satdata.omegao + omgdot * tsince;
xnoddf = satdata.xnodeo + xnodot * tsince;
omega = omgadf;
xmp = xmdf;
tsq = (tsince^2);
xnode = xnoddf + xnodcf * tsq;
tempa = 1.0 - c1 * tsince;
tempe = (satdata.bstar) * c4 * tsince;
templ = t2cof * tsq;
if (isimp == 0)
    delomg = omgcof * tsince;
    delm = xmcof*(((1.0 + eta * cos(xmdf))^ 3.0) - delmo);
    temp = delomg + delm;
    xmp = xmdf + temp;
    omega = omgadf - temp;
    tcube = tsq * tsince;
    tfour = tsince * tcube;
    tempa = tempa - d2 * tsq - d3 * tcube - d4 * tfour;
    tempe = tempe + (satdata.bstar) * c5 * (sin(xmp) - sinmo);
    templ = templ + t3cof * tcube + tfour * (t4cof + tsince * t5cof);
end
a = aodp * (tempa^2);
e = (satdata.eo) - tempe;
xl = xmp + omega + xnode + xnodp*templ;
beta = sqrt(1.0 - (e^2));
xn = xke / (a^1.5);
% Long period periodics
axn = e * cos(omega);
temp = 1.0 / (a * (beta^2));
xll = temp * xlcof * axn;
aynl = temp * aycof;
xlt = xl + xll;
ayn = e * sin(omega) + aynl;
% Solve Kepler's Equation
capu = fmod2p(xlt - xnode);
temp2 = capu;
i=1;
while(1)
    sinepw = sin(temp2);
    cosepw = cos(temp2);
    temp3 = axn * sinepw;
    temp4 = ayn * cosepw;
    temp5 = axn * cosepw;
    temp6 = ayn * sinepw;
    epw = (capu - temp4 + temp3 - temp2) / (1.0 - temp5 - temp6) + temp2;
    temp7 = temp2;
    temp2 = epw;
    i = i+1;
	if ((i>10) || (abs(epw - temp7) <= e6a))
		break
	end
end

% Short period preliminary quantities
ecose = temp5 + temp6;
esine = temp3 - temp4;
elsq = (axn^2) + (ayn^2);
temp = 1.0 - elsq;
pl = a * temp;
r = a * (1.0 - ecose);
temp1 = 1.0 / r;
rdot = xke * sqrt(a) * esine * temp1;
rfdot = xke * sqrt(pl) * temp1;
temp2 = a * temp1;
betal = sqrt(temp);
temp3 = 1.0 / (1.0 + betal);
cosu = temp2 * (cosepw - axn + ayn * esine * temp3);
sinu = temp2 * (sinepw - ayn - axn * esine * temp3);
u = actan(sinu, cosu);
sin2u = 2.0 * sinu * cosu;
cos2u = 2.0 * (cosu^2) - 1.0;
temp = 1.0 / pl;
temp1 = CK2 * temp;
temp2 = temp1 * temp;
% Update for short periodics
rk = r * (1.0 - 1.5 * temp2 * betal * x3thm1) + 0.5 * temp1 * x1mth2 * cos2u;
uk = u - 0.25 * temp2 * x7thm1 * sin2u;
xnodek = xnode + 1.5 * temp2 * cosio * sin2u;
xinck = (satdata.xincl) + 1.5 * temp2 * cosio * sinio * cos2u;
rdotk = rdot - xn * temp1 * x1mth2 * sin2u;
rfdotk = rfdot + xn * temp1 * (x1mth2 * cos2u + 1.5 * x3thm1);
% Orientation vectors
MV.v(1) = -sin(xnodek) * cos(xinck);
MV.v(2) = cos(xnodek) * cos(xinck);
MV.v(3) = sin(xinck);

NV.v(1) = cos(xnodek);
NV.v(2) = sin(xnodek);
NV.v(3) = 0;

for i=1:3
	UV.v(i) = MV.v(i) * sin(uk) + NV.v(i) * cos(uk);
	VV.v(i) = MV.v(i) * cos(uk) - NV.v(i) * sin(uk);
end

% position + velocity
for i=1:3
	pos.v(i) = rk * UV.v(i);
	vel.v(i) = rdotk * UV.v(i) + rfdotk * VV.v(i);
end

[pos, vel] = Convert_Sat_State(pos, vel);

