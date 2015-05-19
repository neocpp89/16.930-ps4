function [psi,velx,vely,Gamma] = potential_trefftz(x,y,V,alpha,tparam)
%POTENTIAL_TREFFTZ Calculates the 2D potential flow for Trefftz airfoil.
%   [PSI,VELX,VELY,GAMMA]=POTENTIAL_TREFFTZ(X,Y,V,ALPHA,TPARAM)
%
%      PSI:       Value of stream function at input points
%      VELX:      x omponent of the velocity at input points
%      VELY:      y component of the velocity at input points
%      GAMMA:     Circultation. Lift Force= V*GAMMA
%      X:         x coordinates of input points 
%      Y:         y coordinates of input points 
%      V:         Free stream velocity magnitude (default=1)
%      ALPHA:     Angle of attack in degrees (default=0)
%      TPARAM:    Trefftz foil parameters
%                 tparam(1) = left x-shift of circle center 
%                             (trailing edge at (1,0)). (default=0.1)
%                 tparam(2) = y-shift of circle center. (default=0.05)
%                 tparam(3) = K-T exponent (=< 2) (2:Jukowski). (default=1.98)                      
%
% - Written by: J. Peraire
%
if nargin<3, V=1; end
if nargin<4, alpha=0; end
if nargin<5, tparam=[0.1,0.05,1.98]; end

x0 = tparam(1);
y0 = tparam(2);
n  = tparam(3);

rot = atan2(y0,1+x0);
r = sqrt((1+x0)^2 + y0^2);

zr = squeeze(x);
zi = squeeze(y);
z = complex(zr,zi);

%First calcualte an approximate camber line (ugly code to ensure we are in the correct branch)
cc = complex(-x0,y0);
th = 0:2*pi/120:2*pi;
xc = (1-cc)*exp(i*th');
wd = cc+xc;
zd = ((wd-1)./(wd+1)).^n;
wd = ((1+zd)./(1-zd))*n;
[xle,ii] = min(real(wd));
sup = spline(real(wd(1:ii)),imag(wd(1:ii)));
slo = spline(real(wd(ii:end)),imag(wd(ii:end)));

%Now K-T inverse
A = ((z-n)./(z+n));
anga = angle(A);
il = find((zr > xle) & (zr <= n) & (zi < 0.5*(ppval(sup,zr)+ppval(slo,zr))) & (anga > 1.5));
iu = find((zr > xle) & (zr <= n) & (zi > 0.5*(ppval(sup,zr)+ppval(slo,zr))) & (anga < -1.5));
anga = reshape(anga,numel(A),1);
anga(il) = anga(il) - 2*pi;
anga(iu) = anga(iu) + 2*pi;
anga = reshape(anga,size(A));

B = (abs(A).^(1/n)).*exp(i*anga/n);
v = (1+B)./(1-B);

%Scale back
w = (1/r)*exp(i*rot)*(v-complex(-x0,y0));

%Now we have a unit circle
alphef = pi*alpha/180+rot;

%Calculate circulation
dphidw = -V*(exp(-i*alphef) - 1./exp(-i*alphef));
dvortdw = i./(2*pi);
Gamma = -real(dphidw/dvortdw);

phi = -V*r*(w*exp(-i*alphef) + 1./(w*exp(-i*alphef)));
vort = i*Gamma*log(w)/(2*pi);
psi = imag(phi+vort);

%Find trailing edge
ii = find(abs(w-1)<1.e-6);
if numel(ii) > 0
    w = reshape(w,numel(w),1);
    w(ii) = complex(2,0);
    w = reshape(w,size(v));
end

dphidw = -V*r*(exp(-i*alphef) - 1./(w.*w*exp(-i*alphef)));
dvortdw = i*Gamma./(2*pi*w);

dwdv = (1/r)*exp(i*rot);
dvdB = 2./(1-B).^2;
dBdz = (1/n)*(abs(A).^((1-n)/n)).*exp(i*anga*(1-n)/n).*(1-A)./(z+n);

dphi = (dphidw+dvortdw).*dwdv.*dvdB.*dBdz;

if numel(ii) > 0
    dphi = reshape(dphi,numel(dphi),1);
    dphi(ii) = complex(0,0);
    dphi = reshape(dphi,size(v));
end

velx = -real(dphi);
vely =  imag(dphi);







