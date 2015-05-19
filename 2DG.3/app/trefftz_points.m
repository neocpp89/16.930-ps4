function [x,y,chord] = trefftz_points(tparam,np)
%TREFFTZ_POINTS Calculates np points on Trefftz airfoil surface.
%   [X,Y,CHORD]=TREFFTZ_POINTS(TPARAM,NP)
%
%      X:         X coordinates of generated points 
%      Y:         Y coordinates of generated points 
%      CHORD:     Foil chord
%      TPARAM:    Trefftz foil parameters
%                 tparam(1) = left x-shift of circle center 
%                             (trailing edge at (1,0)). (default=0.1)
%                 tparam(2) = y-shift of circle center. (default=0.05)
%                 tparam(3) = K-T exponent (=< 2) (2:Jukowski). (default=1.98)
%      NP:        Number of points requested
%
% - Written by: J. Peraire
%
if nargin<1, tparam=[0.1,0.05,1.98]; end
if nargin<2, np=120; end

x0 = tparam(1);
y0 = tparam(2);
n  = tparam(3);

rot = atan2(y0,1+x0);
r = sqrt((1+x0)^2 + y0^2);

%K-T Transform
cc = complex(-x0,y0);
th = 0:2*pi/np:2*pi;
xc = (1-cc)*exp(i*th');
wd = cc+xc;
zd = ((wd-1)./(wd+1)).^n;
wd = ((1+zd)./(1-zd))*n;
xle = min(real(wd));
chord = n-xle;

x=real(wd);
y=imag(wd);





