function fn = wavei_roe(up,um,np,p,param,time)
%WAVEI Calculate Interface Lax Friedrichs Flux for the wave equation.
%   FN=WAVEI(UP,UM,NP,P,PARAM,TIME)
%
%      UP(np,3):   np plus states
%      UR(np,3):   np minus states
%      NP(np,2):   np normal plus vectors 
%      P(np,2):    np x,y coordinates
%      PARAM{1}:   Cell array containing either the wave speed c =PARAM{1}
%      TIME:       Not used
%      FN(np,3):   np normal fluxes (f plus)    
%                          
% - Written by: J. Peraire
%
c = param{1};
ca = abs(c);

fxl = -c*[up(:,3), zeros(size(up(:,1))), up(:,1)];
fyl = -c*[zeros(size(up(:,1))), up(:,3), up(:,2)];
fxr = -c*[um(:,3), zeros(size(um(:,1))), um(:,1)];
fyr = -c*[zeros(size(um(:,1))), um(:,3), um(:,2)];
fav = 0.5*diag(np(:,1))*(fxl+fxr) + 0.5*diag(np(:,2))*(fyl+fyr);

udi = 0.5*(up-um);
fn = fav + ca*udi;
