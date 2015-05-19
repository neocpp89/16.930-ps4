function [fx,fy] = wavev(u,p,param,time)
%WAVEV Calculate the volume flux for the Wave equation.
%   [FX,FY]=WAVEV(U,P,PARAM,TIME)
%
%      U(np,3):    np left (or plus) states
%      P:          Not used
%      PARAM{1}:   Cell array containing either the wave speed c=PARAM{1}
%      TIME:       Not used
%      FX(np,3):   np fluxes in the x direction (f plus)  
%      FY(np,3):   np fluxes in the y direction (f plus)  
%                          
% - Written by: J. Peraire
% 
c = param{1};
fx = -c*[u(:,3), zeros(size(u(:,1))), u(:,1)];
fy = -c*[zeros(size(u(:,1))), u(:,3), u(:,2)];
