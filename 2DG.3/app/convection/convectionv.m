function [fx,fy] = convectionv(u,p,param,time)
%CONVECTIONV Calculate the volume flux for the Linear Convection equation.
%   [FX,FY]=CONVECTIONV(U,P,PARAM,TIME)
%
%      U(np):      np left (or plus) states
%      P(np,2):    np x,y coordinates
%      PARAM{1}:   Cell array containing either 
%                  - a constant velocity field [u,v] = PARAM{1}
%                  - a pointer to a funvtion that returns of velocity
%                    a velocity field as a function of P vvec = PARAM{1}(p)
%      TIME:       Not used
%      FX(np):     np fluxes in the x direction (f plus)  
%      FY(np):     np fluxes in the y direction (f plus)  
%                          
% - Written by: J. Peraire
%
if isa(param{1},'double')
    vel = param{1};
    fx = vel(1)*u;
    fy = vel(2)*u;
else
    vfield = param{1}(p);
    fx = vfield(:,1).*u;
    fy = vfield(:,2).*u;
end

