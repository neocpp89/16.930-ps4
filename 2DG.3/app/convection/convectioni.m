function fn = convectioni(up,um,np,p,param,time)
%CONVECTIONI Calculate Interface Upwind Flux for the linear convection equation.
%   FN=CONVECTIONI(UP,UM,NP,P,PARAM,TIME)
%
%      UP(np):     np plus states
%      UR(np):     np minus states
%      NP(np,2):   np normal plus vectors 
%      P(np,2):    np x,y coordinates
%      PARAM{1}:   Cell array containing either 
%                  - a constant velocity field [u,v] = PARAM{1}
%                  - a pointer to a funvtion that returns of velocity
%                    a velocity field as a function of P vvec = PARAM{1}(p)
%      TIME:       Not used
%      FN(np):     np normal fluxes (f plus)    
%                          
% - Written by: J. Peraire
%

if isa(param{1},'double')
    vel = param{1};
    vn = np*vel';
else
    vfield = param{1}(p);
    vn = sum(np.*vfield,2);
end

avn = abs(vn);
fn = 0.5*vn.*(up+um)+0.5*avn.*(up-um);
