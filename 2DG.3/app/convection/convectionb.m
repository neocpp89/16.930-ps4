function fn = convectionb(up,np,ib,ui,p,param,time)
%CONVECTIONB Calculate the boundary flux for the Linear Convection equation.
%   FN=CONVECTIONB(UP,NP,IB,UI,P,PARAM,TIME)
%
%      UP(np):     np plus states
%      NP(np,2):   np normal plus vectors 
%      IB:         Boundary type
%                  - IB: 1 Far-Field (Radiation)
%      UI(1):      Infinity state associated with IB
%      P(np,2):    np x,y coordinates
%      PARAM{1}:   Cell array containing either the wave speed c=PARAM{1}
%      TIME:       Time
%      FN(np):     np normal fluxes (f plus)  
%                          
% - Written by: J. Peraire
% 
um = repmat(ui,size(up,1),1);
fn = convectioni(up,um,np,p,param,time);

