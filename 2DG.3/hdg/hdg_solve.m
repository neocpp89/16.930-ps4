function [uh,qh,uhath]=hdg_solve(master,mesh,source,dbc,param)

%HDG_SOLVE solves the convection-diffusion equation using the HDG method.
%   [uh,qh,uhath]=hdg_solve(mesh,master,source,dbc,param)
%
%      MASTER:       Master structure
%      MESH:         Mesh structure
%      SOURCE:       Source term
%      DBC:          Dirichlet data 
%      PARAM:        PARAM(1)   = diffusivity coefficient
%                    PARAM(2:3) = convective velocity
%      UH:           Approximate scalar variable
%      QH:           Approximate flux
%      UHATH:        Approximate trace                              


