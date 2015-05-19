function [ustarh]=hdg_postprocess(master,mesh,master1,mesh1,uh,qh)

%HDG_POSTPROCESS postprocesses the HDG solution to obtain a better solution.
%   [ustarh]=hdg_postprocess(mesh,master,uh,qh,uhath)
%
%      MASTER:       Master structure of porder
%      MESH:         Mesh structure of porder
%      MASTER1:      Master structure of porder+1
%      MESH1:        Mesh structure of porder+1
%      UH:           Approximate scalar variable
%      QH:           Approximate flux
%      USTARH:       Postprocessed scalar variable

