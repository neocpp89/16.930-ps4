function [mesh, mesh1] = mkmesh_mycircle(siz,porder)
%MKMESH_SQUARE Creates 2D mesh data structure for unit circle using DISTMESH2D.
%   MESH=MKMESH_SQUARE(SIZ,PORDER)
%
%      MESH:      Mesh structure
%      SIZ:       Desired element size 
%      PORDER:    Polynomial Order of Approximation (default=1)
%
%   See also: DISTMESH2D, FIXMESH (part of DISTMESH), MKT2F, SETBNDNBRS, 
%             UNIFORMLOCALPNTS, CREATENODES
%
if nargin>0 & siz == 0.0, siz=0.4; end
if nargin<1, siz=0.4; end
if nargin<2, porder=1; end

fd=@(p) sqrt(sum(p.^2,2))-1;
fh=@(p) (1.2 - 0.2*(p(:,1)+p(:,2)) - 4*dcircle(p,0,0,1));
[mesh.p,mesh.t] = distmesh2d(fd,fh,siz,[-1,-1;1,1],[]);
[mesh.p,mesh.t] = fixmesh(mesh.p,mesh.t);
[mesh.f,mesh.t2f] = mkt2f(mesh.t);

bndexpr = {'all(sqrt(sum(p.^2,2))>1-1e-3)'};     
mesh.f = setbndnbrs(mesh.p,mesh.f,bndexpr);

mesh.fcurved = (mesh.f(:,4)<0);
ic = find(mesh.fcurved);
mesh.tcurved = repmat(false, size(mesh.t,1),1);
mesh.tcurved(mesh.f(ic,3)) = true;

mesh.porder = porder;
[mesh.plocal,mesh.tlocal] = uniformlocalpnts(mesh.porder);
mesh.dgnodes = createnodes(mesh,fd);

mesh1 = mesh;
mesh1.porder = porder+1;
[mesh1.plocal,mesh1.tlocal] = uniformlocalpnts(mesh1.porder);
mesh1.dgnodes = createnodes(mesh1,fd);

