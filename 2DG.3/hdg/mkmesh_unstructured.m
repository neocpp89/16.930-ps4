function [mesh, meshplus] = mkmesh_unstructured(siz,porder,padditional)
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
if nargin<3, padditional=3; end

fd=@(p) drectangle(p, 0, 1, 0, 1);
% fh=@(p) min(0.01 + 0.3*abs(dcircle(p, 1, 1, 0)), 0.15) ;
fh = @(p) min(0.05 + 0.3*abs(dpoly(p,[0.5, 1.0; 1.0, 1.0; 1.0, 0.5])), 0.25);
[mesh.p,mesh.t] = distmesh2d(fd,fh,0.05,[0,0;1,1],[0,0;1,0;0,1;1,1]);
[mesh.p,mesh.t] = fixmesh(mesh.p,mesh.t);
[mesh.f,mesh.t2f] = mkt2f(mesh.t);

bndexpr = {'all(p(:,2)<1e-3)','all(p(:,1)>1-1e-3)', ...
           'all(p(:,2)>1-1e-3)','all(p(:,1)<1e-3)'};   
mesh.f = setbndnbrs(mesh.p,mesh.f,bndexpr);

mesh.fcurved = (mesh.f(:,4)<0);
ic = find(mesh.fcurved);
mesh.tcurved = repmat(false, size(mesh.t,1),1);
mesh.tcurved(mesh.f(ic,3)) = true;

mesh.porder = porder;
[mesh.plocal,mesh.tlocal] = uniformlocalpnts(mesh.porder);
mesh.dgnodes = createnodes(mesh,fd);

meshplus = cell(1+padditional, 1);
for i=1:padditional
    mesh1 = mesh;
    mesh1.porder = porder+i;
    [mesh1.plocal,mesh1.tlocal] = uniformlocalpnts(mesh1.porder);
    mesh1.dgnodes = createnodes(mesh1,fd);
    meshplus{i+1} = mesh1;
end
meshplus{1} = mesh;
