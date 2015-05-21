clear all; close all; clc;
porder = 1;
ngrid  = 20;

mesh   = mkmesh_square(ngrid,ngrid,porder);
master = mkmaster(mesh,2*porder);

kappa = 1;
c = 0*[10,10];
param = {kappa,c, @(c,n,h) ones(size(n,1), 1)};
source = @(p) 2*pi^2*sin(pi*p(:,1)).*sin(pi*p(:,2));
dbc    = @(p) 0*ones(size(p,1),1);

% HDG solver 
[uh,qh,uhath]=hdg_solve(master,mesh,source,dbc,param);
figure; scaplot(mesh, uh);

% HDG postprocessing 
mesh1   = mkmesh_square(ngrid,ngrid,porder+1);
master1 = mkmaster(mesh1,2*(porder+1));
[ustarh]=hdg_postprocess(master,mesh,master1,mesh1,uh,qh);
figure; scaplot(mesh1, ustarh);
