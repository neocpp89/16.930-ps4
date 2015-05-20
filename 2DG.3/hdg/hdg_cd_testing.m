porder = 2;
ngrid  = 9;

mesh   = mkmesh_square(ngrid,ngrid,porder);
master = mkmaster(mesh,2*porder);

kappa = 1;
c = [10,10];
param = {kappa,c};
source = @(p) 10*ones(size(p,1),1);
dbc    = @(p) 0*ones(size(p,1),1);


% HDG solver 
[uh,qh,uhath]=hdg_solve(master,mesh,source,dbc,param);
scaplot(mesh, uh);

% HDG postprocessing 
mesh1   = mkmesh_square(ngrid,ngrid,porder+1);
master1 = mkmaster(mesh1,2*(porder+1));
[ustarh]=hdg_postprocess(master,mesh,master1,mesh1,uh,qh);

