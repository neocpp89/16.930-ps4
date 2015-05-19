%DRIVER FOR LINEAR CONVECTION PROBLEM                          
%
d0=fileparts([pwd,filesep]);
addpath([d0,'/convection']);   % Add path to application subdirectory

m      = 8;                    
n      = 8;
porder = 3;

time  = 2*pi;
dt    = 2.e-02;
nstep = 5;
ncycl = ceil(time/(nstep*dt));

mesh = mkmesh_square(m,n,porder,1);
% mesh = mkmesh_distort(mesh,0.05);     % Uncomment for mesh distortion 
master = mkmaster(mesh,2*porder);
app = mkapp;

app.bcm = [1,1,1,1];                    % Manually set boundary conditions
app.bcs = [0];

vf = @(p) fliplr(p-0.5)*diag([-1,1]);   % Rotating field
app.arg = {vf};

init  = @(dg) exp(-120*((dg(:,1,:)-0.6).^2 + (dg(:,2,:)-0.5).^2));  % Gaussian hill
u = initu(mesh,app,{init});

tm = 0.0;
for i=1:ncycl   
    scaplot(mesh,u,[-0.1,1.2],[],1); axis off; 
    u = rk4(@rinvexpl,master,mesh,app,u,tm,dt,nstep);
    tm = tm + nstep*dt;
end

while tm < time                         % Complete an exact revolution
    u = rk4(@rinvexpl,master,mesh,app,u,tm,dt,1);
    tm = tm + dt;
end
scaplot(mesh,u,[-0.1,1.2],[],1); axis off;