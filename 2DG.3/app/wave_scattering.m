%DRIVER FOR THE WAVE SCATTERING PROBLEM                          
%
d0=fileparts([pwd,filesep]);
addpath([d0,'/wave']);     % Add path to application subdirectory

m      = 11;
n      = 20;
porder = 5;

time  = 100;
dt    = 0.6e-02;
nstep = 10;
ncycl = ceil(time/(nstep*dt));
c = 1;
k = [3,0];
km = sqrt(k(1)^2+k(2)^2);

mesh = mkmesh_trefftz(m,n,porder,[0,0,1]);
master = mkmaster(mesh,2*porder);
app = mkapp;
app.pg = true;

app.bcm = [2,3];            % Manually set boundary conditions
app.bcs = [0,0,0; 0,0,0; 0,0, 0];   

ub = @(c,k,pg,time)  sin(pg*k'-c*sqrt(k(1)^2+k(2)^2)*time);
app.arg = {c,k,ub};

u = initu(mesh,app,{0,0,0});

tm = 0.0;
for i = 1:ncycl
                            % Incident wave
    ue(:,3,:) = sin(mesh.dgnodes(:,1,:)*k(1)+mesh.dgnodes(:,2,:)*k(2)-c*km*tm);  
    ue(:,1,:) = -k(1)*ue(:,3,:)/km;
    ue(:,2,:) = -k(2)*ue(:,3,:)/km;
    if i == 1
        u = ue;
    end
    
    subplot(2,1,1), scaplot(mesh,u(:,3,:)-ue(:,3,:),[-1.2,1.2],[],1); axis off;
    subplot(2,1,2), scaplot(mesh,u(:,3,:),[-1.2,1.2],[],1); axis off;  
    time=sprintf('t* = %.3f',nstep*dt*i);
    title(time,'Color','white','FontSize',16,'FontName','Courier')

    u = rk4(@rinvexpl,master,mesh,app,u,tm,dt,nstep);
    tm = tm + nstep*dt;
end