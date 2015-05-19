function u = rk4(residexpl,master,mesh,app,u,time,dt,nstep)
%RK4 Time integrator using a 4 stage Runge-Kutta scheme.
%   U=RK4(RESIDEXPL,MASTER,MESH,APP,U,TIME,DT,NSTEP)
%
%      RESIDEXPL:    Pointer to residual evaluation function
%                    R=RESIDEXPL(MASTER,MESH,APP,U,TIME)
%      MASTER:       Master structure
%      MESH:         Mesh structure
%      APP:          Application structure
%      U(NPL,NC,NT): Vector of unknowns
%                    NPL = size(mesh.plocal,1)
%                    NC = app.nc (number of equations in system)
%                    NT = size(mesh.t,1)
%      TIME:         Time
%      DT:           Time step
%      NSTEP:        Number of steps to be performed
%      R(NPL,NC,NT): Residual vector (=dU/dt) (already divided by mass
%                    matrix)                             
%
% - Written by: J. Peraire
%
for i = 1:nstep    
    k1 = dt*residexpl( master, mesh, app, u       , time       );
    k2 = dt*residexpl( master, mesh, app, u+0.5*k1, time+0.5*dt);
    k3 = dt*residexpl( master, mesh, app, u+0.5*k2, time+0.5*dt);
    k4 = dt*residexpl( master, mesh, app, u+    k3, time+    dt);
    u = u + k1/6 + k2/3 + k3/3 + k4/6;
    time = time + dt;
end
