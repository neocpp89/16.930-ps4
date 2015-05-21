clear all; close all; clc;
% common parameters
kappa = 1;
c = 0*[1,1];
source = @(p) 2*pi^2*sin(pi*p(:,1)).*sin(pi*p(:,2));

plot_figures = 1;
% p = [1, 2, 3];
% n = [9, 17, 33]; % number of elements per row is n-1
p = 1;
n = [8];
taufns = {@(c,n,h) h*ones(size(n,1),1), ...
    @(c,n,h) ones(size(n,1),1), ...
    @(c,n,h) (1/h)*ones(size(n,1),1)};
taulabels = {'h', '1', '1/h'};

% I actually made a small error (implemented the equations in the notes)
% so q = kappa * grad(u) instead of q = grad(u).
% Luckily, kappa = 1 in this case...
exact_u = @(xy) sin(pi*xy(:,1)).*sin(pi*xy(:,2));
exact_qx = @(xy) pi*cos(pi*xy(:,1)).*sin(pi*xy(:,2));
exact_qy = @(xy) pi*sin(pi*xy(:,1)).*cos(pi*xy(:,2));

np = numel(p);
nn = numel(n);
nt = numel(taufns);

errs = zeros(np, nn, nt)
for i=1:np
    for j=1:nn
        for k=1:nt
            porder = p(i);
            ngrid  = n(j);

            mesh   = mkmesh_square(ngrid,ngrid,porder);
            master = mkmaster(mesh,2*porder);

            taufn = taufns{k};
            param = {kappa,c, taufn};

            % This is actually a lie, it only does homogeneous solutions
            % (dirichlet does not quite work right now).
            dbc    = @(p) 0*ones(size(p,1),1);

            % HDG solver 
            [uh,qh,uhath]=hdg_solve(master,mesh,source,dbc,param);

            errs(i,j,k) = l2err(mesh, master, uh, exact_u)

            % HDG postprocessing 
            mesh1   = mkmesh_square(ngrid,ngrid,porder+1);
            master1 = mkmaster(mesh1,2*(porder+1));
            [ustarh]=hdg_postprocess(master,mesh,master1,mesh1,uh,qh);

            if (plot_figures == 1)
                figure; scaplot(mesh, uh);
                title(sprintf('u N=%d, p=%d, fn=%d', n(j), p(i), k));
                figure; scaplot(mesh1, ustarh);
                title(sprintf('u^* N=%d, p=%d, fn=%d', n(j), p(i), k));
            end
        end
    end
end

for i=1:np
    for j=1:nn
        for k=1:nt
        end
    end
end
