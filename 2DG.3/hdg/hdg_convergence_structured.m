clear all; close all; clc;
% common parameters
kappa = 1;
c = 0*[1,1];
source = @(p) 2*pi^2*sin(pi*p(:,1)).*sin(pi*p(:,2));

% XXX: should we plot the acutal solution? 
plot_figures = 0;
p = [1, 2, 3];
n = [9, 17, 33, 41]; % number of elements per row is n-1
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

errs = zeros(np, nn, nt);
errsqx = zeros(np, nn, nt);
errsqy = zeros(np, nn, nt);
errspp = zeros(np, nn, nt);
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

            errs(i,j,k) = l2err(mesh, master, uh, exact_u);
            errsqx(i,j,k) = l2err(mesh, master, qh(:,1,:), exact_qx);
            errsqy(i,j,k) = l2err(mesh, master, qh(:,2,:), exact_qy);

            % HDG postprocessing 
            mesh1   = mkmesh_square(ngrid,ngrid,porder+1);
            master1 = mkmaster(mesh1,2*(porder+1));
            [ustarh]=hdg_postprocess(master,mesh,master1,mesh1,uh,qh);

            errspp(i,j,k) = l2err(mesh1, master1, ustarh, exact_u);

            if (plot_figures == 1)
                figure; scaplot(mesh, uh);
                title(sprintf('u N=%d, p=%d, fn=%d', n(j), p(i), k));
                figure; scaplot(mesh1, ustarh);
                title(sprintf('u^* N=%d, p=%d, fn=%d', n(j), p(i), k));
            end
        end
    end
end

% u
fprintf('u\n');
for k=1:nt
    h = figure;
    set(h, 'units', 'inches', 'position', [1 1 4 4])
    set(h, 'PaperUnits','centimeters');
    set(h, 'Units','centimeters');
    pos=get(h,'Position');
    set(h, 'PaperSize', [pos(3) pos(4)]);
    set(h, 'PaperPositionMode', 'manual');
    set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
    hold all;
    for i=1:np
        y = squeeze(errs(i,:,k))';
        loglog(n, y, 'DisplayName', sprintf('p = %d', p(i)));
        pf = polyfit(log(n(end-1:end))', log(y(end-1:end)), 1);
        fprintf('%s - Order %d has rate %g.\n', taulabels{k}, p(i), pf(1));
    end
    hold off;
    title(sprintf('Convergence of u\n(Structured, \\tau = %s)', taulabels{k}));
    xlabel('N');
    ylabel('L_2 error');
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    axis equal square;
    grid on;
    legend(gca, 'show', 'location', 'southwest');
    print(sprintf('../../report/cs_%d.png', k), '-dpng');
end

% qx
fprintf('qx\n');
for k=1:nt
    h = figure;
    set(h, 'units', 'inches', 'position', [1 1 4 4])
    set(h, 'PaperUnits','centimeters');
    set(h, 'Units','centimeters');
    pos=get(h,'Position');
    set(h, 'PaperSize', [pos(3) pos(4)]);
    set(h, 'PaperPositionMode', 'manual');
    set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
    hold all;
    for i=1:np
        y = squeeze(errsqx(i,:,k))';
        loglog(n, y, 'DisplayName', sprintf('p = %d', p(i)));
        pf = polyfit(log(n(end-1:end))', log(y(end-1:end)), 1);
        fprintf('%s - Order %d has rate %g.\n', taulabels{k}, p(i), pf(1));
    end
    hold off;
    title(sprintf('Convergence of q_x\n(Structured, \\tau = %s)', taulabels{k}));
    xlabel('N');
    ylabel('L_2 error');
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    axis equal square;
    grid on;
    legend(gca, 'show', 'location', 'southwest');
    print(sprintf('../../report/csqx_%d.png', k), '-dpng');
end

% qy
fprintf('qy\n');
for k=1:nt
    h = figure;
    set(h, 'units', 'inches', 'position', [1 1 4 4])
    set(h, 'PaperUnits','centimeters');
    set(h, 'Units','centimeters');
    pos=get(h,'Position');
    set(h, 'PaperSize', [pos(3) pos(4)]);
    set(h, 'PaperPositionMode', 'manual');
    set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
    hold all;
    for i=1:np
        y = squeeze(errsqy(i,:,k))';
        loglog(n, y, 'DisplayName', sprintf('p = %d', p(i)));
        pf = polyfit(log(n(end-1:end))', log(y(end-1:end)), 1);
        fprintf('%s - Order %d has rate %g.\n', taulabels{k}, p(i), pf(1));
    end
    hold off;
    title(sprintf('Convergence of q_y\n(Structured, \\tau = %s)', taulabels{k}));
    xlabel('N');
    ylabel('L_2 error');
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    axis equal square;
    grid on;
    legend(gca, 'show', 'location', 'southwest');
    print(sprintf('../../report/csqy_%d.png', k), '-dpng');
end

% u*
fprintf('u*\n');
for k=1:nt
    h = figure;
    set(h, 'units', 'inches', 'position', [1 1 4 4])
    set(h, 'PaperUnits','centimeters');
    set(h, 'Units','centimeters');
    pos=get(h,'Position');
    set(h, 'PaperSize', [pos(3) pos(4)]);
    set(h, 'PaperPositionMode', 'manual');
    set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
    hold all;
    for i=1:np
        y = squeeze(errspp(i,:,k))';
        loglog(n, y, 'DisplayName', sprintf('p = %d', p(i)));
        pf = polyfit(log(n(end-1:end))', log(y(end-1:end)), 1);
        fprintf('%s - Order %d has rate %g.\n', taulabels{k}, p(i), pf(1));
    end
    hold off;
    title(sprintf('Convergence of u^*\n(Structured, \\tau = %s)', taulabels{k}));
    xlabel('N');
    ylabel('L_2 error');
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    axis equal square;
    grid on;
    legend(gca, 'show', 'location', 'southwest');
    print(sprintf('../../report/cspp_%d.png', k), '-dpng');
end
