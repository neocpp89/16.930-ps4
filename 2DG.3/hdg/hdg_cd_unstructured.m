clear all; close all; clc;
% common parameters
c = [1,1];
source = @(p) ones(size(p,1), 1);

% XXX: should we plot the acutal solution? 
plot_figures = 1;
p = [1, 3];
n = [0:2]; % number of refinements
kappas = [1, 1e-1, 1e-2];

np = numel(p);
nn = numel(n);
nk = numel(kappas);

[m,mplus] = mkmesh_mycircle(0.06, min(p), max(p)-min(p)+1);
plist = min(p):(max(p)+1);

meshseries = cell(np, nn);
mesh1series = cell(np, nn);

for i=1:np
    mi = find(plist == p(i));
    meshseries{i, 1} = mplus{mi};
    mesh1series{i, 1} = mplus{mi+1};
    for j=2:nn
        meshseries{i, j} = refinemesh(mplus{mi}, n(j));
        mesh1series{i, j} = refinemesh(mplus{mi+1}, n(j));
    end
end

for i=1:nn
    h = figure;
    set(h, 'units', 'inches', 'position', [1 1 4 4])
    set(h, 'PaperUnits','centimeters');
    set(h, 'Units','centimeters');
    pos=get(h,'Position');
    set(h, 'PaperSize', [pos(3) pos(4)]);
    set(h, 'PaperPositionMode', 'manual');
    set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
    meshplot(meshseries{1,i});
    title(sprintf('Mesh at refinement level %d', n(i)));
    print(sprintf('../../report/um_%d.png', i), '-dpng');
end

errs = zeros(np, nn, nk);
errsqx = zeros(np, nn, nk);
errsqy = zeros(np, nn, nk);
errspp = zeros(np, nn, nk);
for i=1:np
    for j=1:nn
        for k=1:nk
            porder = p(i);
            ngrid  = n(j);

            mesh   = meshseries{i, j};
            master = mkmaster(mesh,2*porder);

            kappa = kappas(k);
            taufn = @(c,n,h) abs(n*c') + kappa / 0.2;
            % taufn = @(c,n,h) 0.5*(abs(n*c') + n*c') + kappa / 0.2;
            param = {kappa,c, taufn};

            % This is actually a lie, it only does homogeneous solutions
            % (dirichlet does not quite work right now).
            dbc    = @(p) 0*ones(size(p,1),1);

            % HDG solver 
            [uh,qh,uhath]=hdg_solve(master,mesh,source,dbc,param);

            % HDG postprocessing 
            mesh1   = mesh1series{i, j};
            master1 = mkmaster(mesh1,2*(porder+1));
            [ustarh]=hdg_postprocess(master,mesh,master1,mesh1,uh,qh);

            if (plot_figures == 1)
                h = figure;
                set(h, 'units', 'inches', 'position', [1 1 4 4])
                set(h, 'PaperUnits','centimeters');
                set(h, 'Units','centimeters');
                pos=get(h,'Position');
                set(h, 'PaperSize', [pos(3) pos(4)]);
                set(h, 'PaperPositionMode', 'manual');
                set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
                scaplot(mesh, uh);
                title(sprintf('Solution u\nrefinement=%d, p=%d, \\kappa=%g', n(j), p(i), kappa));
                print(sprintf('../../report/umu_%d%d%d.png', i, j, k), '-dpng');
                h = figure;
                set(h, 'units', 'inches', 'position', [1 1 4 4])
                set(h, 'PaperUnits','centimeters');
                set(h, 'Units','centimeters');
                pos=get(h,'Position');
                set(h, 'PaperSize', [pos(3) pos(4)]);
                set(h, 'PaperPositionMode', 'manual');
                set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
                scaplot(mesh1, ustarh);
                title(sprintf('Solution u^*\nrefinement=%d, p=%d, \\kappa=%g', n(j), p(i), kappa));
                print(sprintf('../../report/umustar_%d%d%d.png', i, j, k), '-dpng');
            end
        end
    end
end

%{
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
%}
