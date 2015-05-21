function [uh,qh,uhath]=hdg_solve(master,mesh,source,dbc,param)

%HDG_SOLVE solves the convection-diffusion equation using the HDG method.
%   [uh,qh,uhath]=hdg_solve(mesh,master,source,dbc,param)
%
%      MASTER:       Master structure
%      MESH:         Mesh structure
%      SOURCE:       Source term
%      DBC:          Dirichlet data 
%      PARAM:        PARAM{1} = diffusivity coefficient
%                    PARAM{2} = convective velocity [cx, cy]
%                    PARAM{3} = stabilization parameter functionn tau(c,n,h)
%      UH:           Approximate scalar variable
%      QH:           Approximate flux
%      UHATH:        Approximate trace                              

% equation parameters
kappa = param{1};
c = param{2};
if(numel(param) >= 3)
    tau = param{3};
else
    % default
    tau = @(c,n,h) ones(size(n, 1), 1);
end

phi1d(:,:) = master.sh1d(:,1,:);
dphi1d(:,:) = master.sh1d(:,2,:);
phi(:,:) = master.shap(:,1,:);
dphidxi(:,:) = master.shap(:,2,:);
dphideta(:,:) = master.shap(:,3,:);

% components of jacobian
xxi = dphidxi'*squeeze(mesh.dgnodes(:,1,:));
xet = dphideta'*squeeze(mesh.dgnodes(:,1,:));
yxi = dphidxi'*squeeze(mesh.dgnodes(:,2,:));
yet = dphideta'*squeeze(mesh.dgnodes(:,2,:));
detJ = (xxi.*yet - xet.*yxi);

nel = size(mesh.t, 1);
nf = size(mesh.f, 1);
ndof = nf * numel(master.ploc1d);
ns = size(master.mass, 1);
ns1d = size(master.ploc1d, 1);

Ael = cell(nel, 1);
Bel = cell(nel, 1);
Cel = cell(nel, 1);
Del = cell(nel, 1);
Eel = cell(nel, 1);
Fel = cell(nel, 1);
Gel = cell(nel, 1);
Mel = cell(nel, 1);
QUel = cell(nel, 1);
QU0l = cell(nel, 1);

Hel = cell(nel, 1);
Rel = cell(nel, 1);

% build global matrix for trace values (global solve)
for i=1:nel
    % local matrices
    A = zeros(2*ns, 2*ns);
    B = zeros(2*ns, ns);
    C = zeros(2*ns, 3*ns1d);
    D = zeros(ns, ns);
    E = zeros(ns, 3*ns1d);
    M = zeros(3*ns1d, 3*ns1d);
    F = zeros(ns, 1);
    G = zeros(3*ns1d, 1);
    R = zeros(2*ns, 1);

    % components of inverse jacobian
    xix = diag(yet(:,i) ./ detJ(:,i));
    etax = diag(-yxi(:,i) ./ detJ(:,i));
    xiy = diag(-xet(:,i) ./ detJ(:,i));
    etay = diag(xxi(:,i) ./ detJ(:,i));

    % derivatives with respect to global coordinates
    dphidx = (dphidxi*xix + dphideta*etax);
    dphidy = (dphidxi*xiy + dphideta*etay);

    % mass matrix if non-constant jacobians
    svol = master.gwgh .* detJ(:, i);
    mm = phi*diag(svol)*phi';

    elnn = 1:ns;

    % 'A' matrix, (kinv*q_h, v)
    A1d = mm / kappa;
    A(elnn, elnn) = A1d; % xx
    A(ns+elnn, ns+elnn) = A1d; % yy

    % 'B' matrix, (u_h, div(v))
    Bxx = dphidx*diag(svol)*phi';
    Byy = dphidy*diag(svol)*phi';
    B(elnn, elnn) = Bxx;
    B(ns+elnn, elnn) = Byy;

    % 'D' matrix (interior part), ((c.div)*u_h, w)
    Dxx = c(1)*phi*diag(svol)*dphidx';
    Dyy = c(2)*phi*diag(svol)*dphidy';
    D(elnn, elnn) = Dxx + Dyy;

    % 'F' vector, (f, w)
    xg = phi'*mesh.dgnodes(:,:,i);
    fg = source(xg);
    F(elnn, 1) = phi*(svol.*fg);

    % have to go over edges for C, D, E, M matrices
    for j=1:3
        fidx = abs(mesh.t2f(i, j));
        edgenn = master.perm(:, j, 1); % only want ccw direction 
        tracenn = (j-1)*ns1d + (1:ns1d)';
        ep = mesh.dgnodes(edgenn, :, i);

        % tangent vector
        xsg = dphi1d'*ep(:,1);
        ysg = dphi1d'*ep(:,2);
        dsg = sqrt(ysg.*ysg + xsg.*xsg);

        % outward normal
        nepg = [ysg ./ dsg, -xsg ./ dsg];
        scale = (master.gw1d .* dsg);
        S = diag(scale);


        % 'C' matrix, <u_hat, v.n>
        Cx = phi1d*S*diag(nepg(:, 1))*phi1d';
        Cy = phi1d*S*diag(nepg(:, 2))*phi1d';
        C(edgenn, tracenn) = C(edgenn, tracenn) + Cx;
        C(ns + edgenn, tracenn) = C(ns + edgenn, tracenn) + Cy;

        % 'D' matrix, <(tau - c.n)*u_h, w>
        Dx = -phi1d*S*diag(c(1)*nepg(:, 1))*phi1d';
        Dy = -phi1d*S*diag(c(2)*nepg(:, 2))*phi1d';
        Dtau = phi1d*S*diag(tau(c, nepg, 1/nel))*phi1d';
        T = phi1d*diag(tau(c, nepg, 1/nel) .* scale)*phi1d';
        D(edgenn, edgenn) = D(edgenn, edgenn) + (Dtau+Dx+Dy);

        % 'E' matrix, <(tau - c.n)*u_hat, w>
        E(edgenn, tracenn) = E(edgenn, tracenn) + (Dtau+Dx+Dy);

        % 'M' matrix, <tau*u_hat, mu>
        M(tracenn, tracenn) = M(tracenn, tracenn) + Dtau;

        if (mesh.f(fidx, 4) < 0)
            % boundary
            gd = dbc(xg(edgenn));

            % should actually use a param for this
            gn = zeros(size(gd));
            G(tracenn) = G(tracenn) + phi1d*S*(gn + gd);
        end
    end


    % H and R (for global system for u_hat)
    P = [A B; -B' D];
    QU = P \ [C; E];
    QU0 = P \ [R; F];
    H = M + [C' -E']*QU;
    R = G - [C' -E']*QU0;

    % save the matrices in the cell array
    Ael{i} = A;
    Bel{i} = B;
    Cel{i} = C;
    Del{i} = D;
    Eel{i} = E;
    Fel{i} = F;
    Gel{i} = G;
    Mel{i} = M;
    QUel{i} = QU;
    QU0el{i} = QU0;

    Hel{i} = H;
    Rel{i} = R;
end

% generate global node numbering for u_hat
nn = zeros(3*ns1d,nel);
for i=1:nel
    for j=1:3
        fidx = mesh.t2f(i,j);
        tracenn = (j-1)*ns1d + (1:ns1d)';
        if (fidx > 0)
            nn(tracenn,i) = (fidx-1)*ns1d + (1:ns1d)';
        else 
            nn(tracenn,i) = (-fidx-1)*ns1d + (ns1d:-1:1)';
        end
    end
end

% assemble global matrix for u_hat
Hglob = sparse(3*ns1d*nel, 3*ns1d*nel);
Rglob = zeros(3*ns1d*nel, 1);
for i=1:nel
    nnel = nn(:, i);
    Hglob(nnel, nnel) = Hglob(nnel, nnel) + Hel{i};
    Rglob(nnel, 1) = Rglob(nnel, 1) + Rel{i};
end

% remove BC dofs and solve global matrix
uhath = dbc(1)*ones(size(Rglob, 1), 1);
rmidx = [];
keepidx = [];
for i=1:nel
    for j=1:3
        fidx = abs(mesh.t2f(i,j));
        tracenn = (j-1)*ns1d + (1:ns1d)';
        if (mesh.f(fidx, 4) < 0)
            rmidx = [rmidx; (fidx-1)*ns1d + (1:ns1d)';];
        else
            keepidx = [keepidx; (fidx-1)*ns1d + (1:ns1d)';];
        end
    end
end
keepidx = unique(sort(keepidx));
Hhat = Hglob(keepidx, keepidx);
Rhat = Rglob(keepidx);
uhath(keepidx) = Hhat \ Rhat;

% local solve on each element using traces
uh = zeros(ns, 1, nel);
qh = zeros(ns, 2, nel);
for i=1:nel
    nnel = nn(:, i);
    QU = QUel{i};
    QU0 = QU0el{i};
    quvec = QU0 + QU*uhath(nnel, 1);
    qh(:, 1, i) = quvec((1:ns), 1);
    qh(:, 2, i) = quvec((1:ns) + ns, 1);
    uh(:, 1, i) = quvec((1:ns) + 2*ns, 1);
end
