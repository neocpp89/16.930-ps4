function [uh,qh,uhath]=hdg_solve(master,mesh,source,dbc,param)

%HDG_SOLVE solves the convection-diffusion equation using the HDG method.
%   [uh,qh,uhath]=hdg_solve(mesh,master,source,dbc,param)
%
%      MASTER:       Master structure
%      MESH:         Mesh structure
%      SOURCE:       Source term
%      DBC:          Dirichlet data 
%      PARAM:        PARAM(1)   = diffusivity coefficient
%                    PARAM(2:3) = convective velocity
%      UH:           Approximate scalar variable
%      QH:           Approximate flux
%      UHATH:        Approximate trace                              

% equation parameters
kappa = param{1};
c = param{2};

% map node to dof index (position in solution vector)

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

% build global matrix for trace values (global solve)
for i=1:nel
    ns = size(master.mass, 1);
    ns1d = size(master.ploc1d, 1);

    % local matrices
    A = sparse(2*ns, 2*ns);
    B = sparse(2*ns, ns);
    C = sparse(2*ns, 3*ns1d);
    D = sparse(ns, ns);
    E = sparse(ns, 3*ns1d);
    M = sparse(3*ns1d, 3*ns1d);
    F = zeros(3*ns1d, 1);

    % components of inverse jacobian
    xix = diag(yet(:,i) ./ detJ(:,i));
    etax = diag(-yxi(:,i) ./ detJ(:,i));
    xiy = diag(-xet(:,i) ./ detJ(:,i));
    etay = diag(xxi(:,i) ./ detJ(:,i));

    % derivatives with respect to global coordinates
    dphidx = (dphidxi*xix + dphideta*etax);
    dphidy = (dphidxi*xiy + dphideta*etay);

    % mass matrix if non-constant jacobians
    mm = phi*diag(master.gwgh .* detJ(:, i))*phi';

    elnn = 1:ns;

    % 'A' matrix, (kinv*q_h, v)
    A1d = mm / kappa;
    A(elnn, elnn) = A1d; % xx
    A(ns+elnn, ns+elnn) = A1d; % yy

    % 'B' matrix, (u_h, div(v))
    Bxx = dphidx*phi';
    Byy = dphidy*phi';
    B(elnn, elnn) = Bxx;
    B(ns+elnn, elnn) = Byy;

    % have to go over edges for C, D, E, M matrices

    BT = B';
    for j=1:3
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

    end
end

% schur complement and solve global matrix for u_hat


% local solve on each element using traces
for i=1:nel
end
