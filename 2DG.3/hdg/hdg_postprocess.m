function [ustarh]=hdg_postprocess(master,mesh,master1,mesh1,uh,qh)
%HDG_POSTPROCESS postprocesses the HDG solution to obtain a better solution.
%   [ustarh]=hdg_postprocess(mesh,master,uh,qh,uhath)
%
%      MASTER:       Master structure of porder
%      MESH:         Mesh structure of porder
%      MASTER1:      Master structure of porder+1
%      MESH1:        Mesh structure of porder+1
%      UH:           Approximate scalar variable
%      QH:           Approximate flux
%      USTARH:       Postprocessed scalar variable

phi1d(:,:) = master.sh1d(:,1,:);
dphi1d(:,:) = master.sh1d(:,2,:);
phi(:,:) = master.shap(:,1,:);
dphidxi(:,:) = master.shap(:,2,:);
dphideta(:,:) = master.shap(:,3,:);

phi1d1(:,:) = master1.sh1d(:,1,:);
dphi1d1(:,:) = master1.sh1d(:,2,:);
phi1(:,:) = master1.shap(:,1,:);
dphidxi1(:,:) = master1.shap(:,2,:);
dphideta1(:,:) = master1.shap(:,3,:);

% seems messy but I don't know how else to get this
sgw1 = shape2d(master.porder, master.plocal, master1.gpts);
phigw1 = squeeze(sgw1(:,1,:));

% components of jacobian
dgx = squeeze(mesh.dgnodes(:,1,:));
dgy = squeeze(mesh.dgnodes(:,2,:));
xxi = dphidxi'*dgx;
xet = dphideta'*dgx;
yxi = dphidxi'*dgy;
yet = dphideta'*dgy;
detJ = (xxi.*yet - xet.*yxi);

% on new grid
dgx1 = squeeze(mesh1.dgnodes(:,1,:));
dgy1 = squeeze(mesh1.dgnodes(:,2,:));
xxi1 = dphidxi1'*dgx1;
xet1 = dphideta1'*dgx1;
yxi1 = dphidxi1'*dgy1;
yet1 = dphideta1'*dgy1;
detJ1 = (xxi1.*yet1 - xet1.*yxi1);

nel = size(mesh.t, 1);
ns = size(master.mass, 1);
ns1d = size(master.ploc1d, 1);

ustarh = zeros(size(master1.mass, 1), 1, nel);
for i=1:nel
    % components of inverse jacobian
    xix = diag(yet(:,i) ./ detJ(:,i));
    etax = diag(-yxi(:,i) ./ detJ(:,i));
    xiy = diag(-xet(:,i) ./ detJ(:,i));
    etay = diag(xxi(:,i) ./ detJ(:,i));

    % derivatives with respect to global coordinates
    dphidx = (dphidxi*xix + dphideta*etax);
    dphidy = (dphidxi*xiy + dphideta*etay);

    % components of inverse jacobian (p+1 mesh)
    xix1 = diag(yet1(:,i) ./ detJ1(:,i));
    etax1 = diag(-yxi1(:,i) ./ detJ1(:,i));
    xiy1 = diag(-xet1(:,i) ./ detJ1(:,i));
    etay1 = diag(xxi1(:,i) ./ detJ1(:,i));

    % derivatives with respect to global coordinates (p+1 mesh)
    dphidx1 = (dphidxi1*xix1 + dphideta1*etax1);
    dphidy1 = (dphidxi1*xiy1 + dphideta1*etay1);

    S = diag(detJ(:, i) .* master.gwgh);
    S1 = diag(detJ1(:, i) .* master1.gwgh);

    qhx = qh(:, 1, i);
    qhy = qh(:, 2, i);
    u = uh(:, 1, i);

    A = dphidx1*S1*dphidx1' + dphidy1*S1*dphidy1';
    Bxx = dphidx1*S1*phigw1';
    Byy = dphidy1*S1*phigw1';
    F = (Bxx*qhx + Byy*qhy);

    % make sure sum is conserved
    ff = phi*S*ones(size(S,1), 1);
    ff1 = phi1*S1*ones(size(S1,1), 1);
    A(end, :) = ff1';
    F(end, :) = ff'*u;

    ustarel = A \ F;
    ustarh(:, 1, i) = ustarel;
end
