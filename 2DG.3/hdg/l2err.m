function err = l2err(mesh, master, u, exact)
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
    err2 = zeros(nel, 1);
    for i=1:nel
        xyg = phi'*mesh.dgnodes(:, :, i);
        eg2 = (phi'*u(:,1,i) - exact(xyg)).^2; 
        err2(i) = (master.gwgh .* detJ(:, i))'*eg2;
    end

    err = sqrt(sum(err2));
end
