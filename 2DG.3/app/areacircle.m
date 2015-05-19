function [area,perim] = areacircle(siz,porder)
%AREACIRCLE Calcualte the area and perimeter of a unit circle
%   [AREA,PERIM]=AREACIRCLE(SP,PORDER)
%
%      SIZ:       Desired element size 
%      PORDER:    Polynomial Order of Approximation (default=1)
%      AREA:      Area of the circle (\pi)
%      PERIM:     Perimeter of the circumference (2*\pi)
%
mesh = mkmesh_circle(siz,porder);
master = mkmaster(mesh);

xxi = squeeze(master.shap(:,2,:))'*squeeze(mesh.dgnodes(:,1,:));
xet = squeeze(master.shap(:,3,:))'*squeeze(mesh.dgnodes(:,1,:));
yxi = squeeze(master.shap(:,2,:))'*squeeze(mesh.dgnodes(:,2,:));
yet = squeeze(master.shap(:,3,:))'*squeeze(mesh.dgnodes(:,2,:));

jac = xxi.*yet - xet.*yxi;
jcw = diag(master.gwgh)*jac;
area = sum(sum(jcw));

ii=find(mesh.f(:,4)<0);
dgnodes=zeros(master.porder+1,2,length(ii));
for k=1:length(ii)
    j=find(mesh.t2f(mesh.f(ii(k),3),:)==ii(k));
    dgnodes(:,:,k) = mesh.dgnodes(master.perm(:,j,1),:,mesh.f(ii(k),3));
end
xxi = squeeze(master.sh1d(:,2,:))'*squeeze(dgnodes(:,1,:));
yxi = squeeze(master.sh1d(:,2,:))'*squeeze(dgnodes(:,2,:));
jac = sqrt(xxi.^2 + yxi.^2);
jcw = diag(master.gw1d)*jac;
perim = sum(sum(jcw));