function fn = wavei(ul,ur,nl,p,param,time)

c = param{1};
ca = abs(c);

fxl = -c*[ul(:,3), zeros(size(ul(:,1))), ul(:,1)];
fyl = -c*[zeros(size(ul(:,1))), ul(:,3), ul(:,2)];
fxr = -c*[ur(:,3), zeros(size(ur(:,1))), ur(:,1)];
fyr = -c*[zeros(size(ur(:,1))), ur(:,3), ur(:,2)];
fav = 0.5*diag(nl(:,1))*(fxl+fxr) + 0.5*diag(nl(:,2))*(fyl+fyr);

% Roe Flux
qb = 0.5*ca*(ul(:,1)-ur(:,1)).*nl(:,1) + (ul(:,2)-ur(:,2)).*nl(:,2);
ub = 0.5*ca*(ul(:,3)-ur(:,3));
fn(:,1) = fav(:,1) + qb.*nl(:,1);
fn(:,2) = fav(:,2) + qb.*nl(:,2);
fn(:,3) = fav(:,3) + ub;

% Lax Friedrichs' Flux
% vel = param{1};
% udi = 0.5*(ul-ur);
% fn = fav + abs(vel)*udi;
