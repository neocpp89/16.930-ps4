function app = mkapp
%MKAPP Create application structure template for the Linear Convection equation.
%   APP=MKAPP
%
%      APP:   Application Structure
% 
app.nc = 1;
app.pg  = true;
app.arg = {};

app.bcm = [];
app.bcs = [];

app.finvi = @convectioni;
app.finvb = @convectionb;
app.finvv = @convectionv;
app.src = [];


