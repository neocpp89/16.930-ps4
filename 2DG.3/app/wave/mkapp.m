function app = mkapp
%MKAPP Create application structure template for the Wave equation.
%   APP=MKAPP
%
%      APP:   Application Structure
% 
app.nc = 3;
app.pg  = false;
app.arg = {};

app.bcm = [];
app.bcs = [];

app.finvi = @wavei_roe;
app.finvb = @waveb;
app.finvv = @wavev;
app.src = [];

