function app = mkapp
%MKAPP Create application structure template for the Euler equations.
%   APP=MKAPP
%
%      APP:   Application Structure
% 
app.nc = 4;
app.pg  = false;
app.arg = {};

app.bcm = [];  
app.bcs = [];

app.finvi = @euleri_roe;
app.finvb = @eulerb;
app.finvv = @eulerv;
app.src = [];
