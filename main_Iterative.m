function main_Iterative(h,epsilon)
%
%     Author: Jorge Rodriguez Gutierrez
%     Date: 2022-07-28
%     Script that calculates the solution to the problem the steady-state 
%     Fokker-Planck equation the method of simple iterations.
%     To be runned in HPC.
%     Inputs from PBS file:
%         h: grid spacing
%         epsilon: noise intensity
%
addpath('functions')
method = "Iterative";
%
t0cpu = cputime;
t0tic = tic;
mu = 1;
%
parameters.maxiter = 3e6+1;
parameters.tol = 1e-6;
%
dimX = 2;
x = cell(dimX,1);
%
x{1}.max = 2;
x{1}.min = -x{1}.max;
x{1}.d = h;
x{2}.max = x{1}.max;
x{2}.min = -x{2}.max;
x{2}.d = x{1}.d;
%
x = getSpatialDomain(x,dimX);
f = getGovEqn(x,dimX,mu);
phi = getPhi(x,dimX);
%
rho0 = ones(x{2}.N,x{1}.N);
rho0(1,1) = 0;
rho0(x{2}.N,1) = 0;
rho0(1,x{1}.N) = 0;
rho0(x{2}.N,x{1}.N) = 0;
%
ix = 2:x{1}.N-1;
iy = 2:x{2}.N-1;
ir = x{1}.N-2:x{1}.N-1;
il = 2:3;
ib = x{2}.N-2:x{2}.N-1;
it = 2:3;
%
Are = -(2*epsilon/x{1}.d)./(f{1}.val(iy,x{1}.N)-1.5*epsilon/x{1}.d);
Ari = (0.5*epsilon/x{1}.d)./(f{1}.val(iy,x{1}.N)-1.5*epsilon/x{1}.d);
Ale = (2*epsilon/x{1}.d)./(f{1}.val(iy,1)+1.5*epsilon/x{1}.d);
Ali = -(0.5*epsilon/x{1}.d)./(f{1}.val(iy,1)+1.5*epsilon/x{1}.d);
Abe = -(2*epsilon/x{2}.d)./(f{2}.val(x{2}.N,ix)-1.5*epsilon/x{2}.d);
Abi = (0.5*epsilon/x{2}.d)./(f{2}.val(x{2}.N,ix)-1.5*epsilon/x{2}.d);
Ate = (2*epsilon/x{2}.d)./(f{2}.val(1,ix)+1.5*epsilon/x{2}.d);
Ati = -(0.5*epsilon/x{2}.d)./(f{2}.val(1,ix)+1.5*epsilon/x{2}.d);
%
Ar = [Ari Are];
Al = [Ale Ali];
Ab = [Abi;Abe];
At = [Ate;Ati];
%
% Right boundary
rho0(iy,x{1}.N) = sum(Ar.*rho0(iy,ir),2);
% Left boundary
rho0(iy,1) = sum(Al.*rho0(iy,il),2);
% Bottom boundary
rho0(x{2}.N,ix) = sum(Ab.*rho0(ib,ix));
% Top boundary
rho0(1,ix) = sum(At.*rho0(it,ix));
%
rho0 = rho0/trapz2D(x,rho0);
%
[rho,results] = simpleIterations(rho0,x,f,epsilon,parameters);
%
results.phiAvg = trapz2D(x,rho.*phi);
results.timecpu = cputime - t0cpu;
results.timetoc = toc(t0tic);
name = "rho" + method;
save(name,'mu','epsilon','x','f','phi','rho','parameters','results');
%
end