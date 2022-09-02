function main_Direct(h,epsilon)
%
%     Author: Jorge Rodriguez Gutierrez
%     Date: 2022-07-25
%     Script that calculates the solution to the problem the steady-state 
%     Fokker-Planck equation with the direct method for different noise
%     levels. It uses the method of lsqminnorm. To be run in HPC.
%     Inputs from PBS file:
%         h: grid spacing
%         epsilon: noise intensity
%
addpath('functions')
method = "Direct";
%
t0cpu = cputime;
t0tic = tic;
%
mu = 1;
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
Ntot = x{1}.N*x{2}.N;
h = x{1}.d*x{2}.d;
%
A = getAfp(x,f,epsilon);
%
% The problem is A*rho = 0, and we want to find the solution rho ~= 0, with
% min ||rho||^2
%
% decompose rho = v + y, then Ay = -Av.
%
v = ones(Ntot-4,1);
b = -A*v;
%
y = lsqminnorm(A,b);
%
rhoVec = v + y;
%
rho = zeros(x{2}.N,x{1}.N);
rho(1,2:x{1}.N-1) = (reshape(rhoVec(1:x{1}.N-2),x{1}.N-2,1))';
rho(2:x{2}.N-1,:) = (reshape(rhoVec(x{1}.N-1:end-(x{1}.N-2)),x{1}.N,x{2}.N-2))';
rho(x{2}.N,2:x{1}.N-1) = (reshape(rhoVec(end-(x{1}.N-3):end),x{1}.N-2,1))';
%
% absolute value
rho = abs(rho);
% normalize
rho = rho/trapz2D(x,rho);
residual = A*rhoVec;
results.r = zeros(x{2}.N,x{1}.N);
results.r(1,2:x{1}.N-1) = (reshape(residual(1:x{1}.N-2),x{1}.N-2,1))';
results.r(2:x{2}.N-1,:) = (reshape(residual(x{1}.N-1:end-(x{1}.N-2)),x{1}.N,x{2}.N-2))';
results.r(x{2}.N,2:x{1}.N-1) = (reshape(residual(end-(x{1}.N-3):end),x{1}.N-2,1))';
results.rnorm = normGrid(results.r,h);
results.rrel = results.rnorm;
results.phiAvg = trapz2D(x,rho.*phi);
%
results.timecpu = cputime - t0cpu;
results.timetoc = toc(t0tic);
name = "rho" + method;
save(name,'mu','epsilon','x','f','phi','rho','results');
%
end