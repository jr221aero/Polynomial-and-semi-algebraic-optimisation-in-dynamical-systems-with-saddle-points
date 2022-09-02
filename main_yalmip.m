function main_yalmip(degV,epsilon,boundFlag)
%
%     Author: Jorge Rodriguez Gutierrez
%     Date: 2022-07-31
%     Script that calculates lower or upper bound for a deterministic
%     system given by f with the addition of noise using YALMIP.
%     Inputs from PBS file:
%         degV: degree of the polynomial auxiliary function
%         epsilon: noise intensity
%         boundFlag: "u" for upper bound, or "l" for lower bound
%
addpath('functions')
%
dimX = 2;
mu = 1;
%
% Generate the sdp variables, auxiliary function and solver Options
opsSolver = generateOpsSolver(0,0);
%    
[x,c,V,coe,delV,del2V] = generateInputs(dimX,degV);
%
phi = sum(x.*x); % Define the quantity of interest
phifun = sdp2fun(phi,coe);
%    
f = [-mu*x(1); x(2)-x(2)^3-mu*x(1)^2]; % Define RHS of governing equation
ffun = cell(2,1);
ffun{1} = sdp2fun(f(1),coe);
ffun{2} = sdp2fun(f(2),coe);
%
fdelV= delV*f; % Compute the dot product between f and del(V)
%
% Find the bound of the system
epsdel2V = epsilon*del2V; % Compute product of epsilon and del^2(V)
[D,coe,c,genCoe,bound,res] = findbound(c,fdelV,phi,boundFlag,opsSolver,coe,epsdel2V); % Find the bound
%
Dfun = sdp2fun(D,genCoe);
Vfun = sdp2fun(V,coe);
fdelVfun = sdp2fun(fdelV,coe);
epsdel2Vfun = sdp2fun(epsdel2V,coe);
%                            
name = "yalmip" + string(boundFlag);
save(name,'mu','epsilon','ffun','phifun','Dfun','Vfun','fdelVfun','epsdel2Vfun','boundFlag','bound','degV','res')
%
end