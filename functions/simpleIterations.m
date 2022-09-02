function [up,results] = simpleIterations(u0,x,f,epsilon,parameters)
%
%     Author: Jorge Rodriguez Gutierrez
%     Date: 2022-07-28
%     Function that calculates the solution of the steady-state
%     Fokker-Planck Equation by using the method of simple iterations (Gauss-Seidel)
%     Second order of accuracy
%     Inputs:
%         -u0: initial guess for the solution
%         -x: spatial coordinates
%         -f: function giving the dynamics of the system
%         -epsilon: noise intensity
%         -parameters: stopping criteria: maximum number of
%         iterations and tolerance
%     Outputs:
%         -up: solution
%         -results: information about the solution, containing the
%         iterations performed and the error metrics
%
h = x{1}.d*x{2}.d;
maxiter = parameters.maxiter;
tol = parameters.tol;
%
% Initialise solution array
um = u0; % u_{n-1}
dcheck = 1e2; % frequency with which the checks are performed
ncheck = floor(maxiter/dcheck)+1; % number of checks that will be performed
checks = zeros(1,ncheck); % contains the indices corresponding to the iterations at which the checks are performed
dnorm = zeros(1,ncheck); % error estimation
drel = zeros(1,ncheck); % relative error estimation
enorm = zeros(1,ncheck); % true error
erel = zeros(1,ncheck); % true relative error
rnorm = zeros(1,ncheck); % residual
rrel = zeros(1,ncheck); % relative residual
q = zeros(1,ncheck); % true error convergence rate
%
% Spatial derivative: central difference second order
ix = 2:x{1}.N-1;
ixp1 = 3:x{1}.N;
ixm1 = 1:x{1}.N-2;
iy = 2:x{2}.N-1;
iyp1 = 3:x{2}.N;
iym1 = 1:x{2}.N-2;
%
ir = x{1}.N-2:x{1}.N-1;
il = 2:3;
ib = x{2}.N-2:x{2}.N-1;
it = 2:3;
%
Akp1 = zeros(x{2}.N,x{1}.N);
Akm1 = zeros(x{2}.N,x{1}.N);
AkpN = zeros(x{2}.N,x{1}.N);
AkmN = zeros(x{2}.N,x{1}.N);
%
Akp1(iy,ix) = 0.5*(epsilon/(x{1}.d^2)-0.5*f{1}.val(iy,ixp1)/x{1}.d)/(epsilon*(1/(x{1}.d^2)+1/(x{2}.d^2)));
Akm1(iy,ix) = 0.5*(epsilon/(x{1}.d^2)+0.5*f{1}.val(iy,ixm1)/x{1}.d)/(epsilon*(1/(x{1}.d^2)+1/(x{2}.d^2)));
AkpN(iy,ix) = 0.5*(epsilon/(x{2}.d^2)-0.5*f{2}.val(iyp1,ix)/x{2}.d)/(epsilon*(1/(x{1}.d^2)+1/(x{2}.d^2)));
AkmN(iy,ix) = 0.5*(epsilon/(x{2}.d^2)+0.5*f{2}.val(iym1,ix)/x{2}.d)/(epsilon*(1/(x{1}.d^2)+1/(x{2}.d^2)));
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
u = u0; % u_{n}
% Interior points
for j = 2:x{2}.N-1
    u(j,ix) = Akp1(j,ix).*u(j,ixp1) + AkpN(j,ix).*u(j+1,ix) + AkmN(j,ix).*u(j-1,ix);
    for i = 2:x{1}.N-1
        u(j,i) = u(j,i) + Akm1(j,i)*u(j,i-1);
    end
end
% Right boundary
u(iy,x{1}.N) = sum(Ar.*u(iy,ir),2);
% Left boundary
u(iy,1) = sum(Al.*u(iy,il),2);
% Bottom boundary
u(x{2}.N,ix) = sum(Ab.*u(ib,ix));
% Top boundary
u(1,ix) = sum(At.*u(it,ix));
%
u = u/trapz2D(x,u);
%
up = u; % u_{n+1}
% Interior points
for j = 2:x{2}.N-1
    up(j,ix) = Akp1(j,ix).*up(j,ixp1) + AkpN(j,ix).*up(j+1,ix) + AkmN(j,ix).*up(j-1,ix);
    for i = 2:x{1}.N-1
        up(j,i) = up(j,i) + Akm1(j,i)*up(j,i-1);
    end
end
% Right boundary
up(iy,x{1}.N) = sum(Ar.*up(iy,ir),2);
% Left boundary
up(iy,1) = sum(Al.*up(iy,il),2);
% Bottom boundary
up(x{2}.N,ix) = sum(Ab.*up(ib,ix));
% Top boundary
up(1,ix) = sum(At.*up(it,ix));
%
up = up/trapz2D(x,up);
%
k = 1;
ic = 1;
checks(ic) = k;
invunorm = 1/normGrid(u,h);
d = u-um;
dnorm(ic) = normGrid(d,h);
drel(ic) = dnorm(ic)*invunorm;
q(ic) = normGrid(up-u,h)/dnorm(ic);
e = (q(ic)/(q(ic)-1))*d;
enorm(ic) = normGrid(e,h);
erel(ic) = enorm(ic)*invunorm;
r = u(iy,ix) - Akp1(iy,ix).*u(iy,ixp1) - Akm1(iy,ix).*u(iy,ixm1) - AkpN(iy,ix).*u(iyp1,ix) - AkmN(iy,ix).*u(iym1,ix);
rnorm(ic) = normGrid(r,h);
rrel(ic) = rnorm(ic)*invunorm;
%        
for ic = 2:ncheck
    for m = 2:dcheck-1
        k = k + 1;
        % Interior points
        for j = 2:x{2}.N-1
            up(j,ix) = Akp1(j,ix).*up(j,ixp1) + AkpN(j,ix).*up(j+1,ix) + AkmN(j,ix).*up(j-1,ix);
            for i = 2:x{1}.N-1
                up(j,i) = up(j,i) + Akm1(j,i)*up(j,i-1);
            end
        end
        % Right boundary
        up(iy,x{1}.N) = sum(Ar.*up(iy,ir),2);
        % Left boundary
        up(iy,1) = sum(Al.*up(iy,il),2);
        % Bottom boundary
        up(x{2}.N,ix) = sum(Ab.*up(ib,ix));
        % Top boundary
        up(1,ix) = sum(At.*up(it,ix));
        %
        up = up/trapz2D(x,up);
    end
%    
    for t = 1:2
        k = k + 1;
        um = u;
        u = up;
        % Interior points
        for j = 2:x{2}.N-1
            up(j,ix) = Akp1(j,ix).*up(j,ixp1) + AkpN(j,ix).*up(j+1,ix) + AkmN(j,ix).*up(j-1,ix);
            for i = 2:x{1}.N-1
                up(j,i) = up(j,i) + Akm1(j,i)*up(j,i-1);
            end
        end
        % Right bo% Right boundary
        up(iy,x{1}.N) = sum(Ar.*up(iy,ir),2);
        % Left boundary
        up(iy,1) = sum(Al.*up(iy,il),2);
        % Bottom boundary
        up(x{2}.N,ix) = sum(Ab.*up(ib,ix));
        % Top boundary
        up(1,ix) = sum(At.*up(it,ix));
        %
        up = up/trapz2D(x,up);
    end
    checks(ic) = k;
    invunorm = 1/normGrid(u,h);
    d = u-um;
    dnorm(ic) = normGrid(d,h);
    drel(ic) = dnorm(ic)*invunorm;
    q(ic) = normGrid(up-u,h)/dnorm(ic);
    e = (q(ic)/(q(ic)-1))*d;
    enorm(ic) = normGrid(e,h);
    erel(ic) = enorm(ic)*invunorm;
    r = u(iy,ix) - Akp1(iy,ix).*u(iy,ixp1) - Akm1(iy,ix).*u(iy,ixm1) - AkpN(iy,ix).*u(iyp1,ix) - AkmN(iy,ix).*u(iym1,ix);
    rnorm(ic) = normGrid(r,h);
    rrel(ic) = rnorm(ic)*invunorm;
%
    if erel(ic) < tol
        if ic < ncheck
            dnorm(ic+1:end) = NaN;
            drel(ic+1:end) = NaN;
            enorm(ic+1:end) = NaN;
            erel(ic+1:end) = NaN;
            rnorm(ic+1:end) = NaN;
            rrel(ic+1:end) = NaN;
            q(ic+1:end) = NaN;
        end
        results.iter = k;
        results.checks = checks;
        results.dnorm = dnorm;
        results.drel = drel;
        results.enorm = enorm;
        results.erel = erel;
        results.q = q;
        results.rnorm = rnorm;
        results.rrel = rrel;
        return
    end
end
results.iter = k;
results.checks = checks;
results.dnorm = dnorm;
results.drel = drel;
results.enorm = enorm;
results.erel = erel;
results.q = q;
results.rnorm = rnorm;
results.rrel = rrel;
%
end