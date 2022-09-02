function [f] = getGovEqn(x,dimX,mu)
% 
%     Author: Jorge Rodriguez Gutierrez
%     Date: 2022-06-11
%     Function that generates the function f of a dynamical system xdot = f*x
%     Inputs:
%         -x: independent variables
%         -dimX: number of independent variables
%         -mu: parameter of the system
%     Outputs:
%         -f: dynamics of the system
%
    f = cell(dimX,1);
%
    f{1}.fun = @(x) -mu*x{1}.val;
    f{2}.fun = @(x) x{2}.val-x{2}.val.^3-mu*x{1}.val.^2;
%
    for idim = 1:dimX
        f{idim}.val = f{idim}.fun(x);
    end
%
end