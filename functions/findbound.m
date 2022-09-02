function [D,coe,c,genCoe,bound,res] = findbound(c,fdelV,phi,boundFlag,ops,coe,epsdel2V)
% 
%     Author: Jorge Rodriguez Gutierrez
%     Date: 2022-06-16
%     Function that finds the upper or lower bound of the quantity of interest for a
%     system of ODEs using YALMIP with prescribed solver.
%     Inputs:
%         -c: (scalar) sdpvar-type coresponding to the constant for the bound
%         -fdelV: (scalar) dot product between f and the gradient of V
%         -phi: (scalar) average of which is the magnitude of interest
%         -boundFlag: "u" or "l" for upper or lower bound, respectively
%         -ops: options for the solving the problem, mainly the solver that
%         is to be used
%         -coe: coefficients of the monomials of the auxiliary function
%         -epsdel2V: (scalar) product of noise intensity and laplacian of the auxiliary function
%     Outputs:
%         -D: quantitiy that is constrained to be positive
%         -coe: optimized coefficients of the monomials of the auxiliary
%         function
%         -c: optimized bound
%         -genCoe: = [coe' c]
%         -bound: (scalar) estimated bound (upper or lower)
%         -res: (scalar) residual
%
%     Define the default quantities that are optional (case with no noise)
if ~exist('epsdel2V','var')
     epsdel2V = 0; % default value for the product of epsilon and laplacian of the auxiliary function (no noise)
end
%
%     Define the quantity that must be SoS (F) and the quantity to minimize
%     (h) depending on if it is for the upper or lower bound
    switch boundFlag
        case "u" % upper bound
            D = -(epsdel2V+fdelV+phi-c);
            F = sos(D); 
            h = c;
        otherwise % "l" % lower bound
            D = epsdel2V+fdelV+phi-c;
            F = sos(D);
            h = -c;
    end
%
    % Solve the problem
    genCoe = [coe' c];
    [sol,~,~,res] = solvesos(F,h,ops,genCoe);
    % Obtain the bound (upper or lower)
    bound = double(c);
    if sol.problem ~= 0
        fprintf('bound: %d <-- %s\n',bound,sol.info)
    end
%
end