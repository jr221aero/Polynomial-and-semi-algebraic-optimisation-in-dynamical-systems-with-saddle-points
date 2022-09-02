function [x,c,p,coe,delV,del2V] = generateInputs(dimX,degV)
% 
%     Author: Jorge Rodriguez Gutierrez
%     Date: 2022-06-16
%     Function that generates the inputs necessary for the SoS optimization
%     problem.
%     Inputs:
%         -dimX: (scalar) 2 for 2D problem, 3 for 3D problem or N for an N-dimensional
%         problem 
%         -degV: (scalar) degree of the polynomial auxiliary function
%     Outputs:
%         -x: vector of sdpvar-type that contains the independent variables of
%         the problem of dimension [dimX,1]
%         -c: (scalar) sdpvar-type coresponding to the constant for the bound
%         -p: full polynomial of the auxiliary function
%         -coe: coefficients of the monomials of the auxiliary function
%         -m: monomials of the auxiliary function
%         -delV: gradient of the auxiliary function, dimension [1,dimX]
%         -del2V: (scalar) laplacian of the auxiliary function, dimension
%
    x = sdpvar(dimX,1); % Define the independent variables
    c = sdpvar(1); % Define the constant that will correspond to the bound
    [p,coe,~] = polynomial(x,degV); % Obtain the polynomial for the auxiliary function
    delV= jacobian(p,x); % Compute the gradient of the auxiliary function
    del2V= trace(jacobian(delV',x)); % Compute the laplacian of the auxiliary function
%
end

