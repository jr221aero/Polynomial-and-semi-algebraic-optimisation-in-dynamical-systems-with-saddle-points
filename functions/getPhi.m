function phi = getPhi(x,dimX)
% 
%     Author: Jorge Rodriguez Gutierrez
%     Date: 2022-06-11
%     Function that generates the quantity of interest phi of the system
%     Inputs:
%         -x: independent variables
%         -dimX: number of independent variables
%     Outputs:
%         -phi: quantity of interest
%
    phi = 0;
    for idim = 1:dimX
        phi = phi + x{idim}.val.^2;
    end
%
end