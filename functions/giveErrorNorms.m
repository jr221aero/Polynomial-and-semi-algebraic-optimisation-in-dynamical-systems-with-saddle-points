function [L1,L2,LInf] = giveErrorNorms(u,uExact,h,refinement)
%
%     Author: Jorge Rodriguez Gutierrez
%     Date: 2022-07-14
%     Function that generates the error norms for a given numerical grid
%     solution u and the exact grid solution uExact with grid spacing h
%     Inputs:
%         -u: independent variables
%         -uExact: number of independent variables
%         -h: grid spacing
%         -refinement: grid refinement
%     Outputs:
%         -L1: 1-error norm
%         -L2: 2-error norm
%         -LInf: Inf-error norm
%
    uExact = uExact(1:refinement:end,1:refinement:end);
%
    [Ny,Nx] = size(u);
    Ntot = Nx*Ny;
    u = reshape(u,[Ntot,1]);
    uExact = reshape(uExact,[Ntot,1]);
    diff = u - uExact;
%    
    L1 = h^2 * norm(diff,1);
    L2 = h * norm(diff,2);
    LInf = norm(diff,Inf);
%
end