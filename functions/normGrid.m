function unorm = normGrid(u,h)
% 
%     Author: Jorge Rodriguez Gutierrez
%     Date: 2022-07-02
%     Function that generates the 2-norm of 2D grid function u with uniform
%     grid spacing h
%     Inputs:
%         -u: grid function
%         -h: uniform grid spacing
%     Otuputs:
%         -unorm: 2-norm of u
%
    unorm = sqrt(h*sumsqr(u));
%
end