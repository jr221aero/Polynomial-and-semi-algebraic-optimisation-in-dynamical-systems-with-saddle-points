function x = getSpatialDomain(x,dimX)
% 
%     Author: Jorge Rodriguez Gutierrez
%     Date: 2022-06-11
%     Function that generates the spatial domain
%     Inputs:
%         -x: maximum and minimum value, and grid spacing of the space variables
%         -dimX: number of independent variables
%     Outputs:
%         -x: spatial domain
%
    x{1}.val = x{1}.min:x{1}.d:x{1}.max;
    x{2}.val = (x{2}.min:x{2}.d:x{2}.max)';
%    
    for idim = 1:dimX
        x{idim}.N = length(x{idim}.val);
    end
%    
    [x{1}.val,x{2}.val] = meshgrid(x{1}.val,x{2}.val);
%
end