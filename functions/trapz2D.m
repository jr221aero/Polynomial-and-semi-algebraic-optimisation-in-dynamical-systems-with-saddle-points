function int = trapz2D(x,u)
% 
%     Author: Jorge Rodriguez Gutierrez
%     Date: 2022-06-29
%     Function that calculates the 2D integral of a function u along the 2D
%     domain given by the x-coordinates
%     Inputs:
%         -x: independent variables
%         -u: function to integrate
%     Outputs:
%         -int: 2D integral
%
    int = trapz(x{2}.val(:,1),trapz(x{1}.val(1,:),u,2));
%
end