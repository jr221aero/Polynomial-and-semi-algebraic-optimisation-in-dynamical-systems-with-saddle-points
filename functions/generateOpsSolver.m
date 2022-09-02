function ops = generateOpsSolver(verboseFlag,debugFlag)
% 
%     Author: Jorge Rodriguez Gutierrez
%     Date: 2022-06-16
%     Function that generates the options for the solvesos function, mainly the solver
%     Inputs:
%         -solverFlag: (stringer) contains the solver to be used
%     Outputs:
%         -ops: options for the solvesos function
%
    ops = sdpsettings('solver','mosek','verbose',verboseFlag,'debug',debugFlag);
%
end