function A = getAfp(x,f,epsilon)
%
%     Author: Jorge Rodriguez Gutierrez
%     Date: 2022-07-17
%     Function that generates the A matrix for the steady-state Fokker-Planck equation,
%     expressed as A*x = b, including the proper boundary conditions. A
%     distinction is made between interior and exterior points. Sparse definition
%     with index assignation. Second order of accuracy
%     Inputs:
%         -x: spatial coordinates
%         -f: function giving the dynamics of the system
%         -epsilon: noise intensity
%     Outputs:
%         -A: matrix for the Fokker-Planck equation
%
    Ntot = x{1}.N*x{2}.N;
    Ntot = Ntot-4; % remove four corners
%    
    numNonSparse = 5*(x{1}.N-2)*(x{2}.N-2)+2*3*(x{1}.N-2)+2*3*(x{2}.N-2);
    I = zeros(1,numNonSparse);
    J = zeros(1,numNonSparse);
    V = zeros(1,numNonSparse);
%    
    index = 1;
%    
    % Interior points
    j=2;
        for i = 2:x{1}.N-1
            k = i + x{1}.N*(j-1)-2; % remove the two corners at the top
            J(index:index+4) = k;
            I(index) = k;
            I(index+1) = k+1;
            I(index+2) = k-1;
            I(index+3) = k+x{1}.N;
            I(index+4) = k-(x{1}.N-1); % remove the top right corner
            V(index) = - 2*epsilon*(1/x{1}.d^2+1/x{2}.d^2);
            V(index+1) = epsilon/x{1}.d^2-0.5*f{1}.val(j,i+1)/x{1}.d;
            V(index+2) = epsilon/x{1}.d^2+0.5*f{1}.val(j,i-1)/x{1}.d;
            V(index+3) = epsilon/x{2}.d^2-0.5*f{2}.val(j+1,i)/x{2}.d;
            V(index+4) = epsilon/x{2}.d^2+0.5*f{2}.val(j-1,i)/x{2}.d;
            index = index + 5;
        end
%
    for j = 3:x{2}.N-2
        for i = 2:x{1}.N-1
            k = i + x{1}.N*(j-1)-2; % remove the two corners at the top
            J(index:index+4) = k;
            I(index) = k;
            I(index+1) = k+1;
            I(index+2) = k-1;
            I(index+3) = k+x{1}.N;
            I(index+4) = k-x{1}.N;
            V(index) = - 2*epsilon*(1/x{1}.d^2+1/x{2}.d^2);
            V(index+1) = epsilon/x{1}.d^2-0.5*f{1}.val(j,i+1)/x{1}.d;
            V(index+2) = epsilon/x{1}.d^2+0.5*f{1}.val(j,i-1)/x{1}.d;
            V(index+3) = epsilon/x{2}.d^2-0.5*f{2}.val(j+1,i)/x{2}.d;
            V(index+4) = epsilon/x{2}.d^2+0.5*f{2}.val(j-1,i)/x{2}.d;
            index = index + 5;
        end
    end
%    
    j=x{2}.N-1;
        for i = 2:x{1}.N-1
            k = i + x{1}.N*(j-1)-2; % remove the two corners at the top
            J(index:index+4) = k;
            I(index) = k;
            I(index+1) = k+1;
            I(index+2) = k-1;
            I(index+3) = k+x{1}.N-1; % remove the bottom left corner
            I(index+4) = k-x{1}.N;
            V(index) = - 2*epsilon*(1/x{1}.d^2+1/x{2}.d^2);
            V(index+1) = epsilon/x{1}.d^2-0.5*f{1}.val(j,i+1)/x{1}.d;
            V(index+2) = epsilon/x{1}.d^2+0.5*f{1}.val(j,i-1)/x{1}.d;
            V(index+3) = epsilon/x{2}.d^2-0.5*f{2}.val(j+1,i)/x{2}.d;
            V(index+4) = epsilon/x{2}.d^2+0.5*f{2}.val(j-1,i)/x{2}.d;
            index = index + 5;
        end
%    
    % Exterior points
    % Right boundary
    i = x{1}.N;
        for j = 2:x{2}.N-1
            k = i + x{1}.N*(j-1)-2; % remove the two corners at the top
            J(index:index+2) = k;
            I(index) = k;
            I(index+1) = k-1;
            I(index+2) = k-2;
            V(index) = f{1}.val(j,i)-1.5*epsilon/x{1}.d;
            V(index+1) = 2*epsilon/x{1}.d;
            V(index+2) = -0.5*epsilon/x{1}.d;
            index = index + 3;
        end
%
    % Left boundary
    i = 1;
        for j = 2:x{2}.N-1
            k = i + x{1}.N*(j-1)-2; % remove the two corners at the top
            J(index:index+2) = k;
            I(index) = k;
            I(index+1) = k+1;
            I(index+2) = k+2;
            V(index) = f{1}.val(j,i)+1.5*epsilon/x{1}.d;
            V(index+1) = -2*epsilon/x{1}.d;
            V(index+2) = 0.5*epsilon/x{1}.d;
            index = index + 3;
        end
%
    % Bottom boundary
    j = x{2}.N;
        for i = 2:x{1}.N-1
            k = i + x{1}.N*(j-1)-3; % remove the two corners at the top and the bottom left corner
            J(index:index+2) = k;
            I(index) = k;
            I(index+1) = k-(x{1}.N-1); % remove the bottom left corner
            I(index+2) = k-(2*x{1}.N-1); % remove the bottom left corner
            V(index) = f{2}.val(j,i)-1.5*epsilon/x{2}.d;
            V(index+1) = 2*epsilon/x{2}.d;
            V(index+2) = -0.5*epsilon/x{2}.d;
            index = index + 3;
        end
%
    % Top boundary
    j = 1;
        for i = 2:x{1}.N-1
            k = i + x{1}.N*(j-1)-1; % remove the upper left corner
            J(index:index+2) = k;
            I(index) = k;
            I(index+1) = k+x{1}.N-1; % remove the upper right corner
            I(index+2) = k+2*x{1}.N-1; % remove the upper right corner
            V(index) = f{2}.val(j,i)+1.5*epsilon/x{2}.d;
            V(index+1) = -2*epsilon/x{2}.d;
            V(index+2) = 0.5*epsilon/x{2}.d;
            index = index + 3;
        end
%            
    A = sparse(J,I,V,Ntot,Ntot);
%
end