function [f] = getGovEqn(x,dimX,mu)

    f = cell(dimX,1);

    f{1}.fun = @(x) -mu*x{1}.val;
    f{2}.fun = @(x) x{2}.val-x{2}.val.^3-mu*x{1}.val.^2;

    for idim = 1:dimX
        f{idim}.val = f{idim}.fun(x);
    end

end