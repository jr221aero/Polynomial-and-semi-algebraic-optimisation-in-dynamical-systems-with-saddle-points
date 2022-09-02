function ffun = sdp2fun(f,coe)
% 
%     Author: Jorge Rodriguez Gutierrez
%     Date: 2022-08-05
%     Function that convert a sdp variable into a function handle
%     Inputs:
%         -f: function (sdp variable)
%         -coe: coefficients
%     Outputs:
%         -f: function handle
%
    ffun = sdisplay(replace(f,coe,value(coe)));
    ffun = ffun{1};
    ffun = replace(ffun,'*','.*');
    ffun = replace(ffun,'^','.^');
    ffun = eval(['@(internal)' ffun]);
%
end