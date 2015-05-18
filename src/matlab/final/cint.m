function sum = cint(fun, gamma,gammap)
    %---------------------------------------------------------------------
    %%% function cint:  compute contour integral using trapezoidal rule 
    %%%               given a matrix-valued function fun, contour gamma
    %%%               return contour integral sum 
    %---------------------------------------------------------------------
    sum = 0; 
    N = length(gamma); %%contour length 
    for ii=1:N
        sum = sum + fun(gamma(ii))*gammap(ii); 
    end
    sum = sum/(N*i);
end %end circcontour