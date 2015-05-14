function sum = circcontour(fun, g0, rho, N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% function circcontour: a circular contour integral of a given function  
%%%                 center g0, radius rho, sampling number N decides
%%%                 the equally spaced circular contour gamma. 
%%%                 matrix V is chosen as a n x l random matrix
%%%                 initial value of l is chosen as 2*p*n.
%%%                 l should be greater than k. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    theta = linspace(0,2*pi,N+1); 
    gamma = g0 + rho*(cos(theta) + 1i*sin(theta));
    sum = 0; 
    for ii=1:N
        sum = sum + fun(gamma(ii)); 
    end
    sum = sum;
end %end circcontour