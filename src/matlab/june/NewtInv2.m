function [xj, wj, j]= NewtInv2(funA,fundA,x0,w0,nn)
%---------------------------------------------------------------------
%%% NewtInv: iterate to update w and v together using inverse iteration
%%%          input funA, fundA: matrix-valued function and its gradient
%%%                v0, w0: initial guess for eigenvector and eigenvalue
%%%                nn: number of iterations
%%% Last Modified: 2015. 08. 04
%---------------------------------------------------------------------
    %%normalization vector at each step . 
    xj = x0;
    lenx= length(xj); 
    wj = w0;
    for j = 1: nn
        ejH = 1./xj'/lenx;
        uj = funA(wj)\fundA(wj)*xj; 
        temp = 1/(ejH*uj); 
        wj1 = wj - temp; 
        xj1 = temp*uj;
        %if(abs(wj-wj1)<1e-15)
        if(wj==wj1); 
            return;
        else
            xj = xj1;
            wj = wj1;
        end
    end  
end
