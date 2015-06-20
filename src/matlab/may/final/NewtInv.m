function [wj, j]= NewtInv(funA,fundA,w0,nn) 
%---------------------------------------------------------------------
%%% NewtInv: iterate to update w using inverse iteration
%%%          input funA, fundA: matrix-valued function and its gradient
%%%                w0: initial guess for eigenvalue
%%%                nn: number of iterations
%---------------------------------------------------------------------
    wj = w0;
    for j = 1: nn
        wj1 = wj-1/trace(funA(wj)\fundA(wj)); 
        if(abs(wj-wj1)<1e-15)
            return;
        else
            wj = wj1;
        end
    end  
end