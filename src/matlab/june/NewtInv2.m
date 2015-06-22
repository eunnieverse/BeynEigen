function [vj, wj, j]= NewtInv2(funA,fundA,v0,w0,nn)
%---------------------------------------------------------------------
%%% NewtInv: iterate to update w and v together using inverse iteration
%%%          input funA, fundA: matrix-valued function and its gradient
%%%                v0, w0: initial guess for eigenvector and eigenvalue
%%%                nn: number of iterations
%---------------------------------------------------------------------
    eH = 1./v0';%%normalization row vector   
    vj = v0;
    wj = w0;
    for j = 1: nn
        xj1 = (funA(wj)\fundA(wj))*vj; 
        wj1 = wj-eH*vj/(eH*xj1); 
        vj1 = xj1/(eH*xj1); 
        %if(abs(wj-wj1)<1e-15)
        if(wj==wj1); 
            return;
        else
            vj = vj1;
            wj = wj1;
        end
    end  
end
