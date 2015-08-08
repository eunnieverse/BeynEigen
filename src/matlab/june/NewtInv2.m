function [xj, wj, j]= NewtInv2(funA,fundA,x0,w0,nn)
%---------------------------------------------------------------------
%%% NewtInv: iterate to update w and v together using inverse iteration
%%%          input funA, fundA: matrix-valued function and its gradient
%%%                x0, w0: initial guess for eigenvector and eigenvalue
%%%                nn: number of iterations
%%% Last Modified: 2015. 08. 04
%---------------------------------------------------------------------
    %%normalization vector at each step . 
    xj = x0;
    wj = w0;
    
    lenx= length(xj);   
    disp(lenx); 
    
    for j = 1: nn
        eH = conj(xj')/norm(xj,2)^2; 
        
        uj1 = funA(wj)\fundA(wj)*xj; 
        temp = eH*uj1; 
        wj1 = wj - eH*xj/temp; 
        xj1 = uj1/temp;
        
        rj1=norm(funA(wj1)*xj1,2)/norm(xj1,2);  % compute residual 
        %if(abs(wj-wj1)<1e-15)
        xj = xj1;
        wj = wj1;
        if(rj1<1e-13); 
            return;
        end
    end  
end

function [xj, wj, j]= NewtResInv(funA,x0,w0,nn)
%---------------------------------------------------------------------
%%% NewtResInv: Nonlinear Residual Inverse Iteration 
%%%                x0, w0: initial guess for eigenvector and eigenvalue
%%%                nn: number of iterations
%%% Last Modified: 2015. 08. 06
%---------------------------------------------------------------------
    %%normalization vector at each step . 
    xj = x0;
    wj = w0;
    lenx= length(xj);   
   
    for j = 1:nn
        eH = conj(xj')/norm(xj,2)^2; 
        
        uj1 = funA(wj)\fundA(wj)*xj; 
        temp = eH*uj1; 
        wj1 = wj - eH*xj/temp; 
        xj1 = uj1/temp;
        
        rj1=norm(funA(wj1)*xj1,2)/norm(xj1,2);  % compute residual 
        %if(abs(wj-wj1)<1e-15)
        xj = xj1;
        wj = wj1;
        if(rj1<1e-13); 
            return;
        end
    end  
end