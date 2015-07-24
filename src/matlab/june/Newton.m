function [w_Newt,i_Newt,v_Newt] = Newton(w_Beyn,i_Beyn,v_Beyn,funA,fundA)
% Function Newton 
% inputs:  cdouble w_Beyn[]: list of eigenvalues input from Beyn 
%          int i_Beyn[]: list of corresponding eigenvalue index
%          funA, fundA
% outputs: cdouble w_Newt[]: list of converged eigenvalues 
%          int i_Newt[]: list of corresponding eigenvalue index 
% Yoonkyung Eunnie Lee 
% Last Updated 2015.06.22
% Choose method, choose nnmax 
nnmax = 50; % maximum iteration for Newton Run

methd = 1; 
%%% 1: simple Newton-Raphson
%%% 2: Inverse iteration using eigenvector estimate
%%% 3: Residual inverse iteration 
%%% 4: Block-Newton (not implemented)
    if(methd==1)
        i_Newt = zeros(length(i_Beyn),1); 
        w_Newt = zeros(length(i_Beyn),1); 
        mm=0; %index for i_Newt 
        for ll=1:length(i_Beyn) %index for i_Beyn
            [wj, jj]=NewtInv(funA, fundA, w_Beyn(ll), nnmax);
            if(jj<nnmax)%%if converged to machine precision
                mm=mm+1; 
                i_Newt(mm) = i_Beyn(ll); 
                w_Newt(mm) = wj; 
            end
        end
        i_Newt=i_Newt(1:mm); %converged indices
        w_Newt=w_Newt(1:mm); %converged eigenvalues  
        v_Newt=[]; 
    end
    if(methd==2)
        i_Newt = zeros(length(i_Beyn),1); 
        w_Newt = zeros(length(i_Beyn),1); 
        v_Newt = zeros(length(v_Beyn),length(i_Beyn)); 
        mm=0; %index for i_Newt 
        for ll=1:length(i_Beyn) %index for i_Beyn
            %[wj, jj]=NewtInv(funA, fundA, w_Beyn(ll), nnmax);
            [vj, wj, jj]= NewtInv2(funA,fundA,v_Beyn(:,ll),w_Beyn(ll),nnmax);
            if(jj<nnmax)%%if converged to machine precision
                mm=mm+1; 
                i_Newt(mm) = i_Beyn(ll); 
                w_Newt(mm) = wj; 
                v_Newt(:,mm) = vj(:);
            end
        end
        i_Newt=i_Newt(1:mm); %converged indices
        w_Newt=w_Newt(1:mm); %converged eigenvalues 
        v_Newt=v_Newt(:,1:mm); 
    end
end   
    
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
        %if(wj==wj1)
            return;
        else
            wj = wj1;
        end
    end  
end

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
    
    size(vj)
    size(funA(wj)\fundA(wj)); 

    for j = 1:nn
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
