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
nnmax = 20; % maximum iteration for Newton Run
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
