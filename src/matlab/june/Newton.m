function [w_Newt,i_Newt,v_Newt] = Newton(w_Beyn,i_Beyn,v_Beyn,funA,fundA)
% Function Newton 
% inputs:  cdouble w_Beyn[]: list of eigenvalues input from Beyn 
%          int i_Beyn[]: list of corresponding eigenvalue index
%          funA, fundA
% outputs: cdouble w_Newt[]: list of converged eigenvalues 
%          int i_Newt[]: list of corresponding eigenvalue index 
% Yoonkyung Eunnie Lee 
% Last Updated 2015.06.18
    nnmax = 30; % maximum iteration for Newton Run
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
    