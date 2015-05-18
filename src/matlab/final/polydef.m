function [coeffs, funA, fundA] = polydef(matfilebase, p, n)
    %---------------------------------------------------------------------
    %%% polydef.m: create p-th order polynomial NEP  
    %%%               containing the problem and solution of A(w)v=0. 
    %---------------------------------------------------------------------
    %% generate mfiles to store values
    %%% file containig funA, fundA, n, and newtA 
    m = matfile(sprintf('%s_fun.mat',matfilebase),'Writable',true);
    m.n =n;
    %%% file containing X, E
    mE = matfile(sprintf('%s_E',matfilebase),'Writable',true);
    mE.n=n;          %% matrix size
    mE.p=p;          %% polynomial order 

    %---------------------------------------------------------------------
    %% Generate w, A0, A1, A2
    for jj=0:p
         eval(sprintf('A%d=rand(%d); coeffs{%d+1}=A%d;',jj,n,jj,jj));
         %%% A2 = rand(n), coeffs{2}=A2; 
    end
    %---------------------------------------------------------------------
    %% store A(w) and dA(w)/dw as anonymous functions 
    funAstr = 'funA = @(w) A0'; %% string for funA 
    fundAstr = 'fundA = @(w)';  %% string for fundA = dA/dw
    Alist= 'A0'; 
    for k = 1:p %%polynomial order 
        funAstr = sprintf('%s + A%d*w^(%d)',funAstr,k,k); 
        fundAstr = sprintf('%s + %d*A%d*w^(%d)',fundAstr,k,k,k-1); 
        Alist = sprintf('%s,A%d',Alist,k);
    end
    eval(strcat(funAstr,';')); %% funA = @(w) A0 + A1*w^1 + A2*w^2; 
    eval(strcat(fundAstr,';')); %% fundA = @(w) A0 + 1*A1 + 2*A2*w^1; 
%    newtA = @(wn) wn - 1/trace(funA(wn)\fundA(wn));
%    rcondA = @(wn) rcond(funA(wn)\fundA(wn));
%    condA = @(wn) cond(funA(wn)\fundA(wn));
    
    m.funA =funA; 
    m.fundA=fundA; 
%    m.rcondA = rcondA; 
%    m.newtA = newtA;    
    %---------------------------------------------------------------------
    %% Run polyeig: [X,E,Econd] = polyeig(A0,A1,A2)
    eval(strcat('[X,E,Econd] = polyeig(',Alist,');')); 
    mE.X=X;  mE.E=E;  mE.Econd=Econd; 
end