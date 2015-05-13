%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% BeynTest
%%%%% Yoonkyung Eunnie Lee 
%%%%% matlab program to find the eigenvalue of a nonlinear eigenproblem
%%%%% using Newton Method and Beyn's contour integral method together
%%%%% 2015.05.12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Tester for BeynEigen
    clear all; 
    close all; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% run / load polyeigdef 
    matfilebase = 'polyeig20150505_1'; 
    % polyeigdef(matfilebase, 2, 6);  %% matfilebase, pp, mm 
    %%% load the mfile containing A0,A1,A2,...Ap,Alist, pp,mm
    load(strcat(matfilebase,'.mat'));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% define contour g, g', and integrand function
    f_numeig  = @(z) trace(inv(funA(z))*fundA(z)); 
    %%% First test with V = eye(n);
    %%% dimension of an arbitrary matrix M = n x l 
    %%% M is an arbitrary matrix to reduce dimension of computation
    l = round(n/2); 
    M = rand(n,l); 
    f_BeynA0 = @(z) inv(funA(z))*M; 
    f_BeynA1 = @(z) z*f_BeynA0(z); 
    rho = 4; 
    g0 = 0;
    N = 200;  %% N = 50 or 100 gave only the first eigenvalue.  
    %% compute the contour integral using a circle around g0, with r=rho 
    BeynA0 = circcontour(f_BeynA0,g0,rho,N);
    BeynA1 = circcontour(f_BeynA1,g0,rho,N);   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% next , compute the SVD of BeynA0 . 
    [V0,Sigma0,W0] = svd(BeynA0); 
    disp(Sigma0);
    (conj(transpose(V0))*V0);
    %%% Matrix operations are left now. The task is to compute 
    %%% S Lambda inv(S) = adj(V_0) A1 W0 inv(Sigma_0)
    
    B = V0*A1*W0*inv(Sigma0)
    
    