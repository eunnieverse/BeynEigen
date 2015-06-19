%--------------------------------------------------------------------------
%%%%% Benchmark Test for BeynEigen
%%%%% Yoonkyung Eunnie Lee 
%%%%% matlab program to find the eigenvalue of a nonlinear eigenproblem
%%%%% using Newton Method and Beyn's contour integral method together
%%%%% 2015.05.16
%--------------------------------------------------------------------------

%%% strategy
%%% 1. load polynomial eigenproblem
%%% 2. run BeynTest for N = 1:150 and store the error for first two
%%% eigenvalues 
%%% 3. plot error(\omega) vs. N

%%% conditions 
%%%   first, use random matrix M to be identity (not use any random pick)
%%%   next,  use M as rand(n,n) 
%%%   next, use M as  rand(n,l) where l is chosen to be k+1. 

    %% housekeeping
    clear all; 
    close all;
    
    %% create or load funA and newtA
    filebase = 'poly2_60'; 
    n = 60; 
    p = 2; 
    load(strcat(filebase,'_fun')); 

    %---------------------------------------------------------------------
    %% Beyn Step 
    %---------------------------------------------------------------------
    %%% list of beyn parameters: g0, rho, N, tol_rank, M
    %%% Define contour gamma 
    g0 = 0;     %%% center of contour , complex 
    rho =0.6;  %%% radius of contour , real 
    N = 150;    %%% N = 50 or 100 gave only the first eigenvalue, integer
    
    [gamma,gammap] = circcont(g0,rho,N); %%% contour 
    
    %%% define input parameters for Beyn
    tol_rank=1e-3;  %%% rank tolerance to drop near-zero eigenvalues
    
    l = 30;         %%% numcols of arbitrary matrix M  %%%variable
    M = rand(n,l);  %%% dimension of an arbitrary scaling matrix = n x l     
    %M = eye(n); 
    %l = n;
    omegaOut1=Beyn(funA,M,g0,rho,N,gamma,gammap,tol_rank)
    %%% size of vOut = n x k
    %%% size of omegaOut = k x 1
    omegaOut1_re = real(omegaOut1);
    omegaOut1_im = imag(omegaOut1);





function err = errortest1()
    %%% function to test the error of the eigenvalue after Beyn's method
    %%% as the contour integral N increases (exponential accuracy) 
    
    %%% plot the accuracy of the two largest eigenvalues inside the
    %%% contour. omega1 and omega2 
    err1 = zeros(150,1); 
    err2 = zeros(150,1); 
    
    for Nii=1:150
       omega = 
       err1(Nii) = E(1) - omega(1); 
       err2(Nii) = E(2) - omega(2); 
    end
    %%% error plot: 
    cfig = figure()
    
    
end