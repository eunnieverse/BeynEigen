function [k, N, w_Beyn, i_Beyn, w_Beyn_err]=...
        Beyn(k_guess,n,funA,fundA,N_in,w_Newt,i_Newt,w_err_cut)
    %  function Beyn
    % inputs:  int k_guess: guess for number of eigenvalues
    %          int N_in: quadrature pts 
    %          double w_err_cutoff: Beyn error cutoff (what to
    %          consider converged) 
    %          funA, fundA
    % outputs: int k, int N:  # of eigenvalues, # of quadrature pt
    %          int n: size of funA 
    %          cdouble w_Beyn[k]: list of eigenvalues found 
    %          cdouble w_Beyn_err[k]: list of corresponding error
    %          int i_Beyn[] list of corresponding eigenvalue index
    % Yoonkyung Eunnie Lee 
    % Last Updated 2015.06.18    

    %--- define contour 
    if(mod(N_in,2)==0) N=N_in; else N=N_in+1; end %N is even. 
    g0=0.0; rho=1.0; %circular contour 
    [g, gp] = circcont(g0, rho, N); % construct contour N
    %--- calculate dimension of the problem 
    k_calc = cint(trace(trace(funA(wj)\fundA(wj)), g,gp);
    l = max(k_guess,k_calc); 
    M = rand(n,l);  % dimension of initial arbitrary mat. n x k 

    %--- compute BeynA0, BeynA1 for original function 
    BeynA0 = getA(0,funA,M,g,gp); 
    BeynA1 = getA(1,funA,g,gp);  
    BeynA0_half = getA(0,funA,M,g(1:2:N),gp(1:2:N)); 
    BeynA1_half = getA(1,funA,M,g(1:2:N),gp(1:2:N)); 
    % %--- compute BeynA0, BeynA1 for reduced w_Newt 
    % BeynA0 = getAmod(0,funA,w_Newt,M,g,gp);
    % BeynA1 = getAmod(1,funA,w_Newt,M,g,gp); 
    % BeynA0_half = getAmod(0,funA,w_Newt,M,g(1:2:N),gp(1:2:N)); 
    % BeynA1_half = getAmod(1,funA,w_Newt,M,g(1:2:N),gp(1:2:N)); 

    %--- Compute the SVD of BeynA0 
    [V,Sigma,W] = svd(BeynA0); 
    Sigma = Sigma(1:l,1:l); % first l elements, diagonal matrix 
    s = diag(Sigma); % column
    k= sum(s > 1e-15);      %% rank k of BeynA0. 
    if(n<k)
        error('n<rank k , use Beyn2 instead of Beyn1');
    elseif(l<k)
        error('l<rank k, redo Beyn1 with increased l');
        return; 
    end   
    %---------------------------------------------------------------------
    V0 = V(1:n,1:k);          %% cut size of V0, W0, Sigma0 using k. 
    W0 = W(1:l,1:k);
    s0=s(1:k); %%sigma is in decending order. 
    Sinv = diag(1./s0); 
    %---------------------------------------------------------------------
    B = conj(V0')*BeynA1*W0*Sinv; %%linearized matrix 
    [vlist, wlist]=eig(B); %% qz algorithm 
    w_Beyn = diag(wlist); %%convert into vector;
    i_Beyn = (1:length(w_Beyn))';
    w_Beyn_err=w_Beyn-w_Beyn_half; 
    %---------------------------------------------------------------------