function [k, N, BeynA0, BeynA1, w_Beyn, i_Beyn, w_Beyn_err]=...
        Beyn(k_in,N_in,BeynA0_in,BeynA1_in,...
             n,funA,fundA,...
             w_Newt,i_Newt,w_err_cut)
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
<<<<<<< HEAD
    % Yoonkyung Eunnie Lee
    % Last Updated 2015.06.18
    
    %--- define contour 
    q = ceil(log2(N_in)); %N should be 2^q 
    N = 2^q; 

    g0=0.0; rho=1.0;    %circular contour
    [g, gp] = circcont_nest(g0, rho, N); % construct contour N

    %--- determine the size guess for this problem 
    if(N<8) 
        funsize = @(z) trace(funA(z)\fundA(z));
        k_calc = cint(funsize,g,gp); 
        l = ceil(max(k_calc,k_in)); 
    else
        l = ceil(k_in+2);
    end
    M = rand(n,l);      % dimension of initial arbitrary mat. n x k 
    
    %--- compute BeynA0, BeynA1 & store BeynA0_in, BeynA1_in 
    if(N>4)                  % if previous run data exist
        BeynA0 = BeynA0_in*(N/2)*i; % sum from N/2 run, if they are given 
        BeynA1 = BeynA1_in*(N/2)*i; % sum from N/2 run, if they are given
        Nlist = N/2+1:N; 
    else Nlist =1:N; end            % if this is the first run 
    for j=Nlist
        invAj = funA(g(j))\M;       
        BeynA0 = BeynA0 + invAj * gp(j); 
        BeynA1 = BeynA1 + invAj * g(j) * gp(j); 
        if(j==N/2)% store the integral from N/2 run 
            BeynA0_in = BeynA0/(N*i/2);  
            BeynA1_in = BeynA1/(N*i/2); 
        end
    end
    BeynA0 = BeynA0 /(N*i); 
    BeynA1 = BeynA1 /(N*i);    
    
    % w_Beyn, i_Beyn, w_Beyn_err
     
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