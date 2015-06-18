function [vlist, wlist, k]=Beyn1(funA,M,g0,rho,N,tolr)
    %---------------------------------------------------------------------
    %%% function Beyn1: Beyn algorithm 1. for k<n
    %%%                 input funA: NEP A(z)
    %%%                 M: random matrix to reduce problem dimension[n,l]
    %%%                 g0, rho, N: circular contour parametrization
    %%%                 tol: rank tolerance for removing spurious values
    %---------------------------------------------------------------------
    [n,l]=size(M); %% pick up problem dimension 
    [g, gp] = circcont(g0, rho, N); %% construct contour 
    %---------------------------------------------------------------------
    % the following four lines are exactly the same as performing the block
    % below. 
    %    f_BeynA0 = @(z) funA(z)\M; 
    %    f_BeynA1 = @(z) z*funA(z)\M; 
    %    BeynA0 = cint(f_BeynA0,gamma,gammap);
    %    BeynA1 = cint(f_BeynA1,gamma,gammap);     
    %---------------------------------------------------------------------  
    BeynA0 = 0; 
    BeynA1 = 0; 
    for j=1:N
        BeynA0 = BeynA0 + funA(g(j))\M *gp(j); 
        BeynA1 = BeynA1 + funA(g(j))\M *g(j)*gp(j); 
    end
    BeynA0 = BeynA0/(N*i);
    BeynA1 = BeynA1/(N*i);   
    %---------------------------------------------------------------------  
    %%% Compute the SVD of BeynA0 
    [V,Sigma,W] = svd(BeynA0); 
    Sigma = Sigma(1:l,1:l); %% V, W are unitary.
    s = diag(Sigma); 
    %told = max(size(BeynA0))*eps(max(s)); %default tolerance
    %disp(sprintf('default tolerance=%e, given tolr=%e\n',told,tolr));
    k = sum(s > tolr);      %% rank k of BeynA0. 
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
    wlist = diag(wlist); %%convert into vector;
    %---------------------------------------------------------------------   
    %m = matfile('Beyn_result.mat','Writable',true);
    %save('Beyn_result.mat','g0','rho','N','k','l','s','B','vlist','wlist'); 
end