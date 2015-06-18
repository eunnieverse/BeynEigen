function [v,omega]=Beyn(funA,M,gamma,tol)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Step 2: Compute BeynA0, BeynA1 from circcontour 
    [n,l]=size(M); 
    f_BeynA0 = @(z) funA(z)\M; 
    f_BeynA1 = @(z) z*f_BeynA0(z); 
    
    %%% compute the contour integral using a circle around g0, with r=rho 
    BeynA0 = cint(f_BeynA0,gamma);
    BeynA1 = cint(f_BeynA1,gamma);   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Step 3: next , compute the SVD of BeynA0 . 
    [V,Sigma,W] = svd(BeynA0); 
    Sigma = Sigma(1:l,1:l);
    V = V(:,1:l); 
    %%% disp(Sigma);
    %%% (conj(transpose(V0))*V0); (V0 and W0 are unitary. ) 
    %%% V0 : size m x l , Sigma0: l x l, W0: l x l 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Step 4: perform a rank test for Sigma and cut meaningless zeros
    k=0; %%% tolerance size for meaningful eigenvalues 
    for ii=1:length(Sigma); 
        if (Sigma(ii,ii)>tol) 
            k=k+1;
        end
    end
    disp(k); 
    V0 = V(1:n, 1:k); 
    W0 = W(1:l,1:k); 
    Sigma0=Sigma(1:k,1:k); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Step 5: Linearization matrix B , size k x k 
    B = conj(V0')*BeynA1*W0*inv(Sigma0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Step 6: Solve eigenvalue problem for B 
    [v,omega]= eig(B);
    omega = diag(omega); 
end
function sum = cint(fun, gamma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% function contintegral: 
%%%                 given a matrix-valued function fn
%%%                 and contour gamma, 
%%%                 return discretized contour using trapezoidal rule
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sum = 0; 
    N = length(gamma); %%contour length 
    for ii=1:N
        sum = sum + fun(gamma(ii)); 
    end
    sum = sum/(N*i);
end %end circcontour
