function omega=Beyn(funA,M,g0,rho,N,gamma,gammap,tol)
    %---------------------------------------------------------------------
    %%% function Beyn: compute updated guess for omega
    %---------------------------------------------------------------------
    %%% Step 2: Compute BeynA0, BeynA1 from circcontour 
    [n,l]=size(M); 
   f_BeynA0 = @(z) funA(z)\M; 
   f_BeynA1 = @(z) z*f_BeynA0(z); 
    %%% compute the contour integral using a circle around g0, with r=rho 
   BeynA0 = cint(f_BeynA0,gamma,gammap);
   BeynA1 = cint(f_BeynA1,gamma,gammap);     
%    BeynA0 = zeros(n,l);  
%    BeynA1 = zeros(n,l); 
%    N = length(gamma); %%contour length 
%    for ii=1:N
%        invAii = funA(gamma(ii))\M; 
%        Bii = invAii*exp(2*pi*ii/N);
%        BeynA0 = BeynA0 + Bii; 
%        BeynA1 = BeynA1 + Bii * exp(2*pi*ii/N); 
%    end
%    BeynA0 = BeynA0*rho/N;
%    BeynA1 = BeynA1*rho*rho/N + BeynA0*g0; 
    
    %%% Step 3: next , compute the SVD of BeynA0 . 
    [V,Sigma,W] = svd(BeynA0); 
    Sigma = Sigma(1:l,1:l);
    %%V0 = V(:,1:l);    
    %%% (conj(transpose(V0))*V0); (V0 and W0 are unitary. ) 
    %%% V0 : size m x l , Sigma0: l x l, W0: l x l 
    %---------------------------------------------------------------------
    %%% Step 4: perform a rank test for Sigma and cut meaningless zeros
    k=0; %%% tolerance size for meaningful eigenvalues 
    for ii=1:length(Sigma); 
        if (Sigma(ii,ii)>tol) 
            k=k+1;
        end
    end
    disp(k); 
    V0 = V(:,1:k);
    W0 = W(1:l,1:k) ;
    Sigma0=Sigma(1:k,1:k);
    Sigmainv = Sigma0\eye(k);
    %---------------------------------------------------------------------
    %%% Step 5: Linearization matrix B , size k x k 
    B = conj(transpose(V0))*BeynA1*W0*Sigmainv;
    %---------------------------------------------------------------------
    %%% Step 6: Solve linear eigenvalue problem for B 
    omega = eig(B); %% cholsky or qz 
end