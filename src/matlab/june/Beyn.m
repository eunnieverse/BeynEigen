function [k,N,BeynA0,BeynA1,w_Beyn,w_Beyn_err,v_Beyn,M]=...
  Beyn(k,N,g,dg,n,funA,M,BeynA0_h,BeynA1_h,w_Newt)
%% Repeating run of Beyn Cycle 
% inputs:  (k, N, BeynA0, BeynA1, M, w_Beyn_h)_h: Beyn matrix, previous data
%          k: #eigenvalues, N: #quadrature points
%          M: (n x l) random matrix
%          n, funA, fundA: problem definition funA=A(w)
%          w_Newt: converged eigenvalues that need removal
% outputs: k, N, BeynA0, BeynA1: Beyn matrix data for N=N; 
%          cdouble w_Beyn[k]: list of eigenvalues found 
%          cdouble w_Beyn_err[k]: list of corresponding error
    usermw=1; 
    l = size(M,2); 
    N = N*2;
    %%--- define rmw(N): (z-w0)(z-w1).. at every quadrature pt.
    rmw=ones(N,1); %initialize 
    if(usermw==1 & length(w_Newt)>0) 
        for j=1:N
            rmw(j) = ones(1,length(w_Newt))*(g(j) - w_Newt); %% (z-w0)(z-w1)..
        end 
    end;    
    %% --- compute Beyn matrices BeynA0, BeynA1
    [BeynA0,BeynA1]=getBeyn(BeynA0_h,BeynA1_h,N,M,funA,rmw,g,dg);
    
    %% --- compute w_Beyn
    [V,Sigma,W] = svd(BeynA0); 
    s = diag(Sigma(1:l,1:l)); %--- first l elements, column vector s        
    k= sum(s > 1e-15);        %--- rank computation 
    V0 = V(1:n,1:k);              % cut size of V0, W0, Sigma0 using k.
    W0 = W(1:l,1:k);
    s0 = s(1:k);                  % sigma is in decending order.
    Sinv = diag(1./s0);           % inv(Sigma0) 
    B = (V0')*BeynA1*W0*Sinv; % linearized matrix, 
    [s_Beyn, w_diag]=eig(B);    
    w_Beyn = diag(w_diag);    % convert to single column
    
    %% now discard eigenvalues outside the contour w_Beyn by step 6 in Beyn integral algorithm 1. 
    %%% for ii=1:length(w_Beyn)
    %%% end
    
    v_Beyn = V0*s_Beyn; 
    clear V Sigma W s V0 W0 s0 Sinv B w_diag 
    %--- compute w_Beyn_h 
    [V,Sigma,W] = svd(BeynA0_h); 
    s = diag(Sigma(1:l,1:l)); %--- first l elements, column vector s        
    V0 = V(1:n,1:k);          % cut size of V0, W0, Sigma0 using k.
    W0 = W(1:l,1:k);
    s0 = s(1:k);              % sigma is in decending order. 
    Sinv = diag(1./s0);       % inv(Sigma0) 
    B = (V0')*BeynA1_h*W0*Sinv; % linearized matrix
    [s_Beyn_h, w_diag]=eig(B); 
    w_Beyn_h = diag(w_diag);  % convert to single column
    clear V Sigma W s V0 W0 s0 Sinv B w_diag           
    %--- compute error
    w_Beyn_err=zeros(k,1); 
    for ii=1:k
        w_Beyn_err(ii)=min(abs(w_Beyn(ii)-w_Beyn_h));
    end
    %%% output eigenvector should always have column length of A. 
    
 end%%function
    
function [BeynA0,BeynA1]=getBeyn(BeynA0_h,BeynA1_h,N,M,funA,rmw,g,dg)
%% get Beyn matrices for a given N and data from N/2 
    BeynA0 = BeynA0_h*(N/2)*1i; % sum from N/2 run, if they are given 
    BeynA1 = BeynA1_h*(N/2)*1i; % sum from N/2 run, if they are given
    Nlist = (N/2+1):N; 
    for j=Nlist
        invAj = funA(g(j))\M;
        BeynA0 = BeynA0 + invAj * rmw(j) * dg(j); 
        BeynA1 = BeynA1 + invAj * rmw(j) * g(j) * dg(j); 
    end
    BeynA0 = BeynA0 /(N*1i); 
    BeynA1 = BeynA1 /(N*1i);       
end%%function 