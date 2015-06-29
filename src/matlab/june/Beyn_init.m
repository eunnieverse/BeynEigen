function [k,N,BeynA0,BeynA1,w_Beyn,w_Beyn_err,v_Beyn]=Beyn_init(k_in,N,g,dg,n,funA,fundA)
%% The first run of Beyn cycle

%     %% perform size estimation
%     %--- decide size of random matrix 
%     funsize = @(z) trace(funA(z)\fundA(z));
%     k_calc = cint(funsize,g,dg); 
%     disp(sprintf('k_calc=%d',k_calc));
%     l = floor(real(k_calc))+3;
%     disp(sprintf('l=%d',l)); 
%     M = rand(n,l);      % dimension of initial arbitrary mat. n x l 
    

    %--- for initial round, do not use M. use M=eye(n); 
    M = eye(n); 
    l=n; 
    %--- compute Beyn matrices BeynA0, BeynA1
    Nlist=1:N;
    BeynA0=zeros(size(M)); 
    BeynA1=zeros(size(M)); 
    for j=Nlist
        invAj = funA(g(j))\M;
        BeynA0 = BeynA0 + invAj * dg(j); 
        BeynA1 = BeynA1 + invAj * g(j) * dg(j); 
        if(j==N/2)% store the integral from N/2 run 
            BeynA0_h = BeynA0/(N*1i/2);  
            BeynA1_h = BeynA1/(N*1i/2); 
        end
    end
    BeynA0 = BeynA0 /(N*1i); 
    BeynA1 = BeynA1 /(N*1i);               

    %--- compute w_Beyn
    [V,Sigma,W] = svd(BeynA0); 
    s = diag(Sigma(1:l,1:l)); %--- first l elements, column vector s        
    k= sum(s > 1e-15);        %--- rank computation 
    V0 = V(1:n,1:k);              % cut size of V0, W0, Sigma0 using k.
    W0 = W(1:l,1:k);
    s0 = s(1:k);                  % sigma is in decending order. 
    Sinv = diag(1./s0);           % inv(Sigma0) 
    B = conj(V0')*BeynA1*W0*Sinv; % linearized matrix
    [v_Beyn, w_diag]=eig(B);    
    w_Beyn = diag(w_diag);        % convert to single column
    clear V Sigma W s V0 W0 s0 Sinv B w_diag 
    %--- compute w_Beyn_h 
    [V,Sigma,W] = svd(BeynA0_h); 
    s = diag(Sigma(1:l,1:l)); %--- first l elements, column vector s        
    V0 = V(1:n,1:k);          % cut size of V0, W0, Sigma0 using k.
    W0 = W(1:l,1:k);
    s0 = s(1:k);              % sigma is in decending order. 
    Sinv = diag(1./s0);       % inv(Sigma0) 
    B = conj(V0')*BeynA1_h*W0*Sinv; % linearized matrix
    [v_Beyn_h, w_diag]=eig(B);    
    w_Beyn_h = diag(w_diag);  % convert to single column
    clear V Sigma W s V0 W0 s0 Sinv B w_diag                 
    %--- compute error
    w_Beyn_err=zeros(k,1); 
    for ii=1:k
        w_Beyn_err(ii)=min(abs(w_Beyn(ii)-w_Beyn_h));
    end
end%%function 