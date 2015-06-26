function [k, N, BeynA0, BeynA1, M, w_Beyn, w_Beyn_err]=...
  Beyn(N,k_h,BeynA0_h,BeynA1_h,w_Beyn_h,M_h,n,funA,fundA,w_Newt,i_Newt)
    % inputs:  (k, N, BeynA0, BeynA1, M, w_Beyn_h)_h: Beyn matrix, previous data
    %          k: #eigenvalues, N: #quadrature points
    %          M: (n x l) random matrix
    %          n, funA, fundA: problem definition funA=A(w)
    %          w_Newt, i_Newt: converged eigenvalues that need removal
    % outputs: k, N, BeynA0, BeynA1: Beyn matrix data for N=N; 
    %          cdouble w_Beyn[k]: list of eigenvalues found 
    %          cdouble w_Beyn_err[k]: list of corresponding error
    N = N*2;
    l_h = size(M,2);              % store previous size of l  
    %--- define nested contour
    g0=0.0; rho=0.5;              % circular contour
    [g, dg] = circcont_nest(g0, rho, N);  
    
    %--- define rmw(N): (z-w0)(z-w1).. at every quadrature pt.
    rmw=ones(N,1); %initialize 
    if(w_Newt)
        for j=1:N
            rmw(j) = ones(1,length(w_Newt))(g(j) - w_Newt); %% (z-w0)(z-w1)..
        end
    end %%end if(w_Newt)

    %--- compute Beyn matrices BeynA0, BeynA1 using previous M 
    [BeynA0,BeynA1]=getBeyn(BeynA0_h,BeynA1_h,N,M,n,funA,rmw,g,dg);

    l_correct=false;
    while(l_correct==false)
        %--- run svd of BeynA0 and check size to see if l>k 
        [V,Sigma,W] = svd(BeynA0); 
        s = diag(Sigma(1:l,1:l)); %--- first l elements, column vector s
        k= sum(s > 1e-15);        %--- rank computation 
        if(n<k) error('n<rank k , use higher order than BeynA1'); end
        if(k+2<=l)                %--- l should be larger than k+1
            l_correct=true;       %--- exit while loop 
        else                      %--- increase size of M
            M_add = rand(n,1);    % additional column 
            [UpdateA0,UpdateA1,UpdateA0_h,UpdateA1_h]=updateBeyn(funA,M_add,N,rmw,g,dg); 
            BeynA0 = BeynA0 + UpdateA0; 
            BeynA1 = BeynA1 + UpdateA1; 
            BeynA0_h = BeynA0_h + UpdateA0_h; 
            BeynA1_h = BeynA1_h + UpdateA1_h; 
            M = [M M_add];        % update M  
            l = l+1;              % update l 
        end %% end if..else
    end %% end while(l_correct==false)  
        
    %--- compute w_Beyn
    V0 = V(1:n,1:k);              % cut size of V0, W0, Sigma0 using k.
    W0 = W(1:l,1:k);
    s0 = s(1:k);                  % sigma is in decending order. 
    Sinv = diag(1./s0);           % inv(Sigma0) 
    B = conj(V0')*BeynA1*W0*Sinv; % linearized matrix
    [v_Beyn, w_diag]=eig(B);    
    w_Beyn = diag(w_diag);        % convert to single column
    clear V Sigma W s V0 W0 s0 Sinv B w_diag 
    %--- compute error 
    if(l_h ~= l) %if M was updated, recalculate w_Beyn_h 
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
    end %% end if(l_h ~= l)
        
    w_Beyn_err=abs(w_Beyn-w_Beyn_h);  % error    

end%%function 


    
function [BeynA0,BeynA1]=getBeyn(BeynA0_h,BeynA1_h,N,M,funA,rmw,g,dg)
%% get Beyn matrices for a given N and data from N/2 
    BeynA0 = BeynA0_h*(N/2)*i; % sum from N/2 run, if they are given 
    BeynA1 = BeynA1_h*(N/2)*i; % sum from N/2 run, if they are given
    Nlist = (N/2+1):N; 
    for j=Nlist
        invAj = funA(g(j))\M;
        BeynA0 = BeynA0 + invAj * rmw(j) * dg(j); 
        BeynA1 = BeynA1 + invAj * rmw(j) * g(j) * dg(j); 
    end
    BeynA0 = BeynA0 /(N*i); 
    BeynA1 = BeynA1 /(N*i);       
end%%function 
    
function [BeynA0, BeynA1, BeynA0_h, BeynA1_h]=getBeyn2(N,M,n,funA,rmw,g,dg)
%% get Beyn matrices for a given N without data from N/2
    BeynA0 = zeros(size(M)); 
    BeynA1 = zeros(size(M)); 
    Nlist = 1:N; 
    for j=Nlist
        invAj = funA(g(j))\M;
        BeynA0 = BeynA0 + invAj * rmw(j) * dg(j); 
        BeynA1 = BeynA1 + invAj * rmw(j) * g(j) * dg(j); 
    end
    BeynA0 = BeynA0 /(N*i); 
    BeynA1 = BeynA1 /(N*i);
end%%function 

