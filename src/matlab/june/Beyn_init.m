function [k,M,N,BeynA0,BeynA1,w_Beyn,w_Beyn_err]=Beyn_init(k_in,N,n,funA,fundA)
%% The first run of Beyn cycle, performs size estimation
    g0=0.0; rho=0.5;    %circular contour
    [g, dg] = circcont_nest(g0, rho, N); % construct contour N

    %--- decide size of random matrix 
    funsize = @(z) trace(funA(z)\fundA(z));
    k_calc = cint(funsize,g,dg); 
    l = floor(real(k_calc))+2;
    disp(sprintf('l=%d',l)); 
    M = rand(n,l);      % dimension of initial arbitrary mat. n x l 
    
    rmw=ones(N,1); % initialize rmw 
    %--- compute Beyn matrices BeynA0, BeynA1
    Nlist =1:N;
    BeynA0=zeros(size(M)); 
    BeynA1=zeros(size(M)); 
    for j=Nlist 
        invAj = funA(g(j))\M;
        BeynA0 = BeynA0 + invAj * dg(j); 
        BeynA1 = BeynA1 + invAj * g(j) * dg(j); 
        if(j==N/2)% store the integral from N/2 run 
            BeynA0_h = BeynA0/(N*i/2);  
            BeynA1_h = BeynA1/(N*i/2); 
        end
    end
    BeynA0 = BeynA0 /(N*i); 
    BeynA1 = BeynA1 /(N*i);               

    l_correct=false;
    while(l_correct==false)
        %--- run svd of BeynA0 and check size to see if l>k 
        [V,Sigma,W] = svd(BeynA0); 
        disp(size(Sigma)); 
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