%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% BeynTest
%%%%% Yoonkyung Eunnie Lee 
%%%%% matlab program to find the eigenvalue of a nonlinear eigenproblem
%%%%% using Newton Method and Beyn's contour integral method together
%%%%% 2015.05.12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Tester for BeynEigen
clear all; 
close all; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Step 0: Run or Load data for A(z) 
    matfilebase = 'polyeig20150505_1'; 
    % polyeigdef(matfilebase, 2, 6);  %% matfilebase, pp, mm 
    load(strcat(matfilebase,'.mat'));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Step 2: Compute BeynA0, BeynA1 from circcontour 
    l = n-1;  %%% currently this is arbitrary 
    M = rand(n,l);  %%% dimension of an arbitrary scaling matrix = n x l 
    f_BeynA0 = @(z) inv(funA(z))*M; 
    f_BeynA1 = @(z) z*f_BeynA0(z); 
    rho = 3; 
    g0 = 0;
    N = 200;  %% N = 50 or 100 gave only the first eigenvalue.  
    %%% compute the contour integral using a circle around g0, with r=rho 
    BeynA0 = circcontour(f_BeynA0,g0,rho,N);
    BeynA1 = circcontour(f_BeynA1,g0,rho,N);   
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
    tol_rank = 0; 
    k=0; %%% tolerance size for meaningful eigenvalues 
    for ii=1:length(Sigma); 
        if (Sigma(ii,ii)>tol_rank) 
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plot 
    cfig = figure();
    scatter(real(E),imag(E),100,'*'); 
    title(sprintf('Eigenvalues for order-%d polynomial A(w)',p));
    xlabel('Re(w)');ylabel('Im(w)');
    xL = [-max(abs(E))*1.2 max(abs(E))*1.2]; 
    yL = [-max(abs(E))*1.2 max(abs(E))*1.2];
    xlim(xL); ylim(yL); 
    line([0 0], xL,'Color','k','Linewidth',1.5);
    line(yL, [0 0],'Color','k','Linewidth',1.5);
    hold on; 
    scatter(real(omega),imag(omega),90,'ro'); 
    axis square; 
    savefigname=matfilebase;
    saveas(cfig, strcat(savefigname,'.jpg'));
    %%% Matrix operations are left now. The task is to compute 
    %%% S Lambda inv(S) = adj(V_0) A1 W0 inv(Sigma_0)
    