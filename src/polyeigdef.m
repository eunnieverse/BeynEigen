%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% BeynEigen
%%%%% Yoonkyung Eunnie Lee 
%%%%% matlab program to find the eigenvalue of a nonlinear eigenproblem
%%%%% using Beyn's contour integral method 
%%%%% 2015.05.04
function polyeigdef(matfilebase, p, n) 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% polyeigdef.m: create polyeigdef.m 
    %%%               containing the problem and solution of T(e)X=0. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% generate mfile and store size values 
    m = matfile(sprintf('%s.mat',matfilebase),'Writable',true);   
    m.p = p;          %% polynomial order
    m.n = n;          %% matrix size
    Alist= 'A0'; 
    for k=1:p
        Alist = sprintf('%s,A%d',Alist,k);
    end %% A0,A1,A2 
    m.Alist = Alist; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Generate w, A0, A1, A2
    for jj=0:p
         eval(sprintf('A%d=rand(%d);',jj,n));         %% A2 = rand(n)
         eval(sprintf('m.A%d=A%d;',jj,jj));            %% m.A2=A2
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% store A(w) and dA(w)/dw as anonymous functions 
    funAstr = 'funA = @(w) A0'; %% string for funA 
    fundAstr = 'fundA = @(w)';  %% string for fundA = dA/dw
    for k = 1:p %%polynomial order 
        funAstr = sprintf('%s + A%d*w^(%d)',funAstr,k,k); 
        fundAstr = sprintf('%s + %d*A%d*w^(%d)',fundAstr,k,k,k-1); 
    end
    eval(strcat(funAstr,';')); %% funA = @(w) A0 + A1*w^1 + A2*w^2; 
    eval(strcat(fundAstr,';')); %% fundA = @(w) A0 + 1*A1 + 2*A2*w^1; 
    newtA = @(wn) wn - 1/trace(inv(funA(wn))*fundA(wn));
    m.funA =funA; m.fundA=fundA; m.newtA = newtA; 
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Run polyeig: [X,E,Econd] = polyeig(A0,A1,A2)
    eval(strcat('[X,E,Econd] = polyeig(',Alist,');')); 
    m.X=X;  m.E=E;  m.Econd=Econd; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plot 
    cfig = figure();
    scatter(real(E),imag(E),'*'); 
    title(sprintf('Eigenvalues for order-%d polynomial A(w)',p));
    xlabel('Re(w)');ylabel('Im(w)');
    savefigname=matfilebase;
    saveas(cfig, strcat(savefigname,'.jpg'));
end