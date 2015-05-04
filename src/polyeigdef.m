%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% BeynEigen
%%%%% Yoonkyung Eunnie Lee 
%%%%% matlab program to find the eigenvalue of a nonlinear eigenproblem
%%%%% using Beyn's contour integral method 
%%%%% 2015.05.01
function polyeigdef(matfilebase, pp, mm) 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% polyeigdef.m: create polyeigdef.m 
    %%%               containing the problem and solution of T(e)X=0. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% generate mfile and store size values 
    m = matfile(sprintf('%s.mat',matfilebase),'Writable',true);   
    m.pp = pp;          %% polynomial order
    m.mm = mm;          %% matrix size
    Alist= 'A0'; 
    for k=1:pp
        Alist = sprintf('%s,A%d',Alist,k);
    end %% A0,A1,A2 
    m.Alist = Alist; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Generate w, A0, A1, A2
    for jj=0:pp
         eval(sprintf('A%d=rand(%d);',jj,mm));         %% A2 = rand(mm)
         eval(sprintf('m.A%d=A%d;',jj,jj));            %% m.A2=A2
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Run polyeig: [X,E,Econd] = polyeig(A0,A1,A2)
    eval(strcat('[X,E,Econd] = polyeig(',Alist,');')); 
    m.X=X;  m.E=E;  m.Econd=Econd; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plot 
    cfig = figure();
    scatter(real(E),imag(E),'*'); 
    title(sprintf('Eigenvalues for order-%d polynomial A(w)',pp));
    xlabel('Re(w)');ylabel('Im(w)');
    savefigname=matfilebase;
    saveas(cfig, strcat(savefigname,'.jpg'));
end