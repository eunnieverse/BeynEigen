function main()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% NewtonEigen
%%%%% Yoonkyung Eunnie Lee 
%%%%% matlab program to find the eigenvalue of a nonlinear eigenproblem
%%%%% using Newton Method and Beyn's contour integral method together
%%%%% 2015.05.04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% main for Newton Method 
    clear all; 
    close all; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% run / load polyeigdef 
    matfilebase = 'polyeig20150512'; 
    polyeigdef(matfilebase, 1, 6);  %% matfilebase, p, n 
    %%% load the mfile containing A0,A1,A2,...Ap,Alist, p,n
    load(strcat(matfilebase,'.mat'));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Newton Step 
    nn = 25; %% number of iteration for Newton's method 
    w0 =2-i; %% initial guess 

    wnlisttmp = w0*ones(nn+1,1);
    wn = w0; 
    ii= 1;     
    while ii <= nn;
        wn = newtA(wn); 
        wnlisttmp(ii+1) = wn; 
        ii=ii+1;
        if(rcond(funA(wn)*fundA(wn))<10^-15)
            break;
        end
    end
    wnlist = zeros(ii,1);
    wnlist = wnlisttmp(1:ii); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plot together 
    cfig = figure();
    scatter(real(E),imag(E),100,'*'); 
    hold on; 
    scatter(real(w0),imag(w0),70,'r*'); 
    scatter(real(wn),imag(wn),50,'ro'); 
    plot(real(wnlist),imag(wnlist),'r.-','MarkerSize',5); 
    hold off;
    xL = [-max(abs(E))*1.2 max(abs(E))*1.2]; 
    yL = [-max(abs(E))*1.2 max(abs(E))*1.2];
    xlim(xL); ylim(yL); 
    line([0 0], xL,'Color','k','Linewidth',1.5);
    line(yL, [0 0],'Color','k','Linewidth',1.5);
    axis square; 
    title(sprintf('Eigenvalue estimation using Newton method,\n iterated %d times, for order-%d polynomial A(w)',ii-1,p));
    xlabel('Re(w)');ylabel('Im(w)');
    savefigname=strcat(matfilebase,'_newton');
    saveas(cfig, strcat(savefigname,'.jpg'));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end %%function main 

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
    scatter(real(E),imag(E),100,'*'); 
    title(sprintf('Eigenvalues for order-%d polynomial A(w)',p));
    xlabel('Re(w)');ylabel('Im(w)');
    xL = [-max(abs(E))*1.2 max(abs(E))*1.2]; 
    yL = [-max(abs(E))*1.2 max(abs(E))*1.2];
    xlim(xL); ylim(yL); 
    line([0 0], xL,'Color','k','Linewidth',1.5);
    line(yL, [0 0],'Color','k','Linewidth',1.5);
    axis square; 
    savefigname=matfilebase;
    saveas(cfig, strcat(savefigname,'.jpg'));
end


function [matA,matdA]=polymat(w0,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% function polymat: computes the value of A(w0) and dA(w0)/dw 
%%%                    at a given w0, where A(w)=sum(Ai*w^i) is the
%%%                    polynomial function. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p = nargin-2; %% max polynomial order 
    n = length(varargin{2});
    matA = varargin{1}; %% A = A0 
    matdA = zeros(n); %% dA/dw = 0 
    for k = 1:p
        matA = matA + varargin{k+1}*w0^(k); %% A = A0 + A1*w0 + .. 
        matdA = matdA + varargin{k+1}*k*w0^(k-1); %%dA/dw=0+1*A1+2*A2*w0
    end
end

function [wn, wnlist] = newtonpoly(w0,nn,varargin) %%(w0,nn,A0,A1,A2) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% newtonpoly.m : iterate v->v-det(A)/det'(A)
%%%            given a size n polynomial eigenfunction A(w),
%%%            initial guess w0, 
%%%            and the iteration number nn, 
%%%            compute the update of w based on newton iteration 
%%%            wn = w0; 
%%%            wn = wn - 1/trace(inv(A)*dAdw): iterate nn times 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p = nargin-3; %% max polynomial order 
    n = length(varargin{1}); 
    wn = w0; %% initial value as given 
    wnlist = zeros(nn+1,1);
    wnlist(1)=w0; 
    for ii = 1:nn %%for each newton iteration
        matA = varargin{1}; %% A = A0 
        matdA = zeros(n); %% dA/dw = 0           
        for k = 1:p %%polynomial order 
            matA = matA + varargin{k+1}*wn^(k); %% A = A0 + A1*w0 + .. 
            matdA = matdA + varargin{k+1}*k*wn^(k-1); %%dA/dw=0+1*A1+2*A2*w0
        end
        invA = inv(matA);
        wn = wn - 1/trace(invA*matdA); %%% update wn
        wnlist(ii+1) = wn; 
    end  
end