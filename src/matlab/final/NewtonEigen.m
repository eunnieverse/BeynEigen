%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% NewtonEigen
%%%%% Yoonkyung Eunnie Lee 
%%%%% matlab program to find the eigenvalue of a nonlinear eigenproblem
%%%%% using Newton Method and Beyn's contour integral method together
%%%%% 2015.05.04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% main for Newton Method 

%%% Strategy 
%%% For a Newton step, you need following inputs:
%%% funA, size of problem (n), numsteps nn, 
%%% a list of initial guesses omega_guess; 
%%% output would be an updated guess list omega. 

%%% Since I'm not gonna switch back and forth between Newton and Beyn, 
%%% What I need is just a good criteria to stop Beyn and move to Newton.
%%% for this I can fix the numerical problem to 60 x 60 or 100 x 100
%%% problem and then not worry about another problem. Let's do that. Save
%%% polyeig_100
%%% polyeig_60 
%%% for quadratic cases. 
%%% Done. 
 

function NewtonBeyn()
    clear all; 
    close all; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% run / load polyeigdef 
    col=['kbrgmckbrgmckbrgmckbrgmckbrgmckbrgmc'];

    matfilebase = 'poly2_100'; 
    %polyeigdef(matfilebase, 2, 60);  %% matfilebase, p, n 
    %%% load the mfile containing A0,A1,A2,...Ap,Alist, p,n
    load(strcat(matfilebase,'_fun.mat'));
    load(strcat(matfilebase,'_E.mat'));
    newtA =@(wn) wn - 1/trace(funA(wn)\fundA(wn));
    load(strcat(matfilebase,'_Beyn15.mat')); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Newton Step 
    nn = 25; %% number of iteration for Newton's method     
    cfig = figure();
    hold on; 
    for jj=1:19 %%test for 10 initial cases     
        w0 =wlist(jj); %% initial guess 
        wnlisttmp = w0*ones(nn+1,1);
        wn = w0; 
        ii= 1;     
        while ii <= nn;
            wn = newtA(wn); 
            wnlisttmp(ii+1) = wn; 
            ii=ii+1;
        end
        wnlist = zeros(ii,1);
        wnlist = wnlisttmp(1:ii); 
        errlist= log(-log(abs(wnlist(1:nn-1)-wnlist(2:nn))./abs(wnlist(1:nn-1)))); 
        plot(2:nn,errlist,'--','Linewidth',1.5,'Color',col(jj)); %%answer
        hold on; 
        
    end
    hold off ;
    xlim([2,nn]);
    ylim([0,4]);
    xlabel('number of iterations');
    ylabel('log(-log(relative error))');
    title('Beyn N=5'); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        semilogy(2:nn,wnlist(1:nn-1)-wnlist(2:nn),'b*'); %%answer    
        %plot(real(wnlist),imag(wnlist),'r.-','MarkerSize',5); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end %%function main 

function NewtonConvergence()
    clear all; 
    close all; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% run / load polyeigdef 
    col=['kbrgmckbrgmckbrgmc'];
    matfilebase = 'poly2_100'; 
    %polyeigdef(matfilebase, 2, 60);  %% matfilebase, p, n 
    %%% load the mfile containing A0,A1,A2,...Ap,Alist, p,n
    load(strcat(matfilebase,'_fun.mat'));
    load(strcat(matfilebase,'_E.mat'));
    newtA =@(wn) wn - 1/trace(funA(wn)\fundA(wn));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Newton Step 
    nn = 50; %% number of iteration for Newton's method     
    cfig = figure();
    hold on; 
    for jj=1:10 %%test for five initial cases     
        w0 =2*(rand-0.5)+2*(rand-0.5)*i; %% initial guess 
        wnlisttmp = w0*ones(nn+1,1);
        wn = w0; 
        ii= 1;     
        while ii <= nn;
            wn = newtA(wn); 
            wnlisttmp(ii+1) = wn; 
            ii=ii+1;
            %if(rcond(funA(wn)*fundA(wn))<10^-18)
            %    break;
            %end
        end
        wnlist = zeros(ii,1);
        wnlist = wnlisttmp(1:ii); 
        errlist= log(-log(abs(wnlist(1:nn-1)-wnlist(2:nn))./abs(wnlist(1:nn-1)))); 
        plot(2:nn,errlist,'--','Linewidth',1.5,'Color',col(jj)); %%answer
        hold on; 
    end
    hold off ;
     xlim([2,nn]);
    ylim([0,4]);
    xlabel('number of iterations');
    ylabel('log(-log(relative error))');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        semilogy(2:nn,wnlist(1:nn-1)-wnlist(2:nn),'b*'); %%answer    
        %plot(real(wnlist),imag(wnlist),'r.-','MarkerSize',5); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end %%function main


function NewtonPlotSingle()
    clear all; 
    close all; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% run / load polyeigdef 
    matfilebase = 'poly2_100'; 
    %polyeigdef(matfilebase, 2, 60);  %% matfilebase, p, n 
    %%% load the mfile containing A0,A1,A2,...Ap,Alist, p,n
    load(strcat(matfilebase,'_fun.mat'));
    load(strcat(matfilebase,'_E.mat'));
    newtA =@(wn) wn - 1/trace(funA(wn)\fundA(wn));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Newton Step 
    nn = 100; %% number of iteration for Newton's method     
    w0 =-0.5+0.6i; %% initial guess 
    g0 = 0; 
    rho = 0.5; 
 
    wnlisttmp = w0*ones(nn+1,1);
    wn = w0; 
    ii= 1;     
    while ii <= nn;
        wn = newtA(wn); 
        wnlisttmp(ii+1) = wn; 
        ii=ii+1;
        if(rcond(funA(wn)*fundA(wn))<10^-18)
            break;
        end
    end
    wnlist = zeros(ii,1);
    wnlist = wnlisttmp(1:ii); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plot together 
    xLc = [real(g0)-rho*1.5 real(g0)+rho*1.5]; %%xL for contour 
    yLc = [imag(g0)-rho*1.5 imag(g0)+rho*1.5]; %%yL for contour
    xL = [-max(abs(E))*1.2 max(abs(E))*1.2]; %%xL for whole problem
    yL = [-max(abs(E))*1.2 max(abs(E))*1.2]; %%yL for whole problem 
    cfig = figure();
        %scatter(real(gamma),imag(gamma),30,'.'); %% contour  
        hold on; 
        scatter(real(E),imag(E),100,'b*'); %%answer    
        scatter(real(w0),imag(w0),70,'r*'); 
        scatter(real(wn),imag(wn),50,'ro'); 
        plot(real(wnlist),imag(wnlist),'r.-','MarkerSize',5); 
        hold off;
        xlim(xLc); ylim(yLc); 
        line([0 0], xL,'Color','k','Linewidth',1.5);
        line(yL, [0 0],'Color','k','Linewidth',1.5);
        hold off; 
        axis square; 
        xlabel('Re(w)');ylabel('Im(w)');
    title(sprintf('Eigenvalue estimation using Newton method,\n iterated %d times, for order-%d polynomial A(w)',ii-1,p));
    savefigname=strcat(matfilebase,'_newton');
    saveas(cfig, strcat(savefigname,'.jpg'));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end %%function main 


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