%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% BeynEigen
%%%%% Yoonkyung Eunnie Lee 
%%%%% matlab program to find the eigenvalue of a nonlinear eigenproblem
%%%%% using Beyn's contour integral method 
%%%%% 2015.05.04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% main
function main()
    clear all; 
    close all; 
    matfilebase = 'polyeig20150504'; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% run polyeigdef 
    %polyeigdef(matfilebase, 2, 6);  %% matfilebase, pp, mm 
    %%% load the mfile containing A0,A1,A2,...Ap,Alist, pp,mm
    load(strcat(matfilebase,'.mat'));
    
    nn = 20; %% number of iteration for Newton's method 
    %w0 = rand()+i*rand(); 
    w0 = -3+i; %% initial guess 
   
    %%% run newtonpoly to check 
    eval(strcat('[w,wnlist] = newtonpoly(w0,nn,',Alist,');')); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plot together 
    cfig = figure();
    scatter(real(E),imag(E),'*'); 
    hold on; 
    scatter(real(w0),imag(w0),'r*'); 
    scatter(real(w),imag(w),'ro'); 
    plot(real(wnlist),imag(wnlist),'r'); 
    hold off; 
    title(sprintf('Eigenvalue estimation using Newton method,\n iterated %d times, for order-%d polynomial A(w)',nn,pp));
    xlabel('Re(w)');ylabel('Im(w)');
    savefigname=strcat(matfilebase,'_newton');
    saveas(cfig, strcat(savefigname,'.jpg'));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end %%function main 