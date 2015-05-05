%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% BeynEigen
%%%%% Yoonkyung Eunnie Lee 
%%%%% matlab program to find the eigenvalue of a nonlinear eigenproblem
%%%%% using Beyn's contour integral method 
%%%%% 2015.05.04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% main for Newton Method 
function main()
    clear all; 
    close all; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% run / load polyeigdef 
    matfilebase = 'polyeig20150504'; 
    % polyeigdef(matfilebase, 2, 6);  %% matfilebase, pp, mm 
    %%% load the mfile containing A0,A1,A2,...Ap,Alist, pp,mm
    load(strcat(matfilebase,'.mat'));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Newton Step 
    nn = 15; %% number of iteration for Newton's method 
    w0 = -2+i; %% initial guess 

    wnlist = zeros(nn+1,1);
    wnlist(1) = w0 ; 
    wn = w0; 

    ii= 1;     
    while ii <= nn;
        wn = newtA(wn); 
        wnlist(ii+1) = wn; 
        ii=ii+1; 
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Plot together 
    cfig = figure();
    scatter(real(E),imag(E),'*'); 
    hold on; 
    scatter(real(w0),imag(w0),'r*'); 
    scatter(real(wn),imag(wn),'ro'); 
    plot(real(wnlist),imag(wnlist),'r'); 
    hold off; 
    title(sprintf('Eigenvalue estimation using Newton method,\n iterated %d times, for order-%d polynomial A(w)',nn,pp));
    xlabel('Re(w)');ylabel('Im(w)');
    savefigname=strcat(matfilebase,'_newton');
    saveas(cfig, strcat(savefigname,'.jpg'));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end %%function main 