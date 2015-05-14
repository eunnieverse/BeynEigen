%%% Test BeynEigen.m 
    clear all; 
    close all; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% run / load polyeigdef 
    matfilebase = 'polyeig20150505_1'; 
    % polyeigdef(matfilebase, 2, 6);  %% matfilebase, pp, mm 
    %%% load the mfile containing A0,A1,A2,...Ap,Alist, pp,mm
    load(strcat(matfilebase,'.mat'));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% define contour g, g', and integrand function
    
    g0 = 0.0; %center
    rho= 1; %radius 
    g = @(theta) g0 + rho* (cos(theta) + 1i*sin(theta));
    gprime = @(theta) rho* (-sin(theta)+ 1i*cos(theta));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% define integrand function and its parametrization using g, g' 
    f_example = @(z) exp(z)./z;    
    f_numeig  = @(z) trace(inv(funA(z))*fundA(z)); 
    f_eigen   = @(z) inv(funA(z)); 

    fp_example = @(t) f_example(g(t)).*gprime(t);
    fp_numeig  = @(t) f_numeig(g(t)).*gprime(t); 
    fp_eigen   = @(t) f_eigen(g(t)).*gprime(t);  
    figure()
     fplot(f_numeig,[0 2*pi])
%     hold on;
figure()
     fplot(fp_numeig,[0 2*pi])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %% Compute GK Quadrature
    [qgk, errgk] = quadgk(fp_numeig,0,2*pi) %% compute using quadgk 
    q1 = integral(@(t) fp_numeig(t),0,2*pi) %% compute using integral
    %%q1 = integral(@(t) fun(g(t)).*gprime(t),0,2*pi) %% compute using integral
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% Plot 
%     thet = linspace(0,2*pi,100); 
%     cfig = figure();
%     scatter(real(E),imag(E),100,'*'); 
%     hold on; 
%     plot(real(g(thet)),imag(g(thet)),'b.-','Linewidth',1.5); 
%     hold off; 
%     xL = [-max(abs(E))*1.2 max(abs(E))*1.2]; 
%     yL = [-max(abs(E))*1.2 max(abs(E))*1.2];
%     xlim(xL); ylim(yL); 
%     line([0 0], xL,'Color','k','Linewidth',1.5);
%     line(yL, [0 0],'Color','k','Linewidth',1.5);
%     axis square;    
%     title(sprintf('Eigenvalues plotted with the Contour'));
%     xlabel('Re(w)');ylabel('Im(w)');
%     savefigname=strcat(matfilebase,'_beyn');
%     saveas(cfig, strcat(savefigname,'.jpg'));
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%function main()
%end %%function main