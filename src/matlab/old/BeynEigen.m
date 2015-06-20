%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% BeynEigen
%%%%% Yoonkyung Eunnie Lee 
%%%%% matlab program to find the eigenvalue of a nonlinear eigenproblem
%%%%% using Newton Method and Beyn's contour integral method together
%%%%% 2015.05.14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% main for Beyn Method 
function Beynmain()
    clear all; close all; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Step 0: Run or Load data for A(z) 
    matfilebase='polyeig20150515';
    % polyeigdef(matfilebase, 2, 6);  %% matfilebase, pp, mm 
    load(strcat(matfilebase,'.mat'));
    l = 15;  %%% currently this is arbitrary 
    
    rho = 1; 
    g0 = 0; 
    N = 150;   %% N = 50 or 100 gave only the first eigenvalue.  
    tol = 0; 
    
    [v,omega]=Beyn(funA,n,l,g0,rho,N,tol); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% NOT USED YET 
    %%% f_numeig  = @(z) trace(inv(funA(z))*fundA(z)); 
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        
end %% end main
