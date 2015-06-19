%--------------------------------------------------------------------------
%%%%% BeynEigen
%%%%% Yoonkyung Eunnie Lee 
%%%%% matlab program to find the eigenvalue of a nonlinear eigenproblem
%%%%% using Newton Method and Beyn's contour integral method together
%%%%% 2015.05.15
%--------------------------------------------------------------------------
%%% list of functions 
%%% function funA = polydef(matfilebase, p, n)
%%% function gamma = circcont(g0,rho,N)
%%% function sum = cint(fun, gamma)
%%% function omega = Beyn(funA,M,gamma,tol)
%%% function omega = Newt(newtA,omegaIn,nn)

%function main()
    %% housekeeping
    clear all; 
    close all;
    
    %% create or load funA and newtA
    filebase = 'poly2_60'; 
    n = 60;  
    p = 2; 
    [coeffs, funA, fundA] = polydef(filebase, p, n); 
    %[coeffs,funA,fundA,newtA,rcondA] = polydef(filebase,p,n); 
    load(strcat(filebase,'_fun')); 
    %---------------------------------------------------------------------
    %% Beyn Step 
    %---------------------------------------------------------------------
    %%% list of beyn parameters: g0, rho, N, tol_rank, M
    %%% Define contour gamma 
    g0 = 0;     %%% center of contour , complex 
    rho =0.6;  %%% radius of contour , real 
    N = 150;    %%% N = 50 or 100 gave only the first eigenvalue, integer
    
    [gamma,gammap] = circcont(g0,rho,N); %%% contour 
    
    %%% define input parameters for Beyn
    tol_rank=1e-3;  %%% rank tolerance to drop near-zero eigenvalues
    
    l = 30;         %%% numcols of arbitrary matrix M  %%%variable
    M = rand(n,l);  %%% dimension of an arbitrary scaling matrix = n x l     
    %M = eye(n); 
    %l = n;
    omegaOut1=Beyn(funA,M,g0,rho,N,gamma,gammap,tol_rank)
    %%% size of vOut = n x k
    %%% size of omegaOut = k x 1
    omegaOut1_re = real(omegaOut1);
    omegaOut1_im = imag(omegaOut1);
%     %---------------------------------------------------------------------
%     %% Newton Step 
%     %---------------------------------------------------------------------
%     nn = 10; 
%     omega = Newt(newtA,rcondA,omegaOut1,nn); 
%     omega_re = real(omega); 
%     omega_im = imag(omega); 
    %---------------------------------------------------------------------
    %% Plot 
    %---------------------------------------------------------------------
    %%% load E for plotting 
    m = matfile(strcat(filebase,'_E'));
    E = m.E;
    xLc = [real(g0)-rho*1.5 real(g0)+rho*1.5]; %%xL for contour 
    yLc = [imag(g0)-rho*1.5 imag(g0)+rho*1.5]; %%yL for contour
    xL = [-max(abs(E))*1.2 max(abs(E))*1.2]; %%xL for whole problem
    yL = [-max(abs(E))*1.2 max(abs(E))*1.2]; %%yL for whole problem 

    cfig = figure();
        scatter(real(gamma),imag(gamma),30,'.'); %% contour  
        hold on; 
        scatter(real(E),imag(E),100,'b*'); %%answer
        
        scatter(omegaOut1_re,omegaOut1_im,100,'ro');
        
%         scatter(omega_re,omega_im,50,'ro','filled');
        %xlim(xL); ylim(yL); 
        xlim(xLc); ylim(yLc); 
        line([0 0], xL,'Color','k','Linewidth',1.5);
        line(yL, [0 0],'Color','k','Linewidth',1.5);
        hold off; 
        axis square; 
        title(sprintf('Eigenvalues for %s',filebase),'Interpreter','none');
        xlabel('Re(w)');ylabel('Im(w)');
        
        savefigname=filebase;
        saveas(cfig, strcat(savefigname,'.jpg'));
%end

