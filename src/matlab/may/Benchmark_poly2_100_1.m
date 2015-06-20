%--------------------------------------------------------------------------
%%%%% BeynEigen Benchmark Error Check 
%%%%% Yoonkyung Eunnie Lee 
%%%%% matlab program to find the eigenvalue of a nonlinear eigenproblem
%%%%% using Newton Method and Beyn's contour integral method together
%%%%% 2015.05.15
%--------------------------------------------------------------------------
%% housekeeping
clear all; 
close all;
    
%% create or load funA and newtA
filebase = 'poly2_100'; 
load(strcat(filebase,'_fun')); 
%n = 100;  
%p = 2; 
m = matfile(strcat(filebase,'_bench2_r0.3')); %%matfile to save benchmark1
%---------------------------------------------------------------------
%% Beyn Step 
%---------------------------------------------------------------------
%%% list of beyn parameters: g0, rho, N, tol_rank, M
%%% Define contour gamma 
g0 = 0;     %%% center of contour , complex 
rho =0.3;  %%% radius of contour , real 
%N = 150;    %%% N = 50 or 100 gave only the first eigenvalue, integer
%[gamma,gammap] = circcont(g0,rho,N); %%% contour 
    
%%% define input parameters for Beyn
tol_rank=1e-4;  %%% rank tolerance to drop near-zero eigenvalues

%%% choose sampling matrix (simply the identity matrix) 
l = 6;         %%% numcols of arbitrary matrix M  %%%variable
M = rand(n,l);  %%% dimension of an arbitrary scaling matrix = n x l     
%M = eye(n); 
%l = n;

%%% FOR largest 6 eigenvalues, construct a table 
x_N=2:150; 
OmegaVsN = zeros(length(x_N),10); 
ErrorVsN = zeros(length(x_N),10); 
%%% load E for plotting 
m = matfile(strcat(filebase,'_E'));
E = m.E;
Esamp=sort(E(find(rho>abs(E))),'descend'); 

%%% for data containing all Omega obtained from Beyn Method, 
%%% maximum length of Omega 
min_error = zeros(length(x_N),length(Esamp)); 
 for j=x_N
     clear gamma; 
     clear gammap;      
     [gamma,gammap] = circcont(g0,rho,j); %%% contour 
     omegaOut1=Beyn(funA,M,gamma,gammap,tol_rank);
     %% for all items in omegaOut, find the minimum error compared to Esamp(1). 
     for ii=1:length(Esamp)
        min_error(j-1,ii) = min( abs(Esamp(ii) - omegaOut1) ) ; 
     end    
 end
 
 cfig = figure()
    semilogy(x_N,min_error(:,1),'.'); 
    hold on; 
    semilogy(x_N,min_error(:,3),'ro'); 
    semilogy(x_N,min_error(:,5),'b.'); 
    hold off; 
    legend('\omega_1','\omega_2','\omega_3');
    xlabel('N');ylabel('e(\omega_k)');
    
    %---------------------------------------------------------------------
    %% Plot simple BeynEigen Run 
    %---------------------------------------------------------------------
    xLc = [real(g0)-rho*1.5 real(g0)+rho*1.5]; %%xL for contour 
    yLc = [imag(g0)-rho*1.5 imag(g0)+rho*1.5]; %%yL for contour
    xL = [-max(abs(E))*1.2 max(abs(E))*1.2]; %%xL for whole problem
    yL = [-max(abs(E))*1.2 max(abs(E))*1.2]; %%yL for whole problem 

    cfig = figure();
        scatter(real(gamma),imag(gamma),30,'.'); %% contour  
        hold on; 
        scatter(real(E),imag(E),100,'b*'); %%answer
        
        scatter(real(omegaOut1),imag(omegaOut1),100,'ro');
        
%         scatter(omega_re,omega_im,50,'ro','filled');
        %xlim(xL); ylim(yL); 
        xlim(xLc); ylim(yLc); 
        line([0 0], xL,'Color','k','Linewidth',1.5);
        line(yL, [0 0],'Color','k','Linewidth',1.5);
        hold off; 
        axis square; 
        title(sprintf('Eigenvalues for %s',filebase),'Interpreter','none');
        xlabel('Re(w)');ylabel('Im(w)');