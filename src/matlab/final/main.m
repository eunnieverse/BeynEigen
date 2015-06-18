%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% NewtonEigen
%%%%% Yoonkyung Eunnie Lee 
%%%%% matlab program to find the eigenvalue of a nonlinear eigenproblem
%%%%% using Newton Method and Beyn's contour integral method together
%%%%% Last Updated 2015.05.24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% main 
%%%%% allow cputime measurement while calling beyn and newton routines. 
clear all; 
close all; 

%% create or load funA and newtA
filebase = 'poly2_100'; 
[coeffs, funA, fundA]=polydef(filebase,3,60); 
%load(strcat(filebase,'_fun')); 

%---------------------------------------------------------------------
%% Beyn Step 
%---------------------------------------------------------------------
%%% necessary improvements : elliptical contour, 
%%% 
%%% parameters
g0 = 0;     %%% center of contour , complex 
rho =0.51;  %%% radius of contour , real 
tolr=1e-4;  %%% rank tolerance to drop near-zero eigenvalues
l = 25;         %%% sampling matrix dimension  
Nmax=150; 
%l = n;
%M = eye(n); 
%---------------------------------------------------------------------
%% load answers E, X for plotting 
%---------------------------------------------------------------------
m = matfile(strcat(filebase,'_E'));
E = m.E;
X = m.X; 
Esamp=E(find(rho>abs(E))); %exact answer 
Xsamp=X(:,find(rho>abs(E))); 
nE = length(Esamp); 
%---------------------------------------------------------------------
%---------------------------------------------------------------------
%% Error Testing 
%---------------------------------------------------------------------
x_N = 1:Nmax; 
err1=zeros(length(x_N),1); %%initialize error vectors. 
err2=zeros(length(x_N),1); 
 for N=x_N
     M = rand(n,l);  %%% dimension of an arbitrary scaling matrix = n x l    
     [vlist,wlist,k]=Beyn1(funA,M,g0,rho,N,tolr);
     err1(N)=min(abs(Esamp(1)-wlist)); 
     err2(N)=min(abs(Esamp(3)-wlist)); 
 end
 
[gamma,gammap] = circcont(g0,rho,Nmax); %%% contour for plotting
 
cfig = figure()%% plot error 
    semilogy(x_N,err1,'.',x_N,err2,'ro'); 
    legend('e(\omega_1)','e(\omega_2)');
    xlabel('N');ylabel('e(\omega_k)');

%% Plot simple BeynEigen Run 
xLc = [real(g0)-rho*1.5 real(g0)+rho*1.5]; %%xL for contour 
yLc = [imag(g0)-rho*1.5 imag(g0)+rho*1.5]; %%yL for contour
xL = [-max(abs(E))*1.2 max(abs(E))*1.2]; %%xL for whole problem
yL = [-max(abs(E))*1.2 max(abs(E))*1.2]; %%yL for whole problem 
cfig = figure();
    scatter(real(gamma),imag(gamma),30,'.'); %% contour  
    hold on; 
    scatter(real(E),imag(E),100,'b*'); %%answer    
    scatter(real(wlist),imag(wlist),100,'ro');
    xlim(xLc); ylim(yLc); 
    line([0 0], xL,'Color','k','Linewidth',1.5);
    line(yL, [0 0],'Color','k','Linewidth',1.5);
    hold off; 
    axis square; 
    xlabel('Re(w)');ylabel('Im(w)');

%% Store Beyn method result from different N
 for N=10:5:150
     M = rand(n,l);  %%% dimension of an arbitrary scaling matrix = n x l    
     [vlist,wlist,k]=Beyn1(funA,M,g0,rho,N,tolr);
     savefilebase=strcat(filebase,'_Beyn',num2str(N),'.mat');
     save(savefilebase,'g0','rho','N','l','k','vlist','wlist'); 
 end



t0 = cputime;
%%% run task here 
e = cputime - t0 ; 

