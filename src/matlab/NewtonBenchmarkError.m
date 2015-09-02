%--------------------------------------------------------------------------
%%%%% Newton Convergence Check  
%%%%% Yoonkyung Eunnie Lee 
%%%%% matlab program to find the eigenvalue of a nonlinear eigenproblem
%%%%% using Newton Method and Beyn's contour integral method together
%%%%% 2015.05.15
%--------------------------------------------------------------------------
%%% function [vj, wj, j]= NewtInv2(funA,fundA,v0,w0,nn)
%%% function [wj, j]= NewtInv(funA,fundA,w0,nn) 
% %% housekeeping
% clear all; 
% close all;
% %% create or load funA and newtA
% filebase = 'poly2_100'; 
% load(strcat(filebase,'_fun')); 
% %---------------------------------------------------------------------
% %% Load Beyn Step output
% %---------------------------------------------------------------------
% load(strcat(filebase,'_Beyn10')); 
% %m=matfile(strcat(filebase,'_Beyn10')); 
% %w0list = m.wlist; 
% %v0list = m.vlist; 
% %---------------------------------------------------------------------
% %% load answers E, X for plotting 
% %---------------------------------------------------------------------
% m = matfile(strcat(filebase,'_E'));
% E = m.E;
% X = m.X; 
% Esamp=E(find(rho>abs(E))); %exact answer 
% Xsamp=X(:,find(rho>abs(E))); 
% nE = length(Esamp); 
% %% (1) Use only the eigenvalues 
% nn = 50; %% number of maximum iterations 
% %%% for each eigenvector, store the converged value and the 
% %%% number of Newton iterations required
% 
% %% store into three tables, 
% numN=30; 
% k=19
% wtble=zeros(numN,k); %% stores j, w_newton, error
% jtble=zeros(numN,k);
% etble=zeros(numN,k);
% for ii=2:numN
%     m=matfile(strcat(filebase,'_Beyn',num2str(5*ii)));
%     wlist=m.wlist;
%     for kk=1:k %for each eigenvalue inside contour, run Newton Iteration 
%         w0 = wlist(kk); 
%         [w,j]=NewtInv(funA,fundA,w0,nn); 
%         wtble(ii,kk) = w; 
%         jtble(ii,kk) = j; 
%         etble(ii,kk) = min(abs(Esamp-w)); 
%     end
% end

% x = 10:5:150; 
% cfig=figure()
% semilogy(x,etble(2:30,1),'r.',x,etble(2:30,2),'ro', ... 
%          x,etble(2:30,3),'b.',x,etble(2:30,4),'bo', ...
%          x,etble(2:30,5),'g.',x,etble(2:30,6),'go', ...
%          x,etble(2:30,7),'k.',x,etble(2:30,8),'ko', ...    
%          x,etble(2:30,9),'m.',x,etble(2:30,10),'mo');
% legend('e(\omega1)','e(\omega2)',...
%        'e(\omega3)','e(\omega4)',...
%        'e(\omega5)','e(\omega6)',...
%        'e(\omega7)','e(\omega8)',...
%        'e(\omega9)','e(\omega10)'  ); 
%     xlabel('N');ylabel('e(\omega_k)');


x = 10:5:150; 
cfig=figure()
plot(x,jtble(2:30,1),'r.',x,jtble(2:30,2),'ro', ... 
         x,jtble(2:30,3),'b.',x,jtble(2:30,4),'bo', ...
         x,jtble(2:30,6),'go', ...
         x,jtble(2:30,7),'k.',x,jtble(2:30,8),'ko', ...    
         x,jtble(2:30,9),'m.');
% legend('e(\omega1)','e(\omega2)',...
%        'e(\omega3)','e(\omega4)',...
%        'e(\omega5)','e(\omega6)',...
%        'e(\omega7)','e(\omega8)',...
%        'e(\omega9)','e(\omega10)'  ); 
    xlabel('N');ylabel('number of Newton steps');

%%% (2) to include eigenvector benchmark 
