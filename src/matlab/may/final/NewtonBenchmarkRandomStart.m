%--------------------------------------------------------------------------
%%%%% Newton Convergence Check With Random Initial Guess
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
% 
% %% create or load funA and newtA
%  filebase = 'poly2_100'; 
%  load(strcat(filebase,'_fun')); 
% 
%  rho=0.5; 
% %---------------------------------------------------------------------
% %% load answers E, X for plotting 
% %---------------------------------------------------------------------
% m = matfile(strcat(filebase,'_E'));
% E = m.E;
% X = m.X; 
% Esamp=E(find(rho>abs(E))); %exact answer 
% Xsamp=X(:,find(rho>abs(E))); 
% nE = length(Esamp); 
% 
% nn = 50; %% number of maximum iterations 
% %% store into three tables, 
% k=19; %we know k already. 
% wtble=zeros(numN,k); %% stores j, w_newton, error
% jtble=zeros(numN,k);
% etble=zeros(numN,k);
%  x = logspace(1,-3,50); 
% for ii=1:length(x)
%     for kk=1:k %for each eigenvalue inside contour, run Newton Iteration 
%         w0 = Esamp(kk)+rand*x(ii); 
%         w0 = rand*2*ii/100 + Esamp(kk); 
%         [w,j]=NewtInv(funA,fundA,w0,nn); 
%         wtble(ii,kk) = w; 
%         jtble(ii,kk) = j; 
%         etble(ii,kk) = abs(Esamp(kk)-w); 
%     end
% end

%%% plot number of converged eigenvalues 
k=zeros(length(x),1); 
for ii=1:length(x)
    k(ii)=sum(etble(ii,:)<10^-10);
end

cfig=figure()
semilogx(x,k)
xlim([x(1) x(end)]);

cfig=figure()
loglog(x,etble(:,1),'r.',x,etble(:,2),'ro', ... 
         x,etble(:,3),'b.',x,etble(:,4),'bo', ...
         x,etble(:,5),'g.',x,etble(:,6),'go', ...
         x,etble(:,7),'k.',x,etble(:,8),'ko', ...    
         x,etble(:,9),'m.',x,etble(:,10),'mo');
legend('e(\omega1)','e(\omega2)',...
       'e(\omega3)','e(\omega4)',...
       'e(\omega5)','e(\omega6)',...
       'e(\omega7)','e(\omega8)',...
       'e(\omega9)','e(\omega10)'  ); 
    xlim([x(1) x(end)]);
    xlabel('deviation of \omega_0 from \omega');ylabel('e(\omega_k)');


cfig=figure()
semilogx(x,jtble(:,1),'r.',x,jtble(:,2),'ro', ... 
         x,jtble(:,3),'b.',x,jtble(:,4),'bo', ...
         x,jtble(:,6),'go', ...
         x,jtble(:,7),'k.',x,jtble(:,8),'ko', ...    
         x,jtble(:,9),'m.');
     xlim([x(1) x(end)]);
% legend('e(\omega1)','e(\omega2)',...
%        'e(\omega3)','e(\omega4)',...
%        'e(\omega5)','e(\omega6)',...
%        'e(\omega7)','e(\omega8)',...
%        'e(\omega9)','e(\omega10)'  ); 
    xlabel('deviation of \omega_0 from \omega');
    ylabel('number of Newton steps');

%%% (2) to include eigenvector benchmark 