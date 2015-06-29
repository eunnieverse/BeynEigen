% function main()
% Yoonkyung Eunnie Lee 
% last modified on 2015.06.24
clear all; close all; 

t0 = cputime;
%% problem definition 
Nmax = 512; %%(2^9=512)

p=2; 
n=30;
filebase = sprintf('poly%d_%d',p,n); 
%filebase='poly3_100'; 
%%[coeffs, funA, fundA]=polydef(filebase,2,30); 
load(strcat(filebase,'_fun.mat')); 
g0=0.0;
rho=0.7;
[g,dg]=circcont_nest(g0,rho,Nmax); 
% polydefplot(filebase,g(1:2^7)); 

%% starting values 
N=2^2; 
done=-1; %% indicator to stop simulation 
w_Newt = []; 
i_Newt = []; 
k_guess=60; 
%% limits
w_err_cut = 1e-2; %% Beyn error cutoff 
warning off; 
%% First Beyn Step (initialize)
kk=0; %number of w_Beyn_conv
while(kk==0) %run Beyn_init until one eigenvalue converges 
    [k,N,BeynA0,BegynA1,w_Beyn,w_Beyn_err]=Beyn_init(k_guess,N,g(1:N),dg(1:N),n,funA,fundA);
    disp(sprintf('Beyn_init, k=%d,N=%d;',k,N)); 
    i_Beyn_conv=zeros(length(w_Beyn),1); %make integer array 
    for(jj=1:length(w_Beyn_err))
        if(w_Beyn_err(jj)<w_err_cut) kk=kk+1; i_Beyn_conv(kk)=jj; end; 
    end
    N=2*N; 
end
if(kk>0)
    i_Beyn_conv=i_Beyn_conv(1:kk); 
    w_Beyn_conv=w_Beyn(i_Beyn_conv);
end
disp('i_Beyn_conv='); 
disp(i_Beyn_conv'); 
%% First Newton Step 
[w_Newt,i_Newt]=Newton(w_Beyn_conv,i_Beyn_conv,funA,fundA);
disp(sprintf('converged eigenvalues=%d',length(w_Newt)));
disp('i_Newt=');
disp(i_Newt');

%% Next Beyn Step ;
% while(done<0)
%       
%     
%     %% Beyn Step 
%     N_in = N; 
%     [k,N,w_Beyn,i_Beyn,w_Beyn_err]= Beyn(k_guess,n,funA,fundA,N_in, ...
%                                          w_Newt, i_Newt, w_err_cut); 
%     i_Beyn = (1:k)';
%                                      
%     %% Newton Step 
%     [w_Newt,i_Newt]=Newton(w_Beyn,i_Beyn,funA,fundA); 
%     
%     if(length(w_Newt)==k)
%         N=2*N; 
%         k_guess = k; % update the guess for problem size. 
%     else
%         done=1; %yay!
%     end
% end

e = cputime - t0 ; %% store time 
disp(sprintf('elapsed time=%e',e));