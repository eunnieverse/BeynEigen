% function main()
% Yoonkyung Eunnie Lee 
% last modified on 2015.06.24
clear all; close all; 

t0 = cputime;
%% gproblem definition 
n=30;
filebase = 'poly2_30'; 
%%[coeffs, funA, fundA]=polydef(filebase,2,30); 
load(strcat(filebase,'_fun.mat')); 
%% starting values 
N=64; 
done=-1; %% indicator to stop simulation 
w_Newt = []; 
i_Newt = []; 
k_guess=60; 
%% limits
w_err_cut = 1e-5; %% Beyn error cutoff 
%% First Beyn Step (initialize)
[k,M,N,BeynA0,BegynA1,w_Beyn,w_Beyn_err]=Beyn_init(k_guess,N,n,funA,fundA);
disp(sprintf('Beyn_init, k=%d,N=%d;',k,N)); 
kk=0; %number of w_Beyn_conv 
i_Beyn_conv=zeros(length(w_Beyn)); %make integer array 
for(jj=1:length(w_Beyn_err))
    if(w_Beyn_err(jj)<w_err_cut)
        kk=kk+1; i_Beyn_conv(kk)=jj; 
    end
end
w_Beyn_conv=w_Beyn(i_Beyn);

% while(kk==0) % repeat Beyn until at least one converge after Beyn 
%     if(kk==0)
%         [k, N, BeynA0, BeynA1, M, w_Beyn, w_Beyn_err]=...
%             Beyn(N,k,BeynA0,BeynA1,w_Beyn,M,n,funA,fundA,w_Newt, ...
%                  i_Newt);
%         disp(sprintf('Beyn, k=%d,N=%d;',k,N)); 
%     end
% end
% %% First Newton Step 
% 
% if(kk>0)
%     w_Beyn_conv=w_Beyn_conv(1:kk); %reduce size ;    
%     i_Beyn_conv=i_Beyn_conv(1:kk); %reduce size ;    
%     
%     [w_Newt,i_Newt]=Newton(w_Beyn_conv,i_Beyn_conv,funA,fundA); 
% 
%     [k, N, BeynA0, BeynA1, M, w_Beyn, w_Beyn_err]=...
%   Beyn(N,k,BeynA0,BeynA1,w_Beyn,M,n,funA,fundA,w_Newt,i_Newt);
% 
% end
% 
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
disp(sprintf('e=%e',e));