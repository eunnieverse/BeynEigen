% function main()
% Yoonkyung Eunnie Lee 
% last modified on 2015.06.22

clear all; 
close all; 

t0 = cputime;

%% problem definition 
n=30;
filebase = 'poly2_30'; 
[coeffs, funA, fundA]=polydef(filebase,2,30); 
%% starting values 
N=4; 
done=-1; %% indicator to stop simulation 
w_Newt = []; 
i_Newt = []; 
k_guess=60; 
%% limits
w_err_cut = 1e-5; %% Beyn error cutoff 

%% values that need to be stored outside: M, k, 
while(done<0)
    %% Beyn Step 
    N_in = N; 
    [k,N,w_Beyn,i_Beyn,w_Beyn_err]= Beyn(k_guess,n,funA,fundA,N_in, ...
                                         w_Newt, i_Newt, w_err_cut); 
    i_Beyn = (1:k)';
                                     
    %% Newton Step 
    [w_Newt,i_Newt]=Newton(w_Beyn,i_Beyn,funA,fundA); 
    
    if(length(w_Newt)==k)
        N=2*N; 
        k_guess = k; % update the guess for problem size. 
    else
        done=1; %yay!
    end
end

e = cputime - t0 ; 