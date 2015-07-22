% function main()
% Yoonkyung Eunnie Lee 
% last modified on 2015.07.20

% remove eigenvalues outside contour boundary ?
% use rmw 
% w_Beyn, w_Beyn_new
% w_Newt, w_Newt_new

clear all; close all; warning off; 

showplot=1; 
savefig=0;
saveeps=0;
extension='.jpg';

t0 = cputime;

%% problem definition 
Nmax = 512; %%(2^9=512)
p=2; n=30; 
filebase = sprintf('poly%d_%d',p,n); 
%[coeffs, funA, fundA]=polydef(filebase,p,n); 
load(strcat(filebase,'_fun.mat')); 
g0=0.0;
rho=0.5;
[g,dg]=circcont_nest(g0,rho,Nmax); 
if(showplot==1) cfig=polydefplot(filebase,g(1:2^7)); hold on; end; 

%% starting values 
N=2; 
done=-1; %% indicator to stop simulation 
w_Newt = []; 
i_Newt = []; 

%% limits
w_err_cut = 1e-3; %% Beyn error cutoff 

%% First Beyn Step (initialize)
kk=0; % length(w_Beyn)
while(kk==0) %run Beyn_init until one eigenvalue converges
    N=2*N; 
    [k,N,BeynA0,BeynA1,w_Beyn_new,w_Beyn_err,v_Beyn_new]=...
        Beyn_init(N,g(1:N),dg(1:N),n,funA);
    disp(sprintf('Beyn_init, k=%d,N=%d;',k,N)); 
    i_Beyn=zeros(length(w_Beyn_new),1); %make integer array 
    for(jj=1:length(w_Beyn_err))
        if(w_Beyn_err(jj)<w_err_cut) kk=kk+1; i_Beyn(kk)=jj; end; 
    end
end
if(kk>0) %if anything converged in Beyn step
    i_Beyn=i_Beyn(1:kk); % among all calculated by Beyn, list of the converged
    w_Beyn=w_Beyn_new(i_Beyn);
    v_Beyn=v_Beyn_new(:,i_Beyn);
end
disp('i_Beyn='); disp(i_Beyn'); 
if(showplot==1) scatter(real(w_Beyn_conv),imag(w_Beyn_conv),40,'m.'); end

%% First Newton Step 
[w_Newt,i_Newt,v_Newt]=Newton(w_Beyn,i_Beyn,v_Beyn,funA,fundA);
disp(sprintf('converged eigenvalues=%d',length(w_Newt)));
disp('i_Newt='); disp(i_Newt');

if(showplot==1)
    scatter(real(w_Newt),imag(w_Newt),50,'g');
    savefigname=strcat(filebase,'_result'); 
    if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end;
    if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc');end;
end;

%% Loop
%--- define random matrix M; 
M=rand(n,n);

while(done<0)
% you have to run init with M the same size from the previous.
[k,N,BeynA0,BeynA1,w_Beyn_new,w_Beyn_err,v_Beyn_new,M] = ...
      Beyn(k,N,g,dg,n,funA,M,BeynA0,BeynA1,w_Newt,i_Newt);
    disp(sprintf('Beyn, k=%d, N=%d;',k,N));
    i_Beyn_new = 
    disp('i_Beyn_new='); disp(i_Beyn_new'); 
    
    % add w_Beyn_new to the existing list w_Beyn 
    for(jj=1:length(w_Beyn_new))
        for(ll=1:length(w_Beyn))
            if(abs(w_Beyn_new(jj)-w_Beyn(ll))>w_err_cut)
                w_Beyn=[w_Beyn; w_Beyn_new(jj)]; %add to existing list
                v_Beyn=[v_Beyn v_Beyn_new(:,jj)]; 
            end
        end
    end
    i_Beyn = (1:k)';
    w_Beyn_conv = w_Beyn; 
    v_Beyn_conv = v_Beyn; 
        
    %% Newton Step for previously unconverged 
    [w_Newt_temp,i_Newt_temp,v_Newt_temp] = ...
        Newton(w_Beyn_conv,i_Beyn,v_Beyn_conv,funA,fundA);
   i_Newt=[i_Newt;i_Newt_temp];
   w_Newt=[w_Newt;w_Newt_temp];
   v_Newt=[v_Newt v_Newt_temp]; 
disp(sprintf('converged eigenvalues=%d',length(w_Newt)));
disp('i_Newt='); disp(i_Newt');

done=1; 
     if(length(w_Newt)==k)
%         N=2*N; 
%         k_guess = k; % update the guess for problem size. 
%         else
%             done=1; %yay!
     end
end

e = cputime - t0 ; %% store time 
disp(sprintf('elapsed time=%e',e));