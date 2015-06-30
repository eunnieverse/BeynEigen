% function main()
% Yoonkyung Eunnie Lee 
% last modified on 2015.06.29

clear all; close all;
showplot=0; 
savefig=1;
saveeps=1;
extension='.jpg';
t0 = cputime;
%% problem definition 
Nmax = 512; %%(2^9=512)
p=2; n=100; filebase = sprintf('poly%d_%d',p,n); %filebase='poly3_100'; 
%[coeffs, funA, fundA]=polydef(filebase,2,30); 
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
warning off; 
%% First Beyn Step (initialize)
kk=0; %number of w_Beyn_conv
while(kk==0) %run Beyn_init until one eigenvalue converges 
    N=2*N; 
    [k,N,BeynA0,BeynA1,w_Beyn,w_Beyn_err,v_Beyn]=...
        Beyn_init(N,g(1:N),dg(1:N),n,funA);
    disp(sprintf('Beyn_init, k=%d,N=%d;',k,N)); 
    i_Beyn_conv=zeros(length(w_Beyn),1); %make integer array 
    for(jj=1:length(w_Beyn_err))
        if(w_Beyn_err(jj)<w_err_cut) kk=kk+1; i_Beyn_conv(kk)=jj; end; 
    end
end
if(kk>0)
    i_Beyn_conv=i_Beyn_conv(1:kk); 
    w_Beyn_conv=w_Beyn(i_Beyn_conv);
    v_Beyn_conv=v_Beyn(:,i_Beyn_conv);
end
disp('i_Beyn_conv='); 
disp(i_Beyn_conv); 
if(showplot==1) 
    scatter(real(w_Beyn_conv),imag(w_Beyn_conv),40,'m.'); 
end
    %% First Newton Step 
[w_Newt,i_Newt,v_Newt]=Newton(w_Beyn_conv,i_Beyn_conv,v_Beyn_conv,funA,fundA);
disp(sprintf('converged eigenvalues=%d',length(w_Newt)));
disp('i_Newt=');
disp(i_Newt');
if(showplot==1)
    scatter(real(w_Newt),imag(w_Newt),50,'g');
    savefigname=strcat(filebase,'_result'); 
    if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end;
    if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc');end;
end;
%% Next Beyn Step
%--- define random matrix M; 
M=rand(n,kk);
while(done<0)
%---first not using rmw 
% you have to run init with M the same size from the previous. 

[k,N,BeynA0,BeynA1,w_Beyn,w_Beyn_err,v_Beyn,M]=...
      Beyn(k,N,g,dg,n,funA,M,BeynA0,BeynA1,w_Newt,i_Newt);
    disp(sprintf('Beyn, k=%d, N=%d;',k,N));
    disp('i_Beyn_conv='); 
    disp(i_Beyn_conv); 

    i_Beyn_conv = (1:k)';
    %% Newton Step 
    [w_Newt_temp,i_Newt_temp,v_Newt_temp]=Newton(w_Beyn_conv,i_Beyn_conv,v_Beyn_conv,funA,fundA);
%    i_Newt=[i_Newt;i_Newt_temp];
%    w_Newt=[w_Newt;w_Newt_temp];
%    v_Newt=[v_Newt v_Newt_temp]; 

    done=1; 
%     if(length(w_Newt)==k)
%         N=2*N; 
%         k_guess = k; % update the guess for problem size. 
%         else
%             done=1; %yay!
%     end
end

e = cputime - t0 ; %% store time 
disp(sprintf('elapsed time=%e',e));