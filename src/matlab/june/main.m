% function main()
% Yoonkyung Eunnie Lee 
% last modified on 2015.08.03

%% TEST ERROR CONVERGENCE BASED ON BEYN CUTOFF> 
%   k: size of all Beyn output
%  kk: Beyn_converged 
%  jj: input to Newton (from kk)
% jjj: output from Newton, added to kkk 
% kkk: final answer, stored 
%       k,kk,jj,jjj: only stored for the last loop then discarded 

%% Housekeeping 
clear all; close all; 
warning off; 
showplot=1; 
savefig=0;
saveeps=0;
extension='.jpg';
t0 = cputime;

%% problem definition 
Nmax = 512; %%(2^9=512)
% p=3; n=50; 
p=2;n=100;
filebase = sprintf('poly%d_%d',p,n); 
%[coeffs, funA, fundA]=polydef(filebase,p,n); 
load(strcat(filebase,'_fun.mat')); 
g0=0.0; rho=0.5;
[g,dg]=circcont_nest(g0,rho,Nmax); 
if(showplot==1) 
    [xLc,yLc] = getplotlimits(g); 
    cfig=polydefplot(filebase,g(1:2^7)); hold on; 
end; 

% starting values 
N=2; 
done=-1; %% indicator to stop simulation 
w_Newt = []; 
i_Newt = []; 

% limits
w_err_cut = 1e-6; %% Beyn error cutoff 

% define random matrix M or set it to eye(n); 
%M=rand(n,n);
M = eye(n); 

%% First Beyn Step (initialize)
kk=0; % length(w_Beyn); % kk for zeroth step 
while(kk==0) %run Beyn_init until one eigenvalue converges
    N=2*N; %start from N = 4; 

    [k,N,BeynA0,BeynA1,w_Beyn_new,w_Beyn_err,v_Beyn_new]=...
        Beyn_init(N,g(1:N),dg(1:N),n,funA,M);
    
    disp(sprintf('Beyn_init, k=%d,N=%d;',k,N)); 
    i_Beyn=zeros(length(w_Beyn_new),1); %make integer array 
    for j=1:length(w_Beyn_err)  % is length(w_Beyn_err)==k? 
        if(w_Beyn_err(j)<w_err_cut) 
            kk=kk+1; 
            i_Beyn(kk)=j; 
        end; 
    end
end
if(kk>0) %if anything converged in Beyn step
    i_Beyn=i_Beyn(1:kk); % among all calculated by Beyn, list of the converged
    w_Beyn=w_Beyn_new(i_Beyn);
    v_Beyn=v_Beyn_new(:,i_Beyn);
end
disp(sprintf('Initial Beyn, N=%d; total k=%d, converged kk=%d',N,k,kk));
disp('i_Beyn='); disp(i_Beyn'); 
if(showplot==1) scatter(real(w_Beyn),imag(w_Beyn),70,'b'); end %% currently gives best result 

%% Loop: starting from the first Newton Step 
% initialize final list of converged eigenvalues w and eigenvectors v
w=[];
v=[]; 
loopnum=0; 
while(done<0) 
    loopnum=loopnum+1;
    % Newton Step 
    [w_Newt,i_Newt,v_Newt]=Newton(w_Beyn,i_Beyn,v_Beyn,funA,fundA);
    % Display Newton Step 
    disp(sprintf('fully converged eigenvalues=%d',length(w_Newt)));
    disp('i_Newt='); disp(i_Newt');
%     if(showplot==1)
%         scatter(real(w_Newt),imag(w_Newt),50,'g');
%         savefigname=strcat(filebase,'_result'); 
%         if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end;
%         if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc');end;
%     end;

    % Check overlap and append new answers to w,v
    for jjj=1:length(w_Newt)
        add=1; 
        
        if(abs(w_Newt(jjj)-g0)>rho)
            add=0; 
        end
        if(add==1)
            for kkk=1:length(w)
                if(w_Newt==w(kkk))
                    add=0; 
                end
            end
        end
        if(add==1)
            w = [w; w_Newt(jjj)]; 
%            v = [v v_Newt(:,jjj)];
        end   
    end
    % Beyn Step
    kk=0; 
    [k,N,BeynA0,BeynA1,w_Beyn_new,w_Beyn_err,v_Beyn_new,M] = ...
        Beyn(k,N,g,dg,n,funA,M,BeynA0,BeynA1,w);
    
    i_Beyn_loop=zeros(length(w_Beyn_new),1); %make integer array 
    for j=1:length(w_Beyn_err)  % is length(w_Beyn_err)==k? 
        if(w_Beyn_err(j)<w_err_cut) 
            kk=kk+1; 
            i_Beyn_loop(kk)=j; 
        end; 
    end
    
    if(showplot==1) 
        cfig=polydefplot(filebase,g(1:2^7)); hold on; 
        scatter(real(w),imag(w),40,'g.'); 
        scatter(real(w_Beyn),imag(w_Beyn),40,'r'); % from previous Beyn Step; 
        title(sprintf('loop=%d',loopnum)); 
        %legend ('Newton','previous Beyn');
    end
   
    if(kk>0) %if anything converged in Beyn step
        % Check overlap with previous w_Beyn ??
        i_Beyn=i_Beyn_loop(1:kk); % among all calculated by Beyn, list of the converged
        w_Beyn=w_Beyn_new(i_Beyn);
        v_Beyn=v_Beyn_new(:,i_Beyn);
        if(showplot==1)
            scatter(real(w_Beyn),imag(w_Beyn),70,'b'); %from current Beyn Step;
            %legend ('Newton','previous Beyn','current Beyn');
        end
    end   
    disp(sprintf('Beyn, N=%d; total k=%d, converged kk=%d',N,k,kk));
    disp('i_Beyn='); disp(i_Beyn'); 
    
    if(N==512)
        done=1;
    end
end

e = cputime - t0 ; %% store time 
disp(sprintf('elapsed time=%e',e));

% %% Save error file
% % cputime, absolute error, relative error between iterations
% filename = sprintf('error_%s.dat',filebase); 
% fid=fopen(filename,'w');
% formatSpec= '%8.6e    %8.6e    %8.6e \n';
% for(ii=1:length(data))
%     fprintf(fid,formatSpec,...
%         data(ii,1),data(ii,2),data(ii,3)   );
% end
% fclose(fid);