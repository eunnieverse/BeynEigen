clear all; close all; warning off;format short; 
%---------------------------------------------------------------------
%- ContourDistance.m 
%-  plot error as a function of d/s 
%- Yoonkyung Eunnie Lee 
%- last modifled 2015.09.01
%---------------------------------------------------------------------
%- Plot error as a function of d/s for multiple values of s
%- choose s from (N=16, 32, 64, 128, 256) 
% Check for two different distance evaluations: dc, dpt
clist='krgbcmkrgbcmkrgbcm';

showplot=1; 
savefigbase = 'contourdistance'; 
saveeps=0; savejpg=0;       % choose conditions
savemov=0;                  % save movie 
fignum=1;                   % initialize figure number

%---------------------------------------------------------------------
% Construct an NEP (nonlinear eigenvalue problem) 
% funA, fundA, filename, E, X
%---------------------------------------------------------------------
newdef =0 ;                 % run polydef if newdef
n = 10;                    % define problem size 
fA = NEP(2);                % initialize NEP, 1: polynomial, 2: linear
if(newdef==1);
    fA=NEP_linear(fA,n); save(fA);
else
    fA=NEP_load(fA,sprintf(fA.format2,n)); 
end

%---------------------------------------------------------------------
%- Construct BeynData 
%---------------------------------------------------------------------
Nmax =64;                % maximum size of contour 

g0 = 0.0; rho = 0.8; 
[gmax,dgmax,s,dc,isinside]=NestedContour(g0,rho,Nmax);
Mmax = eye(n); 
%Mmax = rand(n);             % square random matrix M0 defined
emax = 1;                % Beyn cutoff error tolerance

%- Initial values for sampling 
l  = 9;                   % initial number of columns for M
N  = 4;                     % quadrature N initialization 
ng = log2(Nmax) - log2(N); 
ng0 = log2(N); 
%---------------------------------------------------------------------
%- Initialize Eigenpairs 
% k, E(k,1), err(k,1), V(n,k), nj(k,1)
%---------------------------------------------------------------------
nd = 20;                     % number of data points on x axis 
S_f  = EigenPairs(4);       % eigenpair list final

for ii = 1:nd; 
    Slist(ii,1)= EigenPairs(2); 
    BDlist(ii,1) = BeynData(N,Mmax,l,gmax,dgmax,emax);
end

%---------------------------------------------------------------------
%- Create wlist and dlist to move a single eigenvalue towards contour 
%---------------------------------------------------------------------
E = fA.S.E; 
% [B,I]=sort(abs(E),'descend'); 
% E=E(I); 
fA=NEP_linear_update(fA,E); 

w0 = E(1);
%-
usedpt=0; % move close to an actual point on the contour at N 
if(usedpt)
    [B,I]=sort(abs(w0-BDlist(1).g)); 
    wmax = BDlist(1).g(I(1));     
else
    r0 = abs(w0-g0);            % radius for w0
    u = (w0-g0)/r0;             % unit vector in w0 direction from g0
    wmax = u*rho;    
end

wlist = w0 + (wmax - w0) * (1-logspace(-2,0,nd)); 
wlist = flip(wlist); 
dlist = dc(wlist); 

%---------------------------------------------------------------------
%- Plot experiment scheme
%---------------------------------------------------------------------
cfig=plot(fA); hold on; plot(BDlist(1));hold on; 
scatter(real(wlist),imag(wlist),40,'r');

disp(wlist); 

%---------------------------------------------------------------------
%- Run Loop to change fA.S.E(1), funA, fundA 
%---------------------------------------------------------------------
% each column of x stands for N = 16, 32, .. 
% row stands for moving w points 
for gg=1:ng
    for ii=1:nd    
        E(1) = wlist(ii) ; 
        fA=NEP_linear_update(fA,E);                 
        [BDlist(ii),Slist(ii)]=Beyn(fA,BDlist(ii),S_f,Slist(ii));        
        if(Slist(ii).k>0)                                % if Beyn has output         
             x(ii,gg) = dlist(ii)/s(N);
             y(ii,gg) = max(Slist(ii).err); 
        end
    end
end

%---------------------------------------------------------------------
%- Plot error 
%---------------------------------------------------------------------
cfig = figure(); 
for gg=1:ng
    loglog (x(:,gg),y(:,gg),'Color',clist(gg),'LineWidth',1.5,'Marker','o');
    legendlist{gg}=sprintf('N=%4d',2^(gg+ng0-1)); 
    hold on; 
end
title('error'); 
legend(legendlist); 
xlabel('d/s'); 
ylabel('max(error)'); 
plotsave(cfig, strcat(savefigbase,num2str(fignum)), fignum, savejpg, saveeps);
fignum = fignum +1; 


%---------------------------------------------------------------------
%- Plot solutions 
%---------------------------------------------------------------------
cfig=plot(fA); plot(BDlist(ii)); plot(Slist(ii)); 
title(sprintf(['N=%5d'],BDlist(ii).N)); 
plotsave(cfig, strcat(savefigbase,num2str(fignum)), fignum, savejpg, saveeps);
%Mov(fignum)=getframe(cfig); 
fignum = fignum +1; 

%%% distance to analytical contour vs. distance to contour point 
%%% d = dc(w); dpt =min(abs(w0-BD.g));  