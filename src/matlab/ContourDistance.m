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
saveeps=1; savejpg=1;       % choose conditions
savemov=0;                  % save movie 
fignum=1;                   % initialize figure number

%---------------------------------------------------------------------
% Construct an NEP (nonlinear eigenvalue problem) 
% funA, fundA, filename, E, X
%---------------------------------------------------------------------
newdef =1;                 % run polydef if newdef
n = 100;                    % define problem size 
savefigbase = sprintf('contourdistance_n%d_',n); 
fA = NEP(2);                % initialize NEP, 1: polynomial, 2: linear
if(newdef==1);
    fA=NEP_linear(fA,n); save(fA);
else
    fA=NEP_load(fA,sprintf(fA.format2,n)); 
end

%---------------------------------------------------------------------
%- Construct BeynData 
%---------------------------------------------------------------------
N  = 16;                     % quadrature N initialization 
Nmax =N*8;                % maximum size of contour 
ng = log2(Nmax) - log2(N); 
ng0 = log2(N); 

g0 = 0.0; rho =0.5; 
[gmax,dgmax,s,dc,isinside]=NestedContour(g0,rho,Nmax);
%Mmax = eye(n); 
Mmax = rand(n);             % square random matrix M0 defined
emax = 1e-4;                % Beyn cutoff error tolerance

%- Initial values for sampling 
l  = n;                   % initial number of columns for M

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
E = fA.S.E;                     % Eigvals list 
Ein = E( find(dc(E)>0) );       % Eigvals inside contour 

g = BDlist(1).g; 
%-
usedpt=1; % move close to an actual point on the contour at N 
if(usedpt)    
    for ii=1:length(Ein)
        dmin(ii) = min(abs(Ein(ii)-g));
    end    
    [B,I]=sort(dmin,'ascend'); 
    w0 = Ein(I(5)); % Eigval closest to a point in g
    i0 = find(E==w0); 
    
    [B,I]=sort(abs(w0-g)); 
    wmax = BDlist(1).g(I(1));     
else
    [B,I]=sort(dc(Ein),'ascend');
    i0 = find(E==Ein(I(1))); 
    w0 = E(i0);               % Eigval closest to the analytical contour 

    r0 = abs(w0-g0);            % radius for w0
    u = (w0-g0)/r0;             % unit vector in w0 direction from g0
    wmax = u*rho;    
end

wlist = w0 + (wmax - w0) * (1-logspace(-10,0,nd)); 
wlist = flip(wlist); 
dlist = dc(wlist); 

%---------------------------------------------------------------------
%- Run Loop to change fA.S.E(i0), funA, fundA 
%---------------------------------------------------------------------
% each column of x stands for N = 16, 32, .. 
% row stands for moving w points 
for gg=1:ng                             % loop for all points on contour
    for ii=1:nd                         % loop for wlist        
        E(i0) = wlist(ii) ; 
        %E(1)=wlist(ii); 
        fA=NEP_linear_update(fA,E);                 
        [BDlist(ii),Slist(ii)]=Beyn(fA,BDlist(ii),S_f,Slist(ii)); 
        Slist(ii)=sample(Slist(ii),find(isinside(Slist(ii).E))); % discard outside gamma        

        if(Slist(ii).k>0)                                % if Beyn has output         
             x(ii,gg) = dlist(ii)/s(N);             
             [B,I]=sort(Slist(ii).err,'descend');
             y1(ii,gg) = B(1); 
%             y2(ii,gg)= B(3); 
%             y3(ii,gg)= B(5); 
             
        end
    end
end
Eerr = Slist(nd).E(I(1)); % eigval with maximum error 
%---------------------------------------------------------------------
%- Plot error 
%---------------------------------------------------------------------
cfig = figure(); 
for gg=1:ng
    loglog (x(:,gg),y1(:,gg),'Color',clist(gg),'LineWidth',2,'Marker','o');
    hold on; 
    %loglog (x(:,gg),y2(:,gg),'Color',clist(gg),'LineWidth',1.5,'Marker','square');
    %loglog (x(:,gg),y3(:,gg),'Color',clist(gg),'LineWidth',1.5,'Marker','diamond');
    legendlist{gg}=sprintf('N=%4d',2^(gg+ng0-1)); 
    hold on; 
end
legend(legendlist); 
xlabel('d/s'); 
ylabel('maximum absolute error'); 
ylim([min(y1(:))*0.1 max(y1(:))*10]);
xlim([min(x(:))*0.1 max(x(:))*10]);

plotsave(cfig, savefigbase, fignum, savejpg, saveeps);
fignum = fignum +1; 

%---------------------------------------------------------------------
%- Plot solutions 
%---------------------------------------------------------------------
cfig=plot(fA); hold on; plot(Slist(ii));  
scatter(real(wlist),imag(wlist),100,'ro');
scatter(real(Eerr),imag(Eerr),100,'g*'); 
h=legend('solution','computed \omega', 'moved \omega', '\omega w. max error'); 
set(h,'Location','northeast','Fontsize',10); 
plot(BDlist(1));
title(sprintf(['N=%5d'],BDlist(ii).N)); 

plotsave(cfig, savefigbase, fignum, savejpg, saveeps);
fignum = fignum +1; 