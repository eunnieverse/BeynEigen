%---------------------------------------------------------------------
%- ContourDistance.m 
%-  plot error as a function of d/s 
%- Yoonkyung Eunnie Lee 
%- last modifled 2015.09.01
%---------------------------------------------------------------------
%- Plot error as a function of d/s for multiple values of s
%- choose s from (N=16, 32, 64, 128, 256) 

% Check for two different distance evaluations: dc, dpt

%- change the first eigenvalue of fA.S.E from w to 
clear all; close all; warning off;
format short; 
showplot=1; 
savefigbase = 'contourdistance'; 
saveeps=0; savejpg=0;       % choose conditions
savemov=0;                  % save movie 
fignum=1;                   % initialize figure number
%---------------------------------------------------------------------
% Construct an NEP (nonlinear eigenvalue problem) 
% funA, fundA, filename, E, X
%---------------------------------------------------------------------
newdef = 0;                 % run polydef if newdef
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
Nmax = 512;                % maximum size of contour 
g0 = 0.0; rho = 1; 
[gmax,dgmax,s,dc,isinside]=NestedContour(g0,rho,Nmax);

Mmax = rand(n);             % square random matrix M0 defined
emax = 1e-3;                % Beyn cutoff error tolerance

%- Initial values for sampling 
l  = n-1;                     % initial number of columns for M
N  = 16;                     % quadrature N initialization 

E = fA.S.E; 
BD = BeynData(N,Mmax,l,gmax,dgmax,emax);


w0 = E(1);
r0 = abs(w0-g0);            % radius for w0
u = (w0-g0)/r0;             % unit vector in w0 direction from g0

nd = 5; 
wlist = linspace(w0,u*rho,nd); 
dlist = dc(wlist); 

cfig=plot(fA); hold on; 
scatter(real(wlist),imag(wlist),40,'r');
plot(BD);


disp(wlist);
disp(dlist); 
%%% distance to analytical contour vs. distance to contour point 
%%% d = dc(w); dpt =min(abs(w0-BD.g));  
% 
% for(ii=1:nd)
%     E(1) = w0 + ii*min((BD.g-w0))/nd;
%     d = abs(min(E(1)-BD.g));     
%     fA.S.E=E; 
% end
% x = ds(w)/s(N);
% y = max(S_bc.err); 
% 
% cfig = scatter(); 
% plotsave(cfig, savefigbase, fignum, savejpg, saveeps)
% 
% % plot 
%--------------------------------------------------------------------------