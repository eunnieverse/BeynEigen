%---------------------------------------------------------------------
%- ContourDistance.m 
%-  plot error as a function of d/s 
%- Yoonkyung Eunnie Lee 
%- last modifled 2015.08.27
%---------------------------------------------------------------------
clear all; close all; warning off;

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
g0 = 0.0; rho = 0.5; 
[gmax,dgmax,s]=NestedContour(g0,rho,Nmax);

Mmax = rand(n);             % square random matrix M0 defined
emax = 1e-3;                % Beyn cutoff error tolerance

%- Initial values for sampling 
l  = n-1;                     % initial number of columns for M
N  = 16;                     % quadrature N initialization 

E = fA.S.E; 
BD = BeynData(N,Mmax,l,gmax,dgmax,emax);

w_0 = E(1); 
d0 =  abs(min(w_0-BD.g)); 
nd = 50; 
dList = linspace(d0,0,nd); 

for(ii=1:nd)
    E(1) = w_0 + ii*min((BD.g-w_0))/nd;
    d = abs(min(E(1)-BD.g));     
    fA.S.E=E; 
end
x = d/s; 

y = max(S_bc.err); 

cfig = scatter(); 

plotsave(cfig, savefigbase, fignum, savejpg, saveeps)

% plot 
%--------------------------------------------------------------------------