%---------------------------------------------------------------------
% BeynEigen project main file
% Yoonkyung Eunnie Lee 
% last modified on 2015.09.02
%---------------------------------------------------------------------
% included classes: NEP, NEPcounter, BeynData, EigenPairs
%---------------------------------------------------------------------
clear all; close all;
showplot=1; saveeps=1; savejpg=1;       % choose conditions
savemov=1;  fignum=1;                   % save movie 
%---------------------------------------------------------------------
%- Initialize Eigenpairs 
% k, E(k,1), err(k,1), V(n,k), nj(k,1)
%---------------------------------------------------------------------
S_bc = EigenPairs(2);       % eigenpair list from Beyn, converged
S_nc = EigenPairs(3);       % eigenpair list from Newton at each new step 
S_f  = EigenPairs(4);       % eigenpair list final

%---------------------------------------------------------------------
%- Construct NEP (nonlinear eigenvalue problem)
% funA, fundA, filename, E, X
%---------------------------------------------------------------------
newdef = 0;                 % run polydef if newdef==1

n = 100;                    % define problem size 
p = 2;                      % polynomial order 
fA = NEP(1);                % initialize NEP, 1: polynomial, 2: linear 

if(newdef==1);
    switch fA.type
        case 1; fA=NEP_poly(fA,n,p);
        case 2; fA=NEP_linear(fA,n); 
    end
    save(fA);
else
    switch fA.type
        case 1; fA=NEP_load(fA,sprintf(fA.format1,p,n)); 
        case 2; fA=NEP_load(fA,sprintf(fA.format2,n)); 
    end
end;
savefigbase=sprintf('%s_randM',fA.filebase);
%---------------------------------------------------------------------
%- Construct BeynData 
%---------------------------------------------------------------------
Nmax = 2^7;                  % maximum size of contour 2^7=128. 2^9=512.
g0 = 0.0; rho =0.5; 
[gmax,dgmax,s,dc,isinside]=NestedContour(g0,rho,Nmax);
 Mmax = rand(n);            % square random matrix M0 defined
% Mmax = eye(n); 
emax = 1e-2;                 % Beyn cutoff error tolerance

%- Initial values for sampling 
l  = n/2;                      % initial number of columns for M
N  = 4;                     % quadrature N initialization 
BD = BeynData(N,Mmax,l,gmax,dgmax,emax);
c = NEPcounter();    % created counter 
%---------------------------------------------------------------------
%- Newton parameters 
%---------------------------------------------------------------------
breakN=@(ej,ej1) ej>=ej1 ;     % break if step size increases    
NewtonType =2; % simple Newton-Raphson       

%---------------------------------------------------------------------
%- Run Cycle
%---------------------------------------------------------------------
while(BD.N<=Nmax/2)
    %-----------------------------------------------------------------
    %- Beyn Step 
    %-----------------------------------------------------------------
    if(S_nc.k>2 && S_bc.k>2)
        BD.N = BD.N/2 ;
    end    
    [BD, S_bc] = Beyn(fA, BD, S_nc, S_bc);     
    if(S_bc.k>0)                                % if Beyn has output 
        c = add(c,l*1.5,max(S_bc.err));         % Record Operation Count                 
        S_bc=sample(S_bc,find(isinside(S_bc.E))); % discard outside gamma
        %-----------------------------------------------------------------
        %- Newton Step  
        %-----------------------------------------------------------------
        S_nc = zero(S_nc); 
        S_nc = Newton(fA.funA, fA.fundA, NewtonType, S_bc, S_nc, breakN);         
        %S_nc = sample(S_nc,find(isinside(S_nc.E))); % discard outside gamma        
        
        for kk=1:S_nc.k                          % if Newton has output 
            if( (kk==1&&S_f.k==0) || min(S_nc.E(kk)-S_f.E)~=0)
                %- update final eigenpair list after checking overlap 
                S_f.k  =  S_f.k+1;               % update k = k+1
                S_f.E(S_f.k,1)  = S_nc.E(kk);    % record eigval
                S_f.err(S_f.k,1)= S_nc.err(kk);  % record step size                
                if(NewtonType==2 || NewtonType==3 || NewtonType==4)
                    S_f.V(:,S_f.k) = S_nc.V(:,kk);% record eigvec
                end
            end
         end
    end 
    cfig=plot(fA); 
    set(cfig,'Position',[50+fignum*5 50+fignum*5 560 420]); 
    plot(BD); plot(S_bc); plot(S_nc); plot(S_f); title(sprintf(['N=' ...
    '%5d'],BD.N)); 
    
    plotsave(cfig, savefigbase, fignum, savejpg, saveeps);
    Mov(fignum)=getframe(cfig); 
    fignum = fignum +1; 
end

%---------------------------------------------------------------------
%- Benchmarking, Save movie if asked 
%---------------------------------------------------------------------
log(c); 
% if(showplot); figure(); plot(c,1); end; % 1 for time, 2 for numSolves, 3 for gflops
if(savemov && fignum>3);
    movie2avi(Mov,strcat(savefigbase,'.avi'),'fps',1); 
end