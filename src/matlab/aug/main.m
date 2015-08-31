%---------------------------------------------------------------------
% BeynEigen project main file
% Yoonkyung Eunnie Lee 
% last modified on 2015.08.27
%---------------------------------------------------------------------
% used classes in this program: 
% NEP
% NEPcounter
% BeynData, NestedContour
%---------------------------------------------------------------------
% Housekeeping 
%---------------------------------------------------------------------
clear all; close all; warning off;

showplot=0; 
saveeps=0; savejpg=0;       % choose conditions
savemov=1;                  % save movie 
fignum=1;                   % initialize figure number

%---------------------------------------------------------------------
% Initialize Eigenpairs 
%---------------------------------------------------------------------
S_bc = EigenPairs(2);       % eigenpair list from Beyn, converged
S_nc = EigenPairs(3);       % eigenpair list from Newton at each new step 
S_f  = EigenPairs(4);       % eigenpair list final

%---------------------------------------------------------------------
% Construct an NEP (nonlinear eigenvalue problem) 
% funA, fundA, filename, E, X
%---------------------------------------------------------------------
newdef = 0;                 % run polydef if newdef==1

n = 10;                    % define problem size 
p = 2;                      % polynomial order 
fA = NEP(2);                % initialize NEP, 1: polynomial, 2: linear 

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

%---------------------------------------------------------------------
%- Construct BeynData 
%---------------------------------------------------------------------
Nmax = 512;                % maximum size of contour 
g0 = 0.0; rho = 1; 
[gmax,dgmax,s]=NestedContour(g0,rho,Nmax);

% Mmax = rand(n);             % square random matrix M0 defined
Mmax = eye(n); 
emax = 1e-3;                % Beyn cutoff error tolerance

%- Initial values for sampling 
l  = 5;                     % initial number of columns for M
N  = 16;                     % quadrature N initialization 

BD = BeynData(N,Mmax,l,gmax,dgmax,emax);

c = NEPcounter(1000);    % created counter 

NewtonType = 1; % simple Newton-Raphson       
while BD.N<=Nmax/2
    [BD, S_bc] = Beyn(fA, BD, S_nc, S_bc); 
    if(length(S_bc.err)>1)
        c=add(c,l*1.5,max(S_bc.err)); % maximum error among converged eigenvalues 
    end
    S_nc  = Newton(fA, NewtonType, S_bc, S_nc); 
    S_f = update(S_f, S_f.k + S_nc.k,[S_f.E; S_nc.E], [S_f.V S_nc.V]); 
    
    plot(fA); plot(BD); plot(S_bc); plot(S_nc);  
end
%- check for overlap between S_nc and BD.E_nc 
% add to S_f and return to beyn cycle 


%     % Check overlap and append new answers to w,v
%     for jjj=1:length(w_Newt)
%         add=1; 
%         
%         if(abs(w_Newt(jjj)-g0)>rho)
%             add=0; 
%         end
%         if(add==1)
%             for kkk=1:length(w)
%                 if(w_Newt==w(kkk))
%                     add=0; 
%                 end
%             end
%         end
%         if(add==1)
%             w = [w; w_Newt(jjj)]; 
% %            v = [v v_Newt(:,jjj)];
%         end   
%     end

%         % Check overlap with previous w_Beyn ??
%         i_Beyn=i_Beyn_loop(1:kk); % among all calculated by Beyn, list of the converged
%         w_Beyn=w_Beyn_new(i_Beyn);
%         v_Beyn=v_Beyn_new(:,i_Beyn);
%         if(showplot==1)
%             plotBeyn(w_Beyn); 
%             plotsave(cfig, savefigbase, fignum, savejpg, saveeps);
%             if(savemov);    Mov(fignum)=getframe(cfig);       end; 
%             fignum=fignum+1; 
%         end
%     end   
%     disp(sprintf('Beyn, N=%d; total k=%d, converged kk=%d',N,k,kk));
%     disp('i_Beyn='); disp(i_Beyn'); 
 
log(c); 
if(showplot); plot(count,1); end; % 1 for time, 2 for numSolves, 3 for gflops
if(savemov && fignum>4);
    movie2avi(Mov,strcat(savefigbase,'.avi'),'fps',1); 
end;