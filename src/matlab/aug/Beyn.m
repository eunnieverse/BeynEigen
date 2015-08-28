function [BD, S_bc] = Beyn(fA, BD, S_nc, S_bc)
%---------------------------------------------------------------------
%- Beyn contour integration, each Beyn run doubles N and changes l
%---------------------------------------------------------------------
% Variables:    BD.N            # of contour pts. 
%               BD.l            # of columns in random mat M = M0(1:n,1:l)
%---------------------------------------------------------------------    
    l0  = BD.l;                % record beginning l value      
    Nj1  = BD.N;                % record beginning N value 
    Nj   = 2 * Nj1;             % double N    
    BD.N = Nj;                  % update BD.N 
    
    S_b  = EigenPairs(1);       % initialize previous list of Beyn eigvals
    %---------------------------------------------------------------------
    %- Update BD.E_nc (Newton-converged list) and decide initial lj
    if(S_nc.k==0)
        lj = l0;               % maintain l if rmw did not increase
    else
        disp('update BD.k ');
        BD.k = S_nc.k + BD.k;   % update BD.k, BD.E_nc to include S_nc
        BD.E_nc = [BD.E_nc; S_nc.E]; 
        lj = l0 - BD.k + 3 ;   % decrease lj according to size of BD.k 
    end
    
    %- Remove eigenvalues from previously converged Newton step S_nc
    usermw=0; 
    rmw=ones(Nj,1);             % initialize     
    if(usermw==1 && S_nc.k>0)
        for jj=1:Nj
            rmw(jj) = ones(1,BD.k)*(BD.g(jj)-BD.E_nc); %% (z-w0)(z-w1)..
        end 
    end
    %---------------------------------------------------------------------
    %- Run with lj and figure out if lj is sufficient. 
    fixl = 0 ;                              % initialize fixl 
    lj1 = l0;                               % temporary lj1 inside loop 
    while( fixl == 0 )
        if(lj1 >= BD.n); break;     end; 
        if(lj  > BD.n); lj = BD.n;  end;
        BD.l = lj;                          % update BD.l 
        %- S_b from half run (1:Nj/2) 
        if(lj == lj1 && BD.NA==Nj1)         % if l=fixed, N=doubled
            S_b = copy(S_b, S_bc);          % reuse previous BeynA as half
            disp('re-used BeynA sum'); 
        else                                % if first run or l changed
            BD = halfBeynA(BD, fA.funA, rmw);  % create new half 
            [S_b, fixl] = BeynSVD(BD, S_b); 
        end
        
        %- S_bc from full run (1:Nj) 
        BD = fullBeynA(BD, fA.funA, rmw);      % update BD
        [S_bc, fixl] = BeynSVD(BD, S_bc);   % get S_bc from current run 
            lj1 = lj;                           % update lj1
            lj = 2 * lj1;                       % double lj 
        %-----------------------------------------------------------------
        %- Check error 
        S_bc = error(S_bc,S_b);     % update error
        S_bc = converged(S_bc,find(S_bc.err<BD.emax)); % discard unconverged

        if(S_bc.k > 0 )
            disp (sprintf('S_bc.k = %d', S_bc.k));
            break; 
        end
    end
    
    %---------------------------------------------------------------------
    % discard eigenvalues outside the contour w_Beyn.

    disp(sprintf('Beyn Run, N:%d->%d, l:%d->%d, k:%d',Nj1,BD.N,l0,BD.l,S_bc.k));    
    
end

function [S_b,fixl] = BeynSVD(BD, S_b)
    BeynA0 = BD.BeynA0sum/(BD.NA * 1i); 
    BeynA1 = BD.BeynA1sum/(BD.NA * 1i);    

    [V,Sigma,W]=svd(BeynA0);            % SVD 
    k=sum(diag(Sigma(1:BD.l,1:BD.l))>BD.tolw); % check rank k 
    if(BD.l > k); fixl = 1;             % determine whether to fix l
    else          fixl = 0;   end;      %       for the next round 
    
    V0 = V(1:BD.n, 1:k);                % cut size down to  V0(n,k)
    S0 = Sigma(1:k,1:k);                %                   S0(k,k) 
    W0 = W(1:BD.l, 1:k);                %                   W0(l,k)

    Sinv = diag(1./diag(S0));           % Sinv = inv(S0) 
    B = conj(V0') * BeynA1 * W0 * Sinv; % linearized matrix 
    [vtemp, wtemp]=eig(B);
    E = diag(wtemp);                    % convert to single column    
    V = V0*vtemp;                       % size of V should be (n,k); 
    S_b = update(S_b, k, E, V); 
end
