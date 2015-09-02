function [BD, S_bc] = Beyn(fA, BD, S_f, S_bc)
%---------------------------------------------------------------------
%- Beyn contour integration, each Beyn run doubles N and modifies l
%---------------------------------------------------------------------
% Variables:    BD.N            # of contour pts. 
%               BD.l            # of columns in random mat M = M0(1:n,1:l)
%---------------------------------------------------------------------    
    l0   = BD.l;                % record beginning l value      
    Nj1  = BD.N;                % record beginning N value 
    Nj   = 2 * Nj1;             % double N    
    BD.N = Nj;                  % update BD.N 
    
    S_b  = EigenPairs(1);       % initialize previous list of Beyn eigvals
    %---------------------------------------------------------------------
    %- Update BD.E (Newton-converged list) and decide initial lj
    if(BD.k==S_f.k)
        lj = l0;                % maintain l if rmw did not increase
        disp(sprintf('   Beyn,   S_nc.k==0, letting lj = l0 = %d',l0)); 
    else
       disp('   Beyn,   New eigenvalues converged from Newton step added to rmw');
        dk = S_f.k-BD.k;        % change in k from last step, + if increased
        BD.k = S_f.k; 
        BD.E = S_f.E;         
        lj = l0 - dk;           % decrease lj by dk 
    end
    
    %- Remove eigenvalues from previously converged Newton step 
    usermw=0; 
    
    if(usermw==1 && BD.k>0) 
        disp('   Beyn,   usermw=1, remove eigenvalues');
        %- Define function rmw(z) 
        rmw = @subfun; 
        function y = subfun(z,k,E)
            %- return y = (z-E(1))*(z-E(2))*...*(z-E(k)) 
            y=1; 
            for kk=1:k
                y=y*(z-E(k)); 
            end
        end % function definition    
        
    end
        
    rmwjj = rmw(BD.g(jj),BD.k,BD.E); 
    %---------------------------------------------------------------------
    %- Run with lj and figure out if lj is sufficient. 
    fixl = 0 ;                              % initialize fixl 
    lj1 = l0 ;                              % temporary lj1 inside loop 
    
    while( fixl == 0 )
        if(lj1 >= BD.n) 
           disp(sprintf('   Beyn,   break since lj1(%d) >= BD.n(%d)',lj1,BD.n));
            break;     
        end
        if(lj  > BD.n)
           disp(sprintf('   Beyn,   set lj(previously %d) = BD.n(%d)',lj,BD.n));
            lj = BD.n;
        end;
        
        BD.l = lj;                          % update BD.l 
        
        %- Set S_b from half run (1:Nj/2) 
        if(lj == lj1 && BD.NA==Nj1 && S_bc.k>0) %condition to reuse previous run
           disp('   Beyn,   re-used BeynA sum'); 
            S_b = copy(S_b, S_bc);            
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
        if(S_bc.k>0)
          S_bc = error(S_bc,S_b);     % update error
          S_bc = sample(S_bc,find(S_bc.err<BD.emax)); % discard unconverged        
          disp (sprintf('   S_bc.k = %d', S_bc.k));
            %break; 
        end
    end        
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
    B = (V0') * BeynA1 * W0 * Sinv; % linearized matrix 
    [vtemp, wtemp]=eig(B);
    E = diag(wtemp);                    % convert to single column    
    V = V0*vtemp;                       % size of V should be (n,k); 
    S_b = update(S_b, k, E, V); 
end