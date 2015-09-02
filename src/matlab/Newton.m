function S_nc = Newton(funA, fundA, type, S_bc, S_nc, breakN)
    %-------------------------------------------------------------
    %- run Newton-based Iteration on every entry of S_bc
    %  type 1: Simple Newton-Raphson
    %       2: Inverse iteration using eigenvector estimate
    %       3: Residual inverse iteration
    %-------------------------------------------------------------    
    jmax = 20;                          % maximum iteration     
    %-------------------------------------------------------------
    %- Start For Loop for every eigenpair in S_bc 
    %-------------------------------------------------------------    
    for kk=1:S_bc.k
        switch type
            case 1; [wj,ej,jj]    = Newton1(funA,fundA,breakN,jmax,S_bc.E(kk));
            case 2; [wj,vj,ej,jj] = Newton2(funA,fundA,breakN,jmax,S_bc.E(kk),S_bc.V(:,kk));
            case 3; error('Nonlinear Residual Inverse Iteration not implemented');
            case 4; error('Block Newton not implemented');
        end % switch        
        if(jj<jmax)                     % record if converged before jmax          
            S_nc.k  =  S_nc.k+1;        % update k = k+1
            S_nc.E(S_nc.k,1) = wj;      % record eigval
            S_nc.err(S_nc.k,1) = ej;    % record step size
            S_nc.nj(S_nc.k,1) = jj;     % record number of iterations
            if(type==2 || type==3 || type==4)
                S_nc.V(:,S_nc.k) = vj;  % record eigvec
            end 
        else
            S_nc=zero(S_nc);
        end
    end % for kk
end % function 

function [wj,ej,jj]= Newton1(funA,fundA,breakN,jmax,w0)
    %-------------------------------------------------------------
    %- Simple Newton-Raphson Iteration
    %- wj1: previous eigenvalue, ej1: previous Newton step size 
    %- wj : current eigenvalue,  ej : current Newton step size 
    %-------------------------------------------------------------
    wj = w0;                            % initialize
    ej = 0;                             % initialize small
    for jj=1:jmax        
        wj1 = wj; 
        ej1 = ej; 
        step = 1/trace(funA(wj1)\fundA(wj1));                    
        if(isnan(step)||isinf(step));   break; end
        wj  = wj1 - step;         
        ej  = abs(wj - wj1);
        if(breakN(ej,ej1));        break; end
    end
end % function 

function [wj,vj,ej,jj]=Newton2(funA,fundA,breakN,jmax,w0,v0)
    %-------------------------------------------------------------
    %- Nonlinear Inverse Iteration 
    %- wj1: previous eigenvalue, ej1: previous Newton step size 
    %- wj : current eigenvalue,  ej : current Newton step size 
    %-------------------------------------------------------------
    wj = w0;                            % initialize
    vj = v0;                            % initialize
    ej = 0;                             % initialize small
    for jj=1:jmax        
        eHj = (vj')/norm(vj,2)^2;       % normalization vector
        wj1 = wj; 
        xj1 = vj; 
        ej1 = ej; 
        uj  = funA(wj1)\fundA(wj1)*xj1; % unnormalized v
        step = 1./(eHj*uj); 
        if(isnan(step)||isinf(step));   break; end
        wj  = wj1 - step;
        ej  = abs(wj - wj1);
        if(breakN(ej,ej1));        break; end
    end
end % function 