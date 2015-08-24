function S_nc = Newton(fA, type, S_bc, S_nc)
    %- run Newton-based Iteration on S_bc 
    %- convergence criterion:
    %       stop when error increases or decreases slowly: if(ej>0.9*ej1) 
    %  type 1: simple Newton-Raphson
    %       2: Inverse iteration using eigenvector estimate
    %       3: Residual inverse iteration 
    %-------------------------------------------------------------
    %- wj, wj1: current & previous eigval
    %- ej, ej1: current & previous error ej=abs(wj-wj1)        
    %-------------------------------------------------------------
    nnmax = 20; % maximum iteration       
    
    switch type 
        case 1 %  type 1: simple Newton-Raphson    
            for jj=1:S_bc.k         % for every eigenpair
                wj = S_bc.E(jj);    % set initial values for wj ej ej1 nj
                ej = 0;          
                ej1= Inf;            
                nj = 0;            
                while(ej<ej1) 
                    if(nj==nnmax)
                        break;      
                    end
                    nj  = nj+1; 
                    wj1 = wj;
                    ej1 = ej;
                    wj  = wj1 - 1/trace(fA.funA(wj1)\fA.fundA(wj1)); 
                    ej  = abs(wj - wj1);
                end % while
                if(nj<nnmax)                    % if converged
                    S_nc.k = S_nc.k+1;          % k = k+1
                    S_nc.E = [S_nc.E; wj];      % record E(k)
                    S_nc.err = [S_nc.err; ej];  % record err(k)                    
                  % S_nc.V(:,S_nc.k) = S_bc.V(:,jj); 
                end
            end %for jj 
        %-----------------------------------------------------------------           
        case 2 %       2: Nonlinear Inverse Iteration, use eigenvector info
            for jj=1:S_bc.k         % for every eigenpair
                wj = S_bc.E(jj);    % set initial values               
                vj = S_bc.V(:,jj);  
                ej = 0;          
                ej1= Inf;            
                nj = 0;                            
                while(ej<ej1) 
                    if(nj==nnmax)
                        break;      
                    end
                    nj  = nj+1; 
                    eHj = conj(vj')/norm(vj,2)^2; % normalization vector
                    wj1 = wj;
                    xj1 = vj; 
                    ej1 = ej;
                    uj  = fA.funA(wj1)\fA.fundA(wj1)*xj1; %unnormalized v
                    wj  = wj1 - 1./eHj*uj;
                    vj  = uj.* (1./eHj*uj); 
                    ej  = abs(wj-wj1); 
                end % while                 
                if(nj<nnmax)                    % if converged
                    S_nc.k = S_nc.k+1;          % k = k+1
                    S_nc.E = [S_nc.E; wj];      % record E(k)
                    S_nc.V = [S_nc.V vj]; 
                    S_nc.err = [S_nc.err; ej];  % record err(k)                    
                end                
            end %for jj 
        %-----------------------------------------------------------------            
        case 3 %       3: Nonlinear Residual Inverse Iteration 
            % Not used 
        case 4 %       4: Block Newton Algorithm 
            % Not used 
    end % switch 
end % function 