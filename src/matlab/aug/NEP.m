classdef NEP
    %- Nonlinear Eigenvalue Problem (NEP) Class
    properties (Constant)
        format1='poly%d_%d'         % string,       format, type 1 poly
        format2='EP_%d'             % string,       format, type 2 linear
    end
    properties (SetAccess=private)
        type                        % int,          eigenpair type 
                                    %               1: polynomial NEP
                                    %               2: linear EP 
    end
    properties
        n                           % int,          problem size A(n,n)
        filebase                    % string,       name of the file 
        funA                        % function definition 
        fundA                       % function derivative definition
        isHermitian                 % int,          1 if Hermitian, 0 if not
        S                           % Solution EigenPairs, contains k,E,V 
                                    % k: int,         # of eigenvalues 
                                    % E: double[k,1]  solution eigenvalues
                                    % V: double[n,k]  solution eigenvectors
    end
    %-------------------------------------------------------------
    %-------------------------------------------------------------
    %-------------------------------------------------------------
    methods
        %- Constructor 
        function obj=NEP(type)
          switch type 
            case 1; disp('Initialized NEP, 1: poly   '); 
            case 2; disp('Initialized NEP, 2: linear '); 
            otherwise; error('EigenPairs type shoud be between 1..2');
          end % switch
            obj.type=type; 
            obj.S=EigenPairs(5); % ask Homer where to declare type of obj.S 
            obj.isHermitian=0;   % default             
        end
        function obj=NEP_load(obj,filebase)
            %- fill NEP data from an existing file 
            obj.filebase=filebase; 
            if(exist(strcat(filebase,'_E.mat'))==2)
                m=load(strcat(obj.filebase,'_fun.mat')); 
                obj.n=m.n; 
                obj.funA=m.funA;
                obj.fundA=m.fundA;
                m=load(strcat(obj.filebase,'_E.mat')); 
                obj.S=update(obj.S,length(m.E),m.E,m.V);
            else
                error('file does not exist: create an NEP instead');
            end
        end        
        %- Utils
        function obj=NEP_poly(obj,n,p)
            %- Create a polynomial eigenproblem 
            if(obj.type~=1)
                error('NEP type should be 1 (polynomial eigenproblem)');
            end
            obj.n=n;
            obj.filebase=sprintf(obj.format1,p,n);
            % Generate A0, A1, A2: random matrices
            for jj=0:p
                 eval(sprintf('A%d=rand(%d); coeffs{%d+1}=A%d;',jj,n,jj,jj));
                 %%% ex) A2 = rand(n), coeffs{2}=A2; 
            end
            % define anonymous fcns. for funA=A(w), fundA=dA/dw
            s0 = 'A0';                  % string for 'A0,A1,A2,...'
            s1 = 'funA = @(w) A0';      % string for defining funA 
            s2 = 'fundA = @(w)';        % string for defining fundA
            for kk = 1:p %%polynomial order 
                s1 = sprintf('%s + A%d*w^(%d)',s1,kk,kk); 
                s2 = sprintf('%s + %d*A%d*w^(%d)',s2,kk,kk,kk-1); 
                s0 = sprintf('%s,A%d',s0,kk);
            end
            eval(strcat(s1,';')); %% funA = @(w) A0 + A1*w^1 + A2*w^2; 
            eval(strcat(s2,';')); %% fundA = @(w) A0 + 1*A1 + 2*A2*w^1; 

            % Run polyeig: [V,E,Econd] = polyeig(A0,A1,A2)
            eval(strcat('[V,E,Econd] = polyeig(',s0,');')); 
            obj.funA= funA;
            obj.fundA=fundA;
            k=length(E);
            obj.S=update(obj.S,k,E,V); 
        end
        
        function obj=NEP_linear(obj,n)
            %- Construct a linear eigenproblem from a given random matrix
            if(obj.type~=2)
                error('NEP type should be 1 (polynomial eigenproblem)');
            end
            obj.n=n; 
            obj.filebase=sprintf(obj.format2,n);
            E = rand(n,1)*1-rand(n,1)*1+rand(n,1)*1i-rand(n,1)*1i;
            obj.funA = @(z) diag(E) - z*eye(n); 
            obj.fundA= @(z) -eye(n); 
            obj.S=update(obj.S,length(E),E,eye(n)); 
        end
        
        function save(obj)
            % generate mfiles to store values
            n=obj.n;
            funA=obj.funA;
            fundA=obj.fundA;
            k=obj.S.k;
            E=obj.S.E;
            V=obj.S.V; 
            save(sprintf('%s_fun.mat',obj.filebase),'n','funA','fundA'); 
            save(sprintf('%s_E.mat',obj.filebase),'n','V','E');
        end
        
        function cfig=plot(obj)
            cfig = figure();
            set(cfig,'Position',[80   540   560   420]); 
            scatter(real(obj.S.E),imag(obj.S.E),100,'b*');   
            xlabel('Re(w)');
            ylabel('Im(w)');
        end
    end % methods 
end