classdef EigenPairs
    %- Eigenpair Data Set Class
    properties (SetAccess=private)
        type        % int,              eigenpair type, used for plotting
                    %                   1: Beyn,   unconverged 
                    %                   2: Beyn,   converged  
                    %                   3: Newton, converged  
                    %                   4: Newton, final
                    %                   5: original answer    
    end
    properties
        k           % int,              number of eigenvalues
        E           % cdouble(k,1)      eigenvalue list 
        V           % cdouble(n,k)      eigenvector list
    
        err         % double(k,1)       Beyn: N/2 error, Newton: step size
        nj          % Newton iteration number, if applicable
    end
    methods
        %- Constructor
        function obj=EigenPairs(type)
            if(type==1||type==2||type==3||type==4||type==5)
                obj.type=type;
                obj.k=0;
            else error('EigenPairs type shoud be between 1..5');
            end
        end
        
        %- Utilities: zero, update, converged, error, copy
        function obj=zero(obj)
            %- re-initialize an eigenpair 
            obj.k=0; obj.E=[]; obj.V=[]; obj.err=[]; obj.nj=[];            
        end
        function obj=update(obj,k,E,V)
            %- set an eigenpair with given k, E, V
            obj.k=k;
            if(size(E,1)==k && size(E,2)==1 && size(V,2)==k) 
                obj.E=E;
                obj.V=V;
            else
                error('size of E and V does not match size of k'); 
            end
        end
        function obj=sample(obj,klist)
            %- Sample original list using index list klist
            obj.k=length(klist);
            obj.E=obj.E(klist);
            if(size(obj.V,2)>=obj.k)
                obj.V=obj.V(:,klist);
            end
            if(size(obj.err,1)>=obj.k)
                obj.err=obj.err(klist); 
            end
            if(size(obj.nj,1)>=obj.k)
                obj.nj=obj.nj(klist); 
            end
        end 
        function obj=error(obj,S)
            if(isa(S,'EigenPairs'))
                %- get error between obj.E and S.E 
                obj.err=zeros(obj.k,1); 
                for kk=1:obj.k
                    obj.err(kk) = min(abs(obj.E(kk)-S.E)); 
                end
            else
                error('both arguments should be EigenPairs'); 
            end                
        end
        function obj=copy(obj,S)
            %- copies two eigenpairs without changing the type 
            if(isa(S,'EigenPairs'))
                obj.k=S.k;
                obj.E=S.E;
                obj.V=S.V;
                obj.err=S.err;                 
            else
                error('both arguments should be EigenPairs'); 
            end           
        end
        function plot(obj)
            hold on; 
            switch obj.type
                case 1; scatter(real(obj.E),imag(obj.E),40,'r');
                case 2; scatter(real(obj.E),imag(obj.E),70,'b');
                case 3; scatter(real(obj.E),imag(obj.E),50,'g*');
                case 4; scatter(real(obj.E),imag(obj.E),50,'g.');
                case 5; scatter(real(obj.E),imag(obj.E),40,'b*'); 
            end
        end
    end % end methods 
  
    methods (Static) 
    end % end methods 
end % end class 