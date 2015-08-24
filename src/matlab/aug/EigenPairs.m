classdef EigenPairs
    %- Eigenpair Data Set Class
    properties (SetAccess=private)
        type        % int,              eigenpair type 
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
    
        err         % double(k,1)       error list 
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
        
        %- Utilities
        function obj=update(obj,k,E,V)
            obj.k = 0; 
            obj.E = []; 
            obj.V = []; 
            obj.k=k;
            obj.E=E;
            obj.V=V;
        end        
        function obj=converged(obj,listconv)
            obj.k=length(listconv);
            obj.E=obj.E(listconv);
            obj.V=obj.V(:,listconv);
            obj.err=obj.err(listconv); 
        end 
        function obj=error(obj,S)
            if(isa(S,'EigenPairs'))
                %- get error between obj.E and S.E 
                obj.err=zeros(obj.k,1); 
                for jj=1:obj.k
                    obj.err(jj) = min(abs(obj.E(jj)-S.E)); 
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
                case 1
                    scatter(real(obj.E),imag(obj.E),40,'r');
                case 2
                    scatter(real(obj.E),imag(obj.E),70,'b');
                case 3
                    scatter(real(obj.E),imag(obj.E),50,'g*');
                case 4
                    scatter(real(obj.E),imag(obj.E),50,'g.');
            end
        end
    end % end methods 
  
%     methods (Static) 
%         function err= geterr(S1,S2)
%             if(isa(S1,'EigenPairs') && isa(S2,'EigenPairs'))
%                 %- get error between obj.E and S.E 
%                 err=zeros(S1.k,1); 
%                 for jj=1:S1.k
%                     err(jj) = min(abs(S1.E(jj)-S2.E(jj))); 
%                 end
%             else
%                 error('both arguments should be EigenPairs'); 
%             end                
%         end
%     end % end methods 
end % end class 