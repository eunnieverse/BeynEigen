classdef BeynData
    %- Data for Beyn run: contour g(N) dg(N), random matrix M(n,l)
    properties (Constant)
        tolw=1e-15;                 % smallest number accepted as eigval
    end
    properties
        
        %- random matrix M 
        n                           % int,          problem size A(n,n)
        l                           % int,          M column size 
        Mmax                        % double(n,n)   full size random matrix
        %- contour g 
        N                           % int,          N of g 
        Nmax                        % int,          max length of g, dg
        gmax                        % cdouble[Nmax] whole contour
        dgmax                       % cdouble[Nmax] whole contour
        emax                        % double,       Beyn cutoff (max error)
        %- previously converged Newton eigenvalues 
        k                           % number of rmw eigenvalues 
        E_nc                        % rmw eigenvalues 
        %- 
        NA                          % int,          N of BeynA0 and BeynA1
        BeynA0sum                   % BeynA0 summed at NA points on g
        BeynA1sum                   % BeynA1 summed at NA points on g
    end
    properties (Dependent)
        M                       % cdouble(n,l)      Mmax(:,1:l) 
        g                       % cdouble[N]  list of N contour points
        dg                      % cdouble[N]  list of N dg/dw        
    end
    %-------------------------------------------------------------
    %-------------------------------------------------------------
    %-------------------------------------------------------------      
    methods
        %- Constructor
        function obj=BeynData(N,Mmax,l,gmax,dgmax,emax)
            obj.N = N; 
            obj.l = l; 
            obj.n = size(Mmax,1);
            obj.Mmax = Mmax; 
            obj.Nmax = length(gmax); 
            obj.gmax = gmax; 
            obj.dgmax= dgmax;
            obj.emax = emax; 
            obj.NA = 0;
            obj.k = 0; 
            obj.E_nc = []; 
        end
        %-------------------------------------------------------------
        %- Get functions 
        function M = get.M(obj)
            %- Get random matrix M(n,l), sampled from Mmax(n,n)
            M= obj.Mmax(:,1:obj.l);
        end        
        function g=get.g(obj)
            %- Get contour points g(N), sampled from gmax(Nmax)
            g=obj.gmax(1:obj.N);
        end
        function dg=get.dg(obj)
            %- Get contour differential dg(N), sampled from dgmax(Nmax)
            dg=obj.dgmax(1:obj.N);
        end
        %-------------------------------------------------------------
        %- Set functions 
        %-------------------------------------------------------------
        %- Utils
        function obj= halfBeynA(obj,funA,rmw)
            %- Compute summation of BeynA0, BeynA1 for half of N 
            obj.NA=obj.N/2; % sync NA = N/2 (half data pts) 
            obj.BeynA0sum = zeros(obj.n,obj.l); 
            obj.BeynA1sum = zeros(obj.n,obj.l); 
            for ii=1:obj.NA; 
                invA = funA(obj.g(ii))\obj.M;
                obj.BeynA0sum = obj.BeynA0sum  + ...
                    invA * rmw(ii) * obj.dg(ii);
                obj.BeynA1sum = obj.BeynA1sum  + ...
                    invA * rmw(ii) * obj.dg(ii) * obj.g(ii) ;
            end 
        end
        
        function obj= fullBeynA(obj,funA,rmw)
            %- update summation of BeynA0, BeynA1 for whole of N, 
            %- using halfBeynA or data from previous run 
            if(obj.NA ~= obj.N/2)
                error('run halfBeynA before proceeding');
            end
            obj.NA = 2 * obj.NA ; % double obj.NA 
            for ii=(obj.NA/2+1):obj.NA
                invA = funA(obj.g(ii))\obj.M;
                obj.BeynA0sum = obj.BeynA0sum + ...
                    invA * rmw(ii) * obj.dg(ii); 
                obj.BeynA1sum = obj.BeynA1sum + ...
                    invA * rmw(ii) * obj.dg(ii) * obj.g(ii); 
            end
        end%%function         

        function plot(obj) 
            %- override plot function for contour
            [xLc,yLc]=plotxylim(obj.g); 
            hold on; 
            scatter(real(obj.g(1:obj.N)),imag(obj.g(1:obj.N)),60,'.');
            xlim(xLc); 
            ylim(yLc); 
            line([0 0], xLc,'Color','k','Linewidth',1.5);
            line(yLc, [0 0],'Color','k','Linewidth',1.5);
            axis square; 
        end
    end % methods 
end