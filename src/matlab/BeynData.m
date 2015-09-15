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
        p                           % number of rmw eigenvalues 
        Ep                          % rmw eigenvalues 
        %- 
        NA                          % int,          N of BeynA0 and BeynA1
        BeynA0sum                   % BeynA0 summed at NA points on g
        BeynA1sum                   % BeynA1 summed at NA points on g
    end
    properties (Dependent)
        M                           % cdouble(n,l)      Mmax(:,1:l) 
        g                           % cdouble[N]  list of N contour points
        dg                          % cdouble[N]  list of N dg/dw        
        rmw                         % function handle, used to remove 
                                    % fully converged eigenvalues 
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
            obj.p  = 0; 
            obj.Ep = []; 
        end
        %-------------------------------------------------------------
        %- Get functions 
        %-------------------------------------------------------------
        function rmw = get.rmw(obj)                        
            if(obj.p==0)
                error('BenyData k should be larger than 0 to call rmw');
            end            
            rmw = @subfun; 
            function y = subfun(z)
                %- return y = (z-E(1))*(z-E(2))*...*(z-E(k)) 
                y=1; 
                for kk=1:obj.p
                    y=y*(z-obj.Ep(kk)); 
                end
            end % function definition    
        end
               
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
        
        %- Utils
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