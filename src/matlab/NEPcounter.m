classdef NEPcounter
    properties (Constant)
        header=sprintf('%12s    %12s    %12s    %12s \n','time elapsed','nsolves','gflops','error')
        formatSpec='%12.5e    %12.5e    %12.5e    %12.5e \n'        
    end
    properties
        nmax                    %int            length of ErrList
        nn                      %int            current index 
        nsolves                 %int            # of solves
        t0                      %double         initial cputime 
        ErrList                 %double(nn,4)   [elapsed, nsolves, gflops, error]        
        fid
    end
    properties (Dependent)
        gflops                  %double         # of giga flops
        elapsed                 %double         recorded cputime         
    end
    methods
        %- Constructor
        function obj = NEPcounter()
            obj.nmax = 500;     %initialize
            obj.nn=0; 
            obj.t0 = cputime; 
            obj.ErrList=zeros(obj.nmax,4); 
        end
        
        %- Get functions 
        function elapsed=get.elapsed(obj)
            elapsed=cputime-obj.t0; 
        end      
        function gflops=get.gflops(obj)
                    % gflops = ns.^2 ./ t0 * 1e-9
            gflops= (obj.nsolves)^2./obj.elapsed * 1e-9; 
        end
        
        %- Utilities
        function obj=add(obj,nsolves,error)
            if(obj.nn >= obj.nmax)
                obj.ErrList = [obj.ErrList; zeros(obj.nmax,4)]; 
                disp(sprintf('NEPcounter nmax=%d exceeded, changed to %d',obj.nmax, 2 * obj.nmax )); 
                obj.nmax = 2 * obj.nmax; 
            end
            obj.nn = obj.nn + 1; 
            obj.nsolves=nsolves; 
            obj.ErrList(obj.nn,:) = ...
                [obj.elapsed, obj.nsolves, obj.gflops, error]; 
        end

        function log(obj)    
            obj.fid = fopen('counter.log','w');
            fprintf(obj.fid,obj.header);
            for ii = 1:obj.nn
                fprintf(obj.fid, obj.formatSpec, obj.ErrList(ii,:));
            end
            fclose(obj.fid);
        end
        
        function plot(obj,type)
            switch type
                case 1 
                    x = obj.ErrList(:,1); %obj.elapsed
                    xlab = 'time elapsed (s)';
                case 2
                    x = obj.ErrList(:,2); %obj.nsolves                     
                    xlab = 'number of solves';
                case 3 
                    x = obj.ErrList(:,3); %obj.gflops
                    xlab = 'gflops (n^2/t)';
            end
            y = obj.ErrList(:,4); 
            
            semilogy(x,y,'ro');
            xlabel(xlab);
            ylabel('maximum error'); 
        end
        
    end % methods
end