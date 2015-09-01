function [gmax,dgmax,s,dc,isinside]=NestedContour(g0,rho,Nmax) 
    %- create a circular nested contour centered at g0, with radius rho 
    theta    = zeros(Nmax,1);
    theta(1) = 2*pi;  % qq = 0; 
    for qq = 1:log2(Nmax)
        for jj = 1:2^(qq-1) 
            theta(2^(qq-1)+jj) = 2*pi*jj/2^(qq-1);
        end
    end
    gmax  = g0 + rho*( cos(theta) + 1i*sin(theta));
    dgmax =      rho*(-sin(theta) + 1i*cos(theta));
    %----------------------------------------------------------------
    s = @(N) 2*pi*rho/N;        %- arc distance between adjacent pts
                                %  as a function of N   
    %-----------------------------------------------------------------                            
    dc = @subfun1; 
    function y=subfun1(w)
        %- anaytical distance from w to the circular contour
        y=zeros(size(w));
        y=rho-abs(w-g0);
    end
    %-----------------------------------------------------------------
    isinside = @subfun2; 
    function y=subfun2(w)
        y=zeros(size(w)); 
        y=abs(w-g0)<rho;
    end
    %-----------------------------------------------------------------                            
end