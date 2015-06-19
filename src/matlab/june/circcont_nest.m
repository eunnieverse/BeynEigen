function [gamma,gammap] = circcont_nest(g0,rho,N)
%---------------------------------------------------------------------
%%% function circcont_nest: compute nested circular contour given N 
%---------------------------------------------------------------------
    theta=[]; 
    for qq=1:log2(N) 
        for rr = 1:2:2^qq
            theta = [theta 2*pi*rr/(2^qq)]; 
        end
    end   
    gamma = g0 + rho*(cos(theta) + 1i*sin(theta));
    gammap = rho*(-sin(theta) + 1i*cos(theta)); 
end