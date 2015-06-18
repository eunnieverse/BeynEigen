function [gamma,gammap] = circcont(g0,rho,N)
    %---------------------------------------------------------------------
    %%% function circcontour: compute circular contour gamma 
    %%%                 given the center g0, radius rho, numpts N
    %%%                 return contour gamma 
    %---------------------------------------------------------------------
    theta = linspace(0,2*pi,N+1); 
    gamma = g0 + rho*(cos(theta) + 1i*sin(theta));
    gammap = rho*(-sin(theta) + 1i*cos(theta)); 
end