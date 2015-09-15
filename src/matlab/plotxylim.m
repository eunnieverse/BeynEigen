function [xLc,yLc]=plotxylim(g)
    %- getplotlimits: get plot limits([xL,yL]); 
    minx = min(real(g)); 
    maxx = max(real(g)); 
    miny = min(imag(g)); 
    maxy = max(imag(g)); 
    xLc = [minx-0.5*(maxx-minx) maxx+0.5*(maxx-minx)];  
    yLc = [miny-0.5*(maxy-miny) maxy+0.5*(maxy-miny)];             
end
