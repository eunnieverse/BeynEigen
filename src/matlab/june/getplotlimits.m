function [xLc,yLc] = getplotlimits(gamma)
    minx = min(real(gamma)); 
    maxx = max(real(gamma)); 
    miny = min(imag(gamma)); 
    maxy = max(imag(gamma)); 
    xLc = [minx-0.5*(maxx-minx) maxx+0.5*(maxx-minx)]; %%xL for contour 
    yLc = [miny-0.5*(maxy-miny) maxy+0.5*(maxy-miny)]; %%yL for contour
end