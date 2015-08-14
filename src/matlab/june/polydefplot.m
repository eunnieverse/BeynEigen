function cfig=polydefplot(filebase,gamma,N)
    %% create or load funA and newtA
    load(strcat(filebase,'_fun')); 
    load(strcat(filebase,'_E')); 
    
    %% Plot simple BeynEigen Run 
    [xLc,yLc] = getplotlimits(gamma); 
    xL = [-max(abs(E))*1.2 max(abs(E))*1.2]; %%xL for whole problem
    yL = [-max(abs(E))*1.2 max(abs(E))*1.2]; %%yL for whole problem 

    cfig = figure();
        scatter(real(E),imag(E),100,'b*'); %%answer    
        scatterContour(gamma,N); 
        xlim(xLc); ylim(yLc); axis square; 
        line([0 0], xLc,'Color','k','Linewidth',1.5);
        line(yLc, [0 0],'Color','k','Linewidth',1.5);
        xlabel('Re(w)');ylabel('Im(w)');
        title(' ');
end

function scatterContour(g,N)
    hold on; 
    scatter(real(g(1:N)),imag(g(1:N)),30,'.'); %% contour 
end
