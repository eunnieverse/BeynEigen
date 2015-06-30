function cfig=polydefplot(filebase,gamma)
    extension='.jpg';
    savefig=1; 
    saveeps=1; 
    
    %% create or load funA and newtA
    load(strcat(filebase,'_fun')); 
    load(strcat(filebase,'_E')); 
    
    %% Plot simple BeynEigen Run 
    minx = min(real(gamma)); 
    maxx = max(real(gamma)); 
    miny = min(imag(gamma)); 
    maxy = max(imag(gamma)); 
    
    xLc = [minx-0.5*(maxx-minx) maxx+0.5*(maxx-minx)]; %%xL for contour 
    yLc = [miny-0.5*(maxy-miny) maxy+0.5*(maxy-miny)]; %%yL for contour
    xL = [-max(abs(E))*1.2 max(abs(E))*1.2]; %%xL for whole problem
    yL = [-max(abs(E))*1.2 max(abs(E))*1.2]; %%yL for whole problem 

    cfig = figure();
        scatter(real(gamma),imag(gamma),30,'.'); %% contour  
        hold on; scatter(real(E),imag(E),100,'b*'); %%answer    
        xlim(xLc); ylim(yLc); 
        line([0 0], xLc,'Color','k','Linewidth',1.5);
        line(yLc, [0 0],'Color','k','Linewidth',1.5);
        hold off; 
        axis square; 
        xlabel('Re(w)');ylabel('Im(w)');
        savefigname=strcat(filebase,'_polydef');
        if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end;
        if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc');end;