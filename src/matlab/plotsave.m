function plotsave(cfig, savefigbase, fignum, savejpg, saveeps)
    set(cfig,'PaperPositionMode','auto'); 
    set(findobj('type', 'text'), 'color', 'black');
    set(findobj('type', 'line'), 'linewidth', 2.0);
    savefigname=sprintf('%s_%03d',savefigbase,fignum);      
    if(savejpg==1) 
        saveas(cfig, strcat(savefigname,'.jpg'));       
    end
    if(saveeps==1)
        print(cfig,strcat(savefigname,'.eps'),'-r600', '-depsc', '-tiff', '-cmyk'); 
    end
end
