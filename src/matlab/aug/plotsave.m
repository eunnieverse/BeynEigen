function plotsave(cfig, savefigbase, fignum, savejpg, saveeps)
    savefigname=sprintf('%s_%03d',savefigbase,fignum); 
    if(savejpg==1); saveas(cfig, strcat(savefigname,'.jpg'));       end; 
    if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end;
end
