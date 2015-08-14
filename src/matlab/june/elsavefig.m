function elsavefig(cfig, savefigbase, savefigindex, savejpg, saveeps)
    savefigname=sprintf('%s_%03d',savefigbase,savefigindex); 
    if(savejpg==1) saveas(cfig, strcat(savefigname,'.jpg')); end; 
    if(saveeps==1) saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 
end