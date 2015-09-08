close all; clear all; 
addpath('/home/eunnie12/Work/BeynEigen/src/matlab'); 
savefigbase = 'contourdistance_n100_Mrand';
load(strcat(savefigbase,'.mat')); 

clist='krgbcmkrgbcmkrgbcm';
showplot=1; 
saveeps=1; savejpg=1;       % choose conditions
savemov=0;                  % save movie 
fignum=1;                   % initialize figure number


[B,I]=sort(Slist(length(Slist)).err,'descend');
Eerr = Slist(length(Slist)).E(I(1));                    % eigval with maximum error

%---------------------------------------------------------------------
%- Plot error 
%---------------------------------------------------------------------
cfig = figure(); 
for gg=1:size(x,2)
    loglog (x(:,gg),y(:,gg),'Color',clist(gg),'LineWidth',2,'Marker','o'); 
    hold on;     
    legendlist{gg}=sprintf('N=%4d',2^(gg+ng0)); 
end
legend(legendlist); 
xlabel('d/s'); 
ylabel('maximum absolute error'); 
ylim([min(y(:))*0.1 max(y(:))*10]);
xlim([min(x(:))*0.1 max(x(:))*10]);

plotsave(cfig, savefigbase, fignum, savejpg, saveeps);
fignum = fignum +1; 
%---------------------------------------------------------------------
%- Plot solutions 
%---------------------------------------------------------------------
cfig=plot(fA); hold on; plot(Slist(length(Slist)));  
plot(real(wlist(:,end)),imag(wlist(:,end)),'r-x','LineWidth',2);
scatter(real(Eerr),imag(Eerr),100,'g*');
h=legend('solution','computed', 'moved \omega', '\omega w. max error'); 
set(h,'Location','northeastoutside'); 
plot(BDlist(length(BDlist)));
%title(sprintf(['N=%5d'],BDlist(length(BDlist)).N)); 

plotsave(cfig, savefigbase, fignum, savejpg, saveeps);
fignum = fignum +1; 