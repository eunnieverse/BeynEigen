%% housekeeping
clear all; 
close all;
%% create or load funA and newtA
filebase = 'poly2_100'; 
load(strcat(filebase,'_fun')); 
load(strcat(filebase,'_E')); 
g0 = 0; 
rho = 0.5; 
N = 150; 
[gamma,gammap] = circcont(g0,rho,N);

%% Plot simple BeynEigen Run 
xLc = [real(g0)-rho*1.5 real(g0)+rho*1.5]; %%xL for contour 
yLc = [imag(g0)-rho*1.5 imag(g0)+rho*1.5]; %%yL for contour
xL = [-max(abs(E))*1.2 max(abs(E))*1.2]; %%xL for whole problem
yL = [-max(abs(E))*1.2 max(abs(E))*1.2]; %%yL for whole problem 

cfig = figure();
    scatter(real(gamma),imag(gamma),30,'.'); %% contour  
    hold on; 
    scatter(real(E),imag(E),100,'b*'); %%answer    
    xlim(xL); ylim(yL); 
    line([0 0], xL,'Color','k','Linewidth',1.5);
    line(yL, [0 0],'Color','k','Linewidth',1.5);
    hold off; 
    axis square; 
    xlabel('Re(w)');ylabel('Im(w)');