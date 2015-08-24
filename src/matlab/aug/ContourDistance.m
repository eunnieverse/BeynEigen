%--------------------------------------------------------------------------
%%%%% ContourDistance.m 
%%%%% Yoonkyung Eunnie Lee 
%%%%% created       2015.06.23
%%%%% last modifled 2015.08.12
%--------------------------------------------------------------------------
%%% randomly generated matrix A is created, saved, and plotted here 
%%% a_rand: matrix A 
%%% v_rand: eigenvectors (Av=wv) 
%%% w_rand: eigenvalues, where w_rand(30) is tuned from 0.3 to 0.5 
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% %% create the file 
% n_rand = 30;   % size of matrix A
% w_rand = rand(n_rand,1)*1-rand(n_rand,1)*1+rand(n_rand,1)*i-rand(n_rand,1)*i;
% a_rand = diag(w_rand); % matrix A 
% v_rand = eye(n_rand);  % eigenvectors 
% w_rand(30)=0.3; 
% save('ContourDistance.mat','n_rand','w_rand','a_rand','v_rand'); 
%--------------------------------------------------------------------------

% %% check 
% [V,E]=eig(a_rand); 
% E=diag(E); 
% scatter(real(E),imag(E)); 

%--------------------------------------------------------------------------
%% load file
load('ContourDistance.mat'); 
[g,dg]=circcont_nest(0,0.5,100); 
% each column of w_matrix is the list of eigenvalues
% the last eigenvalue moves towards 0.5 
w_matrix = [w_rand w_rand w_rand w_rand w_rand]; 
w_matrix(30,:)=0.3:0.05:0.5; 
save('ContourDistance.mat','n_rand','w_rand','a_rand','v_rand','w_matrix'); 
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
%% plot 
savefig = 1;
saveeps = 1; 
extension = '.png';
cfig =figure(); 
scatter(real(w_rand),imag(w_rand),'kx');  hold on; 
scatter(real(g),imag(g),'b.'); axis equal;
plot(real(w_matrix(30,:)),imag(w_matrix(30,:)),'r.-'); 
h=legend('all eigvals','contour','moved eigval','location','northeastoutside');
savefigname='ContourDistance';
    if(savefig==1); saveas(cfig, strcat(savefigname,extension)); end; 
    if(saveeps==1); saveas(cfig,strcat(savefigname,'.eps'),'epsc'); end; 
%--------------------------------------------------------------------------