n_rand = 30; 
w_rand = rand(n,1)*1-rand(n,1)*1+rand(n,1)*i-rand(n,1)*i;
a_rand = diag(w_rand); 
v_rand = eye(n); 
scatter(real(w_rand),imag(w_rand)); 

%[V,E]=eig(a_rand); 
%E=diag(E); 
%scatter(real(E),imag(E)); 

save('ContourDistance.mat','n_rand','w_rand','a_rand','v_rand'); 
