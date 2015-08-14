function scatterBeyn(w_Beyn,N)
    hold on; 
    scatter(real(w_Beyn),imag(w_Beyn),70,'b'); 
    title(sprintf('N=%d',N));
end
