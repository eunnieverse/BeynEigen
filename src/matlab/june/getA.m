function [BeynA] = getA(p, funA, M, g, gp)
    f_BeynA = @(z) (z^p)*funA(z)\M;
    BeynA = cint(f_BeynA, g, gp);
    
function [BeynA] = getAmod(p, funA, w_Newt, M, g, gp)
    f_BeynA = @(z) (z^p)*fz(z,w_Newt)*funA(z)\M; 
    BeynA = cint(f_BeynA, g, gp);
    
function f = fz(z,w_Newt) %(z-w1)(z-w2)... 
    f=1; 
    if(length(w_Newt)>0)
        for(j=1:length(w_Newt)) 
            f=f*(z-w_Newt(j)); 
        end
    end
    
function [] = BeynSVD(BeynA0,BeynA1)
    [V,Sigma,W]=svd(BeynA0);
    Sigma= Sigma(1:l,1:l); 
    s = diag(Sigma); 
    k = sum(s>1e-15); %%actual rank 
    
    V0 = V(1:n, 1:k);