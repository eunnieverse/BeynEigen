function [UpdateA0,UpdateA1,UpdateA0_h,UpdateA1_h]=updateBeyn(funA,M_add,N,rmw,g,dg) 
    for j=1:N; %should add the change for all quadrature points
        invAj = funA(g(j))\M_add;
        UpdateA0 = 
        UpdateA0 = UpdateA0 + invAj * rmw(j) * dg(j); 
        UpdateA1 = UpdateA1 + invAj * rmw(j) * g(j) * dg(j); 
        if(j==N/2)
            UpdateA0_h = UpdateA0/(N*1i/2); 
            UpdateA1_h = UpdateA1/(N*1i/2); 
        end
    end
    UpdateA0 = UpdateA0 /(N*1i);
    UpdateA1 = UpdateA1 /(N*1i);
end%%function 