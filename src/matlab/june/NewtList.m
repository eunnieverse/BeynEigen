function wjlist= NewtList(funA,fundA,w0list,nn)
%---------------------------------------------------------------------
%%% Newt.m : iterate to update v using newtA(v)=v-det(A)/det'(A)
%%%          omegaIn: a list of initial guesses 
%%%          nn: iteration number
%---------------------------------------------------------------------
    %%% define iteration function and condition number 
    newtA = @(wn) wn - 1/trace(funA(wn)\fundA(wn));
    condA = @(wn) cond(funA(wn)\fundA(wn));
    %rcondA =@(wn) rcond(funA(wn)\fundA(wn));
    wjlist = w0list ; %%data size is k x 1 
    for jj = 1:length(w0list); %%for every eigenvalue 
        for ii = 1:nn %% repeat for specified number of iterations nn
            cond=condA(wjlist(jj)); %%condition number for A\(dA/dw)
            if(cond<1e-10)
                printf('Newton iteration completed at cond = %f \n',cond); 
                break;
            else
                wjlist(jj) = newtA(wjlist(jj));    %%iterate
            end
        end
    end
end
