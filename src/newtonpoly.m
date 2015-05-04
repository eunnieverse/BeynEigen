%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% BeynEigen
%%%%% Yoonkyung Eunnie Lee 
%%%%% matlab program to find the eigenvalue of a nonlinear eigenproblem
%%%%% using Beyn's contour integral method 
%%%%% 2015.05.04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% newtonpoly.m : iterate v->v-det(A)/det'(A)
%%%            given a size mm polynomial eigenfunction A(w),
%%%            initial guess w0, 
%%%            and the iteration number nn, 
%%%            compute the update of w based on newton iteration 
%%%            wn = w0; 
%%%            wn = wn - 1/trace(inv(A)*dAdw): iterate nn times 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [wn, wnlist] = newtonpoly(w0,nn,varargin) %%(w0,nn,A0,A1,A2) 
    p = nargin-3; %% max polynomial order 
    n = length(varargin{1}); 
    wn = w0; %% initial value as given 
    wnlist = zeros(nn+1,1);
    wnlist(1)=w0; 
    for ii = 1:nn %%for each newton iteration
        matA = varargin{1}; %% A = A0 
        matdA = zeros(n); %% dA/dw = 0           
        for k = 1:p %%polynomial order 
            matA = matA + varargin{k+1}*wn^(k); %% A = A0 + A1*w0 + .. 
            matdA = matdA + varargin{k+1}*k*wn^(k-1); %%dA/dw=0+1*A1+2*A2*w0
        end
        invA = inv(matA);
        wn = wn - 1/trace(invA*matdA); %%% update wn
        wnlist(ii+1) = wn; 
    end  
end