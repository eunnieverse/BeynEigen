%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% function polymat: computes the value of A(w0) and dA(w0)/dw 
%%%                    at a given w0, where A(w)=sum(Ai*w^i) is the
%%%                    polynomial function. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [matA,matdA]=polymat(w0,varargin)
    p = nargin-2; %% max polynomial order 
    n = length(varargin{2});
    matA = varargin{1}; %% A = A0 
    matdA = zeros(n); %% dA/dw = 0 
    for k = 1:p
        matA = matA + varargin{k+1}*w0^(k); %% A = A0 + A1*w0 + .. 
        matdA = matdA + varargin{k+1}*k*w0^(k-1); %%dA/dw=0+1*A1+2*A2*w0
    end
end