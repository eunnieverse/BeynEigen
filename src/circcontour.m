%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% BeynEigen
%%%%% Yoonkyung Eunnie Lee 
%%%%% matlab program to find the eigenvalue of a nonlinear eigenproblem
%%%%% using Beyn's contour integral method 
%%%%% 2015.05.04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% main for Beyn's Method 
%function main()
%end %%function main

    clear all; 
    close all; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% run / load polyeigdef 
    matfilebase = 'polyeig20150504'; 
    % polyeigdef(matfilebase, 2, 6);  %% matfilebase, pp, mm 
    %%% load the mfile containing A0,A1,A2,...Ap,Alist, pp,mm
    load(strcat(matfilebase,'.mat'));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% define contour g, g', and integrand function fun
    
    g0 = 0.0; %center
    rho= 1.0; %radius 
    g = @(theta) g0 + rho* (cos(theta) + 1i*sin(theta));
    gprime = @(theta) rho* (-sin(theta)+ 1i*cos(theta));

    fun = @(z) exp(z)./z;     %% define an anonymous function for the integrand 



integrandf = @(t) fun(g(t)).*gprime(t); 

q1 = integral(@(t) fun(g(t)).*gprime(t),0,2*pi) %% compute using integral
[qgk, errgk] = quadgk(integrandf,0,2*pi) %% compute using quadgk 
