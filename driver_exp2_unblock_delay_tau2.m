clear,clc;
n = 100;        % Number of columns
m = 400; % Number of rows
s = 10; % Number of nonzeros in solution
maxiter = 12000; % Number of iterations 
num_repeats = 25;
lambda = 1;
alpha = 2;
exp2_unblock_delay_tau2(n,m,s,alpha,lambda,maxiter,num_repeats);