clear,clc;
n = 100;        % Number of columns
m = 400; % Number of rows
s = 10; % Number of nonzeros in solution
maxiter = 16000; % Number of iterations 
num_repeats = 20;
lambda = 1;
exp2_unblock_undelay(n,m,s,lambda,maxiter,num_repeats);