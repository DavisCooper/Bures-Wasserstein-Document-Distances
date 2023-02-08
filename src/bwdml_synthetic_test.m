clear all
warning('off','MATLAB:sqrtm:SingularMatrix');
addpath('util/')

% define data parameters
d = 2;
n_c = 2;
m = 50;
alpha = 2;
beta = 2;
[y,mu,S] = synthetic_data(d,m,n_c,alpha,beta);

% define problem parameters
tol = 1e-6;
maxItrs = 5e3;
numConst = 1e2;
gamma = 1e3;
knn_size = 5;

M_0 = eye(d);
X_0 = eye(d);
    
profile on
[l,u] = ComputeBWDExtremes(mu,S,5,95,M_0,X_0);
C = GetConstraints(y,numConst,l,u);
[M_sol,X_sol,b_sol] = bwdml(M_0,X_0,mu,S,C,gamma,maxItrs,tol);
profile off
profile report
   
%%

figure(1)
plot_class_distances(y,mu,S,M_0,X_0);   
plot_class_distances(y,mu,S,M_sol,X_0);


