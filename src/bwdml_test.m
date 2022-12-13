clear all
warning('off','MATLAB:sqrtm:SingularMatrix');

% define grid parameter search space
tol = 1e-6;
maxItrs = 1e5;
numConst = 1e2;
gamma = 1e3;
knn_size = 5;

% load datasets
load('../data/bbcsport_cov_tr_te_split.mat','u','S','Y','TR','TE','indices')
    
u_train = u(:,TR(1,:));
S_train = S(:,:,TR(1,:));
Y_train = Y(TR(1,:));
N_train = size(Y_train,1);

u_test = u(:,TE(1,:));
S_test = S(:,:,TE(1,:));
Y_test = Y(TE(1,:));
N_test = size(Y_test,1);

dim = size(u,1);
M_0 = eye(dim);
X_0 = eye(dim);
    
tic
[l,u] = ComputeBWDExtremes(u_train,S_train,5,95);
C = GetConstraints(Y_train,numConst,l,u);
[M_sol,X_sol,b_sol] = bwdml(M_0,X_0,u_train,S_train,C,gamma,maxItrs,tol);
toc
   
tic
pred = GBW_KNN(Y_train,u_train,S_train,M_sol,X_sol,knn_size,S_test);
acc = sum(pred==Y_test)/size(Y_test,1);
toc

