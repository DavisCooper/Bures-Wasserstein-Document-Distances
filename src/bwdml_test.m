clear all
addpath('../util')
addpath('../algs')
warning('off','MATLAB:sqrtm:SingularMatrix');

% define grid parameter search space
tol = 1e-6;
maxItrs = 1e5;
numConst = 1e2;
gamma = 1e-3;
knn_size = 4;

% load datasets
load('../data/ETH-80_RCD.mat','N','N_train','N_test','dim','D','D_train','D_test',...
    'y','y_train','y_test','datasetName','knnSize');

% select training and test data if necessary
if ~exist('D_train','var') 
    c = cvpartition(N,'Holdout',0.2);
    
    D_train = D(:,:,training(c));
    y_train = y(training(c));
    N_train = size(y_train,1);
    
    D_test = D(:,:,test(c));
    y_test = y(test(c));
    N_test = size(y_test,1);
end

if N_test > 1e2
    I = randperm(N_test,1e2);
    D_test = D_test(:,:,I);
    y_test = y_test(I,:);
    N_test = 1e2;
end

X_train = nan*ones(size(D_train));
for i=1:N_train
    X_train(:,:,i) = sqrtm(D_train(:,:,i));
end

X_test = nan*ones(size(D_test));
for i=1:N_test
    X_test(:,:,i) = sqrtm(D_test(:,:,i));
end

Q_0 = eye(dim);
    
tic
[l,u] = ComputeBWDExtremes(D_train,5,95);
C = GetConstraints(y_train,numConst,l,u);
[Q_sol,b_sol] = gbwml_seq_proj(Q_0,X_train,C,gamma,maxItrs,tol);
toc
   
% train_idx = unique([C(:,1);C(:,2)]);
% D_train = D_train(:,:,train_idx);
% X_train = X_train(:,:,train_idx);
tic
pred = GBW_KNN(y_train,D_train,Q_sol,knn_size,D_test);
acc = sum(pred==y_test)/size(y_test,1);
toc

