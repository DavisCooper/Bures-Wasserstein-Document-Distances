load('../data/bbcsport_tr_te_split.mat')
n = size(X,2);
d = size(X{1},1);
u = zeros(d,n);
S = zeros(d,d,n);
for i=1:n

    u(:,i) = sum(X{i}*BOW_X{i},2)/sum(BOW_X{i});
    S(:,:,i) = (X{i} - u)*(X{i} - u)';
end