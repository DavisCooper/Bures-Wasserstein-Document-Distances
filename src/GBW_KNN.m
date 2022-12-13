function preds = GBW_KNN(Y, u, S, M, X, k, St)
% function preds = KNN(y, X, M, k, Xt)
% perform knn classification on each row of Xt
% M is the factor matrix : A = M*M'


add1 = 0;
if (min(Y) == 0),
    Y = Y + 1;
    add1 = 1;
end

n = size(S,3);
nt = size(St,3);
D = inf*ones(n,nt);
for i=1:n
    for j=1:nt
%         D(i,j) = trace(M*X(:,:,i) + M*Xt(:,:,j) - 2*sqrtm(M*X(:,:,i)*M*Xt(:,:,j)));
        [U,~,V] = svd(St(:,:,j)*M*S(:,:,i)); % SVD calculation is faster than matrix square root
        W = U*V';
        V = S(:,:,i) - St(:,:,j)*W;
        D(i,j) = (u(:,i) - u(:,j))'*M*(u(:,i) - u(:,j)) + trace(V'*X*V);
    end
end

[V, Inds] = sort(D);

preds = zeros(nt,1);
for (i=1:nt),
    counts = [];
    for (j=1:k),        
        if (Y(Inds(j,i)) > length(counts)),
            counts(Y(Inds(j,i))) = 1;
        else
            counts(Y(Inds(j,i))) = counts(Y(Inds(j,i))) + 1;
        end
    end
    [v, preds(i)] = max(counts);
end
if (add1 == 1),
    preds = preds - 1;
end

end