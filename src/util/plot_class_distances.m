function [] = plot_class_distances(y,u,S,M,X)

n_c = size(unique(y),1);
m = size(u,2);

for i=1:n_c
    for j=1:n_c

        h = [];
        idx = 1;
        for k=1:m
            for l=1:m
                if y(k) == i & y(l) == j
                    [U,~,V] = svd(S(:,:,k)*X*S(:,:,l)); % SVD calculation is faster than matrix square root
                    W = U*V';
                    V = S(:,:,k) - S(:,:,l)*W;
                    h(idx) = (u(:,k) - u(:,l))'*M*(u(:,k) - u(:,l)) + trace(V'*X*V);
                    idx = idx + 1;
                end

            end
        end
        subplot(n_c,n_c, j + (i-1)*n_c)
        hold on
        histogram(h)
    end
end

end