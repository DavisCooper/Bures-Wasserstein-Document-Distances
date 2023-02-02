function [y,u,S] = synthetic_data(d,m,n_c,alpha,beta)

    y = zeros(m,1);
    u = zeros(d,m*n_c);
    S = zeros(d,d,m*n_c);
    for i=1:n_c
        % set mean and covariance
        u_c = alpha*randn(d,1);
        V_c = beta*randn(d,d);
        V_c = V_c*V_c';
        for j=1:m
            k = j + (i-1)*m;
            y(k) = i;
            u(:,k) = u_c + randn(d,1);
            S(:,:,k) = randn(d,d);
            S(:,:,k) = V_c + S(:,:,k)*S(:,:,k)';
        end
    end

end