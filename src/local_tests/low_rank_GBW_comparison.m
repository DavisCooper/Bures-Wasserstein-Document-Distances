n = 300;
p = 250;

A = randn(n,n);
B = randn(n,n);
L = randn(n,p);
L = L/sqrtm(L'*L);

A = A*A';
B = B*B';
% M = M/sqrtm(M'*M);
M = L*L';

sqrt_A = sqrtm(A);
sqrt_B = sqrtm(B);

% compute polar factor 
tic;,
[U,S,V] = svds(sqrt_B*M*sqrt_A,p); 
t_svd_s = toc;
W_svd_p = U*V';

tic;
[U,S,V] = svd(sqrt_B*M*sqrt_A); 
t_svd = toc;
W_svd = U*V';
fprintf('SVDS Elapsed Time: %f, SVD Elapsed Time: %f\n', t_svd_s, t_svd)

W = sqrt_B*M*sqrt_A/sqrtm(sqrt_A*M*B*M*sqrt_A);

% problem.M = stiefelfactory(n,n,1);
% problem.cost = @(W) norm(L'*(sqrt_A - sqrt_B*W),'fro')^2;
% problem.egrad = @(W) 2*sqrt_B*M*(sqrt_B*W - sqrt_A);
% W0 = eye(n);
% options = struct;
% options.tolgradnorm = 1e-9;
% checkgradient(problem);
% [W_, Wcost, info, options] = trustregions(problem,W0,options);


f1 = trace((sqrt_A - sqrt_B*W_svd_p)'*M*(sqrt_A - sqrt_B*W_svd_p));
f2 = trace((sqrt_A - sqrt_B*W_svd)'*M*(sqrt_A - sqrt_B*W_svd));
f3 = trace((sqrt_A - sqrt_B*W)'*M*(sqrt_A - sqrt_B*W));
f = trace(M*A + M*B - 2*sqrtm(sqrt_A*M*B*M*sqrt_A));
f_ = trace(L'*A*L + L'*B*L - 2*sqrtm(L'*B*L*L'*A*L));

% slns = [f1, f2, f3, Wcost, f];