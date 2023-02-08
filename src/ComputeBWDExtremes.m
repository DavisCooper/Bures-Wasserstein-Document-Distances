function [l, u] = ComputeBWDExtremes(u, S, a, b, M_0, X_0)
% [l, u] = ComputeDistanceExtremes(X, a, b, M)
%
% Computes sample histogram of the distances between rows of X and returns
% the value of these distances at the a^th and b^th percentils.  This
% method is used to determine the upper and lower bounds for
% similarity / dissimilarity constraints.  
%
% X: (n x m) data matrix 
% a: lower bound percentile between 1 and 100
% b: upper bound percentile between 1 and 100
% M: Mahalanobis matrix to compute distances 
%
% Returns l: distance corresponding to a^th percentile
% u: distance corresponding the b^th percentile

if (a < 1 || a > 100),
    error('a must be between 1 and 100')
end
if (b < 1 || b > 100),
    error('b must be between 1 and 100')
end

n = size(S, 3);

num_trials = min(100, n*(n-1)/2);

% we will sample with replacement
dists = zeros(num_trials, 1);
for i=1:num_trials
    j1 = ceil(rand(1)*n);
    j2 = ceil(rand(1)*n);    

    [U,~,V] = svd(S(:,:,j2)*X_0*S(:,:,j1)); % SVD calculation is faster than matrix square root
    W = U*V';
    
    u_ij = u(:,j1) - u(:,j2);
    S_ij = S(:,:,j1) - S(:,:,j2)*W; 
    
    dists(i) = u_ij'*M_0*u_ij + trace(S_ij*S_ij'*X_0);
end

[~, c] = hist(dists, 100);
l = c(floor(a));
u = c(floor(b));