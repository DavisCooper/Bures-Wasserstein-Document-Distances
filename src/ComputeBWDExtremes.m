function [l, u] = ComputeBWDExtremes(u, S, a, b)
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
    [U,~,V] = svd(S(:,:,j2)*S(:,:,j1)); % SVD calculation is faster than matrix square root
    W = U*V';
    dists(i) = norm(u(:,j1) - u(:,j2))^2 + norm(S(:,:,j1) - S(:,:,j2)*W)^2;
end

[~, c] = hist(dists, 100);
l = c(floor(a));
u = c(floor(b));