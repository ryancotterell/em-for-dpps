function S = sample_dpp_via_K(K)
% Draws a single sample from a DPP with the given marginal kernel.
%   K = N x N marginal kernel matrix, decomposed (see decompose.m)
  
% Choose eigenvectors randomly.
v = find(rand(length(K.D), 1) <= K.D);  

% Sample points according to their support under v.
S = sample_points_via_K(K.V, v);
