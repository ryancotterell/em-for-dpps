function S = sample_dpp(L, k)
% Draws a single sample from a DPP with the given (non-marginal) kernel.
%   L = N x N kernel matrix, decomposed (see decompose.m)
%   k = desired sample size (optional)
  
if ~exist('k', 'var')
  % Choose eigenvectors randomly.
  D = L.D ./ (1 + L.D);
  v = find(rand(length(D), 1) <= D);
else
  % Select k eigenvectors.
  E = elem_symm_polys(L.D, k);
  v = sample_k(L.D, k, E);
end

% Sample points according to their support under v.
S = sample_points_via_K(L.V, v);
