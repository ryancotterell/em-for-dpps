function S = sample_points_via_K(V_full, v)
% Used internally by sample_dpp as step-two of the sampling process.
%   V_full = all eigenvectors (one per column)
%   v = indices of selected eigenvectors

k = numel(v);
N = size(V_full, 1);
V = V_full(:, v);

% Compute initial marginals.
K_diag = sum(V .* V, 2);

% Iterate, selecting one new item each iteration.
prev_Kcs = zeros(size(V));
S = zeros(k, 1);
for i = 1:k
  % Compute probability for each item.
  P = K_diag / sum(K_diag);

  % Choose a new item to include.
  S(i) = find(rand <= cumsum(P), 1);

  % Get row of chosen point from current conditional K.
  Kcs = V * V(S(i), :)';
  if i > 1
    J = 1:i-1;
    Kcs = Kcs - sum(bsxfun(@times, prev_Kcs(S(i), J) ./ ...
      prev_Kcs(S(J)' + (J - 1) * N), prev_Kcs(:, J)), 2);
  end
  prev_Kcs(:, i) = Kcs;

  % Update K's diagonal.
  K_diag = K_diag - Kcs.^2 ./ K_diag(S(i));
  %assert(all(K_diag(S(1:i)) < 1e-10));
  %assert(all(K_diag(S(1:i)) > -1e-10));
  
  % Prevent numerical errors from propagating.
  K_diag(S(1:i)) = 0;
  
  % Note: For faster sampling, could decrease size of K over iterations:
  %         K_diag([1:S(i)-1, S(i)+1:end]);
  %       Just need bookkeeping for the selections.
end

S = sort(S);
