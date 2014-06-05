function C = k_greedy_map(L, k)
% Attempts to find the highest-probability size-k set under L's DPP.
%   L = a (non-marginal) DPP kernel (no need for decomposition)
%   k = desired set size

% Greedily add elements one at a time.
C = [];
N = size(L, 1);
U = 1:N;
num_left = N;
while numel(C) < k
  % Compute the determinant with each remaining unused element added to
  % the chosen set.
  scores = diag(L);
  
  % Select the max-scoring addition to the chosen set.
  [max_score, max_loc] = max(scores);
  C = [C; U(max_loc)];
      
  % Switch row max_loc with row 1; switch col max_loc with col 1.
  elim_row = L(max_loc, :);
  if max_loc ~= 1
    elim_row(1, max_loc) = elim_row(1, 1);
    L(max_loc, :) = L(1, :);
    L(:, max_loc) = L(:, 1);
    U(max_loc) = U(1);
  end
    
  % Remove first row and col.
  elim_row(1) = [];
  U(1) = [];
  L = L(2:end, 2:end);
    
  % For each remaining item, clear its first col.
  for i = 1:num_left-1
    L(i, :) = L(i, :) - (elim_row(i) / max_score) * elim_row;
  end
 
   num_left = num_left - 1;
end
