function E = elem_symm_polys(lambda, k)
% Builds a matrix of elementary symmetric polynomials:
%   E(i + 1, j) =
%     \sum_{|R| = i : m \notin R if m > j} \prod_{r \in R} lambda_r.
% The k-th order polynomial for N lambdas is entry E(k + 1, N).

N = numel(lambda);
assert(N == length(lambda));
assert(N >= k);
assert(k >= 0);
E = zeros(k + 1, N);

E(1, :) = 1;
if N > 0
  E(2, 1) = lambda(1);
  for i = 2:N
    j = 2:min(i + 1, k + 1);
    E(j, i) = E(j, i - 1) + lambda(i) * E(j - 1, i - 1);
  end
end