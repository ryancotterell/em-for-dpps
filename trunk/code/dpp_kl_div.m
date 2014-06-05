function res = dpp_kl_div(K1, K2, num_samples)
% Approximates the KL divergence between the DPP with marginal
% kernel K1 and the DPP with marginal kernel K2 by drawing
% num_samples samples from K1's distribution.  Assumes K1 and K2 are
% already decomposed (see decompose.m).

N = size(K1.M, 1);

% Converting to L in not an option, as the K could have eigenvalues of
% exatly 1, making L infinite.
diag_idxs = 1:N+1:N^2;
K1_diag = diag(K1.M);
K2_diag = diag(K2.M);

res = 0;
for j = 1:num_samples
  S = sample_dpp_via_K(K1);
  S_bar = ones(N, 1);
  S_bar(S) = 0;
  
  K1.M(diag_idxs) = K1_diag - S_bar;
  p1 = abs(det(K1.M));
  if p1 == 0
    continue;
  end

  K2.M(diag_idxs) = K2_diag - S_bar;
  p2 = abs(det(K2.M));
    
  if (p2 == 0) && (p1 > 0)
    res = Inf;
    break;
  end
  
  res = res + log(p1 / p2) / num_samples;
end
