function S = sample_k(lambda, remaining, E)
% Used internally by sampl_dpp to select remaining eigenvalues from
% i downwards using elem symm polys E.
  
% Iterate from i down to 1.
i = numel(lambda);
S = zeros(remaining, 1);
while remaining > 0

  % Compute marginal of i given that we choose remaining values from 1:i.
  if i == remaining
    marg = 1;
  else
    if lambda(i) == 0
        marg = 0;
    else
      marg = lambda(i) * E(remaining, i - 1) / E(remaining + 1, i);
      %assert(marg <= 1, ['Marginal too large: ' num2str(marg)]);
      %assert(marg >= 0, ['Marginal too small: ' num2str(marg)]);
    end
  end
    
  % Sample marginal.
  if rand < marg
    S(remaining) = i;
    remaining = remaining - 1;        
  end
  i = i-1;
end
