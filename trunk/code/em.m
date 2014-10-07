function [K, obj_vals] = em(T, opts, utils, K_start)
% Expectation maximization for a DPP.
%   T = training set (see opt_utils.m's build_training_set)
%   opts = convergence options (see opt_utils.m's get_K_ascent_opts)
%   utils = access to generic optimization code (typically opt_util.m's utils) 
%   K_start = initial K from which to begin gradient ascent

if opts.max_em_iters == 0
  K = K_start;
  return;
end

% Initialize q, p.
V = K_start.V;
D = K_start.D;

% Iterate E and M steps.
R = sqrt(D ./ (1 - D));
prev_log_like = -Inf;
em_iter = 0;
step_size = 1;
obj_vals = [];
while 1
  em_iter = em_iter + 1;
  fprintf('EM iteration %d\n', em_iter);
  
  % Iterate through the examples once.
  log_like = 0;
  VRV = bsxfun(@times, V, R'.^2) * V';
  VRV = (VRV + VRV') / 2;
  D_new = zeros(T.N, 1);
  Vg = zeros(T.N);
  for i = 1:T.n_dedup
    k = T.Y_sizes(i);
    if k == 0
      continue;
    end
    
    % Get the eigendecomposition associated with q(J | Y_i).
    Y = T.Ys{i};
    qYC = decompose(VRV(Y, Y));
    qYLCh = bsxfun(@times, R, V(Y, :)');
    qYLV = qYLCh * bsxfun(@times, qYC.V, 1 ./ sqrt(qYC.D'));
    
    % Eigenvalue update.
    D_new = D_new + T.Y_fracs(i) * sum(qYLV.^2, 2);
    
    % V gradient accumulation.
    Vg(Y, :) = Vg(Y, :) + 2 * T.Y_fracs(i) * inv(VRV(Y, Y)) * ...
      bsxfun(@times, R, qYLCh)';
    
    % Objective value accumulation.
    log_like = log_like + T.Y_fracs(i) * sum(log(qYC.D));
  end
  
  % Compute the overall objective at the beginning of the M-step; when
  % q = p_K, this is just the log likelihood.
  log_like = log_like + sum(log(1 - D));
  
  % Make sure the log likelihood has increased since the last iteration.
  fprintf('Log likelihood: %f\n', log_like);
  obj_vals(end + 1) = log_like;
  if em_iter > 1
    obj_diff = log_like - prev_log_like;
    assert(obj_diff >= -1e-5 * abs(prev_log_like));
    obj_change = obj_diff / max(abs(prev_log_like), opts.zero_limit);
    if obj_change < opts.min_obj_change
      fprintf('Exiting EM due to negligible objective change.\n');
      break;
    end
  end
  prev_log_like = log_like;
  
  % Sightly decrease any eigenvalues that exceed max_eig.
  D_new(D_new > opts.max_eig) = opts.max_eig;
  
  % Transform the standard V gradient so that it lies in the tangent space
  % of the Stiefel manifold.
  Vg = V' * Vg - Vg' * V;
  
  % Take just one gradient step for V.
  R_new = sqrt(D_new ./ (1 - D_new));
  [V_new, ~, step_size] = utils.take_one_step('V', V, Vg, ...
    step_size, @V_step_func, V_obj_func(T, V, D_new, R_new), ...
    @(V) V_obj_func(T, V, D_new, R_new), ...
    opts.min_step_size);
  
  % E-step: Set q = p_K.
  D = D_new;
  R = R_new;
  V = V_new;
  
  if em_iter == opts.max_em_iters
    fprintf('Exiting EM after executing maximum # of iterations.\n');
    break;
  end
end

% Construct K from the learned D and V.
K.D = D;
K.V = V;
K.M = K.V * diag(K.D) * K.V';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Helpers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function V = V_step_func(step_size, V, Vg)
V = V * expm(step_size * Vg);


function log_like = V_obj_func(T, V, D, R)
log_like = 0;
VRV = bsxfun(@times, V, R'.^2) * V';
for i = 1:T.n_dedup
  if T.Y_sizes(i) == 0
    continue;
  end
  
  Y = T.Ys{i};
  log_like = log_like + T.Y_fracs(i) * log(det(VRV(Y, Y)));
end
log_like = log_like + sum(log(1 - D));
