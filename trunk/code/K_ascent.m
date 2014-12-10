function [K, obj_vals] = K_ascent(T, opts, utils, K_start, type)
% Gradient ascent on a DPP's marginal kernel, K.
%   T = training set (see opt_utils.m's build_training_set)
%   opts = convergence options (see opt_utils.m's get_K_ascent_opts)
%   utils = access to generic optimization code (typically opt_util.m's utils) 
%   K_start = initial K from which to begin gradient ascent

if strcmp(type, 'pg')
  K_step_func = @K_pg_step_func;
else
  if strcmp(type, 'eg')
    % Enforce min eig constraint if running exponentitated gradient.
    K_start.D = max(K_start.D, opts.min_eig);
    K_start.M = K_start.V * diag(K_start.D) * K_start.V';
    K_start.M = (K_start.M + K_start.M') / 2;
    
    K_step_func = @(step_size, K, Kg) ...
      K_eg_step_func(step_size, K, Kg, opts.min_eig);
  else
    throw(MException('K_ascent:BadType', ...
      ['Recognized types are projected gradient (pg) and' ...
       'exponentiated gradient (eg).']));
  end
end

% Run the optimization.
K_start = K_start.M;
initial_step_size = 1;
[K, obj_vals] = utils.optimize_param('K', K_start, ...
  @(K) K_grad_func(T, K), ...
  initial_step_size, ...
  K_step_func, ...
  utils.K_log_likelihood(T, K_start), ...
  @(K) utils.K_log_likelihood(T, K), ...
  opts.min_step_size, ...
  @(old_vals, new_vals) utils.get_max_change(old_vals, new_vals, ...
    opts.zero_limit), ...
  opts.min_obj_change, opts.max_iters);

% Make sure symmetry is perfect, so that call to decompose
% doesn't result in imaginary eigenvectors.
K = (K + K') ./ 2;
K = decompose(K);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Helpers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
function Kg = K_grad_func(T, K)
% Compute the gradient of the likelihood with respect to K.
Kg = zeros(T.N);
K_diag = diag(K);
diag_idxs = 1:T.N+1:T.N^2;
for i = 1:T.n_dedup
  K(diag_idxs) = K_diag - T.Y_bar_inds{i};
  if cond(K) > 1 / eps
    % Matrix will not be invertible.
    continue;
  end
  
  Kg = Kg + T.Y_fracs(i) * inv(K);
end


function K = K_pg_step_func(step_size, K, Kg)
% Take a step in the direction of the gradient.
K = K + step_size * Kg;

% Project back to the constraint space.
[V, D] = eig(K);
D = diag(D);
D = min(max(real(D), 0), 1);
K = V * diag(D) * V';


function K = K_eg_step_func(step_size, K, Kg, min_eig)
% Take an exponentiated gradient step, truncating eigenvalues at 1 and (so
% that logm(K) is feasible) min_eig.
[V, D] = eig(logm(K) + step_size * Kg);
D = diag(D);
D = max(min(exp(real(D)), 1), min_eig);
K = V * diag(D) * V';
