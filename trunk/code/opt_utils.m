function utils = opt_utils()
% Contains dataset generators, K initializers, likelihood
% calculators, and generic parameter optimization functions.  See
% synth_test.m for example usage.

utils.build_training_set = @build_training_set;
utils.sort_training_set = @sort_training_set;

utils.K_wishart_init = @K_wishart_init;
utils.K_clusters_init = @K_clusters_init;
utils.get_low_moments = @get_low_moments;
utils.K_moments_init = @K_moments_init;
utils.K_perturbation_init = @K_perturbation_init;
utils.K_Dshaken_init = @K_Dshaken_init;
utils.K_Vrotated_init = @K_Vrotated_init;

utils.K_log_likelihood = @K_log_likelihood;
utils.L_log_likelihood = @L_log_likelihood;

utils.take_one_step = @take_one_step;
utils.optimize_param = @optimize_param;
utils.get_max_change = @get_max_change;

utils.get_em_opts = @get_em_opts;
utils.get_K_ascent_opts = @get_K_ascent_opts;
utils.get_K_expgrad_opts = @get_K_expgrad_opts;


%%%%%%%%%%%%%%%%%%%%%%%%%%% Dataset Generators %%%%%%%%%%%%%%%%%%%%%%%%%%%


function [T, Ys_map] = build_training_set(L, n, min_Y_size)
% Record sizes for convenience.
T.N = numel(L.D);
T.n = n;

% Use a map from sets to counts to build up a de-duplicated training set.
import EqualByValueIntArray;
Ys_map = java.util.HashMap;
for i = 1:n
  % Sample a set.
  while 1
    Y = sample_dpp(L);
    if numel(Y) >= min_Y_size
      break
    end
  end
  Y = EqualByValueIntArray(Y);
  
  % Check to see if it has been sampled before.
  Y_count = Ys_map.get(Y);
  if isempty(Y_count)
    % If not, add the new set to the map.
    Ys_map.put(Y, 1);
  else
    % Else, increment the count of this set.
    Ys_map.put(Y, Y_count + 1);
  end
end
T.n_dedup = Ys_map.size();

% Dump the map into an array; we won't need random access anymore.
map_iter = Ys_map.entrySet().iterator();
T.Ys = cell(T.n_dedup, 1);
T.Y_counts = zeros(T.n_dedup, 1);
T.Y_fracs = zeros(T.n_dedup, 1);
i = 0;
while map_iter.hasNext()
  map_entry = map_iter.next();
  
  i = i + 1;
  T.Ys{i} = map_entry.getKey().getArray();
  T.Y_counts(i) = map_entry.getValue();
  T.Y_fracs(i) = T.Y_counts(i) / T.n;
end

% Record a few other useful shortcuts for the data.
T.Y_sizes = zeros(T.n_dedup, 1);
T.Y_bar_inds = cell(T.n_dedup, 1);
T.Y_bars = cell(T.n_dedup, 1);
for i = 1:T.n_dedup
  Y = T.Ys{i};
  
  % Save the set sizes.
  T.Y_sizes(i) = numel(Y);
  
  % Save the complement set.
  bar_ind = true(T.N, 1);
  bar_ind(Y) = false;
  T.Y_bar_inds{i} = bar_ind;
  T.Y_bars{i} = find(bar_ind);
end


function T = sort_training_set(T)
% Sort the sets by size.
[T.Y_sizes, order] = sort(T.Y_sizes);
T.Ys = T.Ys(order);
T.Y_counts = T.Y_counts(order);
T.Y_fracs = T.Y_fracs(order);
T.Y_bar_inds = T.Y_bar_inds(order);
T.Y_bars = T.Y_bars(order);

% Mark the start of a new size by placing its index in Y_size_starts.
T.Y_size_starts = [-1; T.Y_sizes - circshift(T.Y_sizes, [-1, -1])];
T.Y_size_starts = find(T.Y_size_starts(1:end-1) ~= 0);

% Record the max, min, and average set size.
T.avg_Y_size = sum(T.Y_sizes .* T.Y_fracs);
T.max_Y_size = T.Y_sizes(end);
T.min_Y_size = T.Y_sizes(1);

% Store a vector of unique set sizes.
T.unique_Y_sizes = T.Y_sizes(T.Y_size_starts);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% K Initializers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function K = K_wishart_init(N, L_scaling, max_L_eig)
% Draw a random L from Wishart(N, I).
B = randn(N);
L = decompose((L_scaling / N) * B * B');

% Enforce max_L_eig constraint.
L.D = min(L.D, max_L_eig);
L.M = L.V * diag(L.D) * L.V';

% Transform from L to K.
K = eye(N) - inv(L.M + eye(N));
K = (K + K') / 2;
K = decompose(K);


function K = K_clusters_init(N, cluster_closeness, kernel_closeness)
% Model N/10 clusters in a 2D space by drawing from a Gaussian
% mixture model with N/10 components.
num_clusters = ceil(N / 10);
Nf = 2;
mus = rand(num_clusters, 2);
sigmas = cluster_closeness * ones(1, Nf);
props = ones(1, num_clusters) / num_clusters;
gmd = gmdistribution(mus, sigmas, props);
[B, cluster_idxs] = random(gmd, N);

% Build an RBF kernel based on points drawn from the mixture model.
gamma = 1 / (0.5 * (kernel_closeness)^2);
P = B * B';
dP = diag(P);
L = exp(-gamma * bsxfun(@plus, bsxfun(@plus, dP, dP'), -2 * P));

% Convert from L to K.
L = decompose(L);
K = eye(N) - inv(L.M + eye(N));
K = (K + K') / 2;
K = decompose(K);


function [Q, Ys_ind] = get_low_moments(T)
% Represent each Y as a 0-1 bit string, multiplied by its count.
Ys_ind = zeros(T.n_dedup, T.N);
for i = 1:T.n_dedup
  Ys_ind(i, T.Ys{i}) = sqrt(T.Y_counts(i));
end

% Set diagonal (i, i) to the fraction of sets in which element i occurs,
% and set off-diagonal (i, j) to the fraction of sets in which both
% elements i and j occur.
Q = zeros(T.N);
for i = 1:T.N
  Q(i, :) = sum(bsxfun(@times, Ys_ind(:, i), Ys_ind), 1);
end
Q = Q / T.n;


function [K, Ys_ind] = K_moments_init(T, max_eig)
[K, Ys_ind] = get_low_moments(T);

% Transform off-diagonal such that: det(K_{ij}) = K_{ii}K_{jj} - K_{ij}^2.
pairs = nchoosek(1:T.N, 2);
num_pairs = size(pairs, 1);
for i = 1:num_pairs
  pair = pairs(i, :);
  a = pair(1);
  b = pair(2);
  K(a, b) = K(a, a) * K(b, b) - K(a, b);
  if K(a, b) < 0
    % The best we can do while keeping singleton probabilities intact is to 
    % set these to zero.
    K(a, b) = 0;
    K(b, a) = 0;
  else
    K(a, b) = sqrt(K(a, b));
    K(b, a) = K(a, b);
  end
end

% Project K to the PSD cone and clip it's eigenvalues at 1.
[V, D] = eig(K);
D = min(max(diag(D), 0), max_eig);
zdiff = T.max_Y_size - sum(D > 0);
if zdiff > 0
  Dz = find(D == 0);
  min_eig = 1 - max_eig;
  D(Dz(1:zdiff)) = min_eig;
end
K = [];
K.M = V * diag(D) * V';
K.M = (K.M + K.M') / 2;
K.V = V;
K.D = D;


function K = K_Dshaken_init(K, shake_strength)
% Change eigenvalues by at most shake_strength, truncating to keep a small
% margin between their values and 0, 1.
lower_limit = 0.01;
upper_limit = 0.99;
N = size(K.M, 1);
K.D = max(min(K.D + 2 * shake_strength * (rand(N, 1) - 0.5), ...
  upper_limit), lower_limit);
K.M = K.V * diag(K.D) * K.V';


function K = K_Vrotated_init(K, rotation_strength)
% Change the eigenvectors by rotating them by a random rotation matrix.
% Rotation strength should be a number in [0, 1].
N = size(K.M, 1);
[Q, R] = qr(randn(N));
Q = Q * diag(sign(diag(R)));
if det(Q) < 0
  Q(:, 1) = -Q(:, 1);
end
K.V = real(Q^(rotation_strength) * K.V);
K.M = K.V * diag(K.D) * K.V';


function K = K_perturbation_init(K, shake_strength, rotation_strength)
% Perturb eigenvalues and eigenvectors.
K = K_Dshaken_init(K, shake_strength);
K = K_Vrotated_init(K, rotation_strength);

% Ensure perfect symmetry of the final full matrix.
K.M = (K.M + K.M') / 2;


%%%%%%%%%%%%%%%%%%%%%%%%%%% Objective Computers %%%%%%%%%%%%%%%%%%%%%%%%%%%


function obj = K_log_likelihood(T, K)
obj = 0;
K_diag = diag(K);
diag_idxs = 1:T.N+1:T.N^2;
for i = 1:T.n_dedup
  K(diag_idxs) = K_diag - T.Y_bar_inds{i};
  obj = obj + T.Y_fracs(i) * log(abs(det(K)));
end


function obj = L_log_likelihood(T, L)
obj = 0;
for i = 1:T.n_dedup
  Y = T.Ys{i};
  obj = obj + T.Y_fracs(i) * log(det(L(Y, Y)));
end
obj = obj - log(det(L + eye(T.N)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Gradient Ascender %%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [A, obj, step_size] = take_one_step(param_name, A_old, Ag, ...
  step_size, step_func, obj_old, obj_func, min_step_size)
% Note: Wolfe rule line search would be nice to use, but it only applies
%       when the update is of the form A_old + step_size * Ag.

while 1
  % Take a step in the direction of the gradient.
  A = step_func(step_size, A_old, Ag);
  
  % If the objective value has not increased, decrease the step size.
  obj = obj_func(A);
  fprintf('Objective at step size %f for %s update: %f\n', ...
    step_size, param_name, obj);
  
  if obj > obj_old
    break;
  else
    step_size = step_size / 2;
    if step_size < min_step_size
      % This can happen if either: 1) the estimate of the objective is
      % poor, or 2) the gradient is directly perpendicular to the
      % constraint set and A is already right at the edge of it.
      warning('Minimum step size (%f) reached; %s left as is.', ...
        min_step_size, param_name);
      A = A_old;
      obj = obj_old;
      break;
    end
  end
end


function [A, obj_vals, iter] = optimize_param(param_name, A, ...
  grad_func, step_size, step_func, obj, obj_func, min_step_size, ...
  short_get_max_change, min_obj_change, max_iters)
iter = 0;
if max_iters == 0
  return;
end

obj_vals = obj;
while 1
  fprintf('%s update iteration %d\n', param_name, iter + 1);
  A_old = A;
  obj_old = obj;
  Ag = grad_func(A);

  [A, obj, step_size] = take_one_step(param_name, A_old, Ag, step_size, ...
    step_func, obj_old, obj_func, min_step_size);
  iter = iter + 1;
  obj_vals(end + 1) = obj;

  % Check optimization-ending conditions.  
  obj_change = short_get_max_change(obj_old, obj);
  if obj_change < min_obj_change
    fprintf(['Objective change (%f) negligible; exiting %s ' ...
      'optimization.\n'], obj_change, param_name);
    break;
  end
  if iter >= max_iters
    fprintf(['Max # of iterations (%d) reached; exiting %s ' ...
      'optimization.\n'], max_iters, param_name);
    break;
  end
end


function max_change = get_max_change(old_vals, new_vals, zero_limit)
max_change = max(abs(old_vals(:) - new_vals(:)) ./ ...
  max(abs(old_vals(:)), zero_limit));


%%%%%%%%%%%%%%%%%%%%%%%% Parameter Option Setters %%%%%%%%%%%%%%%%%%%%%%%%


% Sets up all the optimization parameters for EM.
function opts = get_em_opts(N)
% Stopping criteria: stop if the max change in the objective's value is
% less than (100 * min_obj_change)% of the objective's previous value.
opts.min_obj_change = 1e-4;
% In the case that the objective's previous value is (near) zero, the
% difference between it and the current value will be considered relative
% to zero_limit instead.
opts.zero_limit = 1e-5;
% Absolute maximum number of iterations.
opts.max_em_iters = 100 * N^2;
% Below this step size, the V-updating code will stop trying to change V.
opts.min_step_size = 1e-8;
% Eigenvalue cap (so that q's L-kernels won't be infinite).
opts.max_eig = 0.99;


% Sets up all the optimization parameters for projected gradient ascent
% on K.
function opts = get_K_ascent_opts(N)
% Same as for EM.
opts.min_obj_change = 1e-4;
opts.zero_limit = 1e-5;
opts.max_iters = 100 * N^2;
% Below this step size, the K-updating code will stop trying to change K.
opts.min_step_size = 1e-8;
% Eigenvalue lower bound (only used with exponentiated gradient).
opts.min_eig = 0.01;

