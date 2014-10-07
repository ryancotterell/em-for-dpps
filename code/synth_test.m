function synth_test(Ns, start_type, result_file)
% Tests learning methods on synthetic data.  See README.txt for an
% example usage.  Lines marked with "% *" here indicate internal
% parameters.

save(result_file, 'Ns', '-v7.3');

% Iterate over the requested kernel sizes.
results = [];
utils = opt_utils();
for N_num = 1:numel(Ns)
  N = Ns(N_num);

  % Check a variety of training set sizes.
  T_sizes = 5 * N;  % *
  while T_sizes(end) < 10 * N  % *
    T_sizes(end+1) = T_sizes(end) * 2;  % *
  end
  T_sizes(end) = 10 * N;  % *
  results(N_num).T_sizes = T_sizes;

  % K_start = moment matching.
  if strcmp(start_type, 'mm')
    max_eig = 0.99;  % *
    results(N_num).mm_start = ...
      test_all_target_types(N, T_sizes, ...
        @(T, K) utils.K_moments_init(T, max_eig), 0);
  end
  
  % K_start = moment matching where training data does not contain
  % small set sizes.
  if strcmp(start_type, 'mc')
    max_eig = 0.99;  % *
    min_Y_sizes = [2, 5];
    results(N_num).min_Y_sizes = min_Y_sizes;
    for i = 1:numel(min_Y_sizes)
      results(N_num).mc_start(i) = ...
      test_all_target_types(N, T_sizes, ...
			    @(T, K) utils.K_moments_init(T, max_eig), min_Y_sizes(i));
      save(result_file, 'results', '-append', '-v7.3');
    end
  end
  
  % K_start = Wishart(N, I).
  if strcmp(start_type, 'wi')
    max_L_eig = 0.99 / (1 - 0.99);  % *
    L_scalings = [0.5, 1, 2];  % *
    results(N_num).L_scalings = L_scalings;
    for i = 1:numel(L_scalings)
      results(N_num).wi_start(i) = ...
        test_all_target_types(N, T_sizes, ...
          @(T, K) utils.K_wishart_init(N, L_scalings(i), max_L_eig), 0);
      save(result_file, 'results', '-append', '-v7.3');
    end
  end

  % K_start = K_true with perturbed eigenvalues.
  if strcmp(start_type, 'pd')
    shake_strengths = [0.25, 0.5, 1];  % *
    results(N_num).shake_strengths = shake_strengths;
    for i = 1:numel(shake_strengths)
      results(N_num).Dshaken_start(i) = ...
        test_all_target_types(N, T_sizes, ...
          @(T, K) ...
            utils.K_Dshaken_init(K, shake_strengths(i)), 0);
      save(result_file, 'results', '-append', '-v7.3');
    end
  end
  
  % K_start = K_true with perturbed eigenvectors.
  if strcmp(start_type, 'pv')
    rotation_strengths = [0.25, 0.5, 1];  % *
    results(N_num).rotation_strengths = rotation_strengths;
    for i = 1:numel(rotation_strengths)
      results(N_num).Vrotated_start(i) = ...
        test_all_target_types(N, T_sizes, ...
          @(T, K) ...
            utils.K_Vrotated_init(K, rotation_strengths(i)), 0);
      save(result_file, 'results', '-append', '-v7.3');
    end
  end

  save(result_file, 'results', '-append', '-v7.3');
end

save(result_file, 'results', '-append', '-v7.3');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Helpers %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function results = test_all_target_types(N, T_sizes, initializer, ...
  min_Y_size)
fprintf('Beginning test_all_target_types ...\n');
utils = opt_utils();

% Tests with target L ~ Wishart(N, I).
max_L_eig = 0.99 / (1 - 0.99);  % *
L_scalings = [0.5, 1, 2];  % *
wishart_tt.L_scalings = L_scalings;
for i = 1:numel(L_scalings)
  wishart_tt.results(i) = ...
    test_target_type(N, T_sizes, ...
		     @() utils.K_wishart_init(N, L_scalings(i), max_L_eig), ...
		     initializer, min_Y_size);
end

% Tests with target L an RBF kernel built from clustered points.
cluster_closenesses = [0.0005, 0.005];  % *
kernel_closenesses = [0.25, 0.5];  % *
cluster_tt.cluster_closenesses = cluster_closenesses;
cluster_tt.kernel_closenesses = kernel_closenesses;
for i = 1:numel(cluster_closenesses)
  for j = 1:numel(kernel_closenesses)
    cluster_tt.results(i, j) = ...
      test_target_type(N, T_sizes, ...
		       @() utils.K_clusters_init(N, cluster_closenesses(i), ...
						 kernel_closenesses(j)), ...
		       initializer, min_Y_size);
  end
end

results.wishart_tt = wishart_tt;
results.cluster_tt = cluster_tt;


function results = test_target_type(N, T_sizes, target_generator, ...
  initializer, min_Y_size)
fprintf('Beginning test_target_type ...\n');
utils = opt_utils();
eutils = eval_utils();

% Repeat this setting num_trials times.
num_trials = 10;  % *

% Get learning method options.
em_opts = utils.get_em_opts(N);
ka_opts = utils.get_K_ascent_opts(N);

% Run the trials and collect a few basic statistics.
num_stats = 39;
results = [];
num_T_sizes = numel(T_sizes);
results.stats = zeros(num_T_sizes, num_stats, num_trials);
results.quantiles = zeros(num_T_sizes, num_stats, 3);
for t = 1:num_trials
  % Build target K and L.  
  K = target_generator();
  L = K.M * inv(eye(N) - K.M);
  L = (L + L') / 2;
  L = decompose(L);
  results.targets(t).K = K;

  T_num = 0;
  for T_size = T_sizes
    T_num = T_num + 1;
    
    % Build train and test sets.
    T_train = utils.build_training_set(L, T_size, min_Y_size);
    T_train = utils.sort_training_set(T_train);
    T_test = utils.build_training_set(L, 10 * N, min_Y_size);  % *
    T_test = utils.sort_training_set(T_test);
    
    % Obtain an optimization starting point.
    K_start = initializer(T_train, K);
    results.Ks(T_num, t).start = K_start;

    % Run learning.
    tic;
    K_em = em(T_train, em_opts, utils, K_start);
    em_time = toc;
    results.Ks(T_num, t).em = K_em;
    tic;
    K_ka = K_ascent(T_train, ka_opts, utils, K_start, 'pg');  %'eg'
    ka_time = toc;
    results.Ks(T_num, t).ka = K_ka;
     
    % Compare methods.
    [kls, kl_diffs] = eutils.compare_learned_Ks(K, K_start, K_em, K_ka);
    [lls, ll_diffs] = ...
      eutils.compare_learned_Ks_likelihoods(K, K_start, K_em, K_ka, ...
      T_train, T_test);
    results.stats(T_num, :, t) = [kls, lls, em_time, ka_time, ...
      ll_diffs, kl_diffs, (em_time - ka_time) / em_time];
  end
end

% Compile summary statistics.
T_num = 0;
for T_size = T_sizes
  T_num = T_num + 1;
  results.quantiles(T_num, :, :) = ...
    prctile(reshape(results.stats(T_num, :, :), num_stats, num_trials), ...
      [25, 50, 75], 2);
end
