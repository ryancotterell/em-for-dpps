function registries_test(sets_file, names_file, result_file, ...
			 train_frac, num_trials, start_type, max_numN_train)
% Tests learning methods on baby registry data.  See README.txt for
% an example usage.  Lines marked with "% *" here indicate internal
% parameters.

s = dir(sets_file);
if s.bytes == 0
  return;
end

% Read in the registries.
M = csvread(sets_file);
N = max(M(:));
num_ps = size(M, 1);

% Filter extreme cases.
if N > 750 || N < 25  % *
  fprintf('N = %d is too extreme.\n', N)
  return
end
num_train = floor(train_frac * num_ps);
if (num_train > 100 * N^2) || (num_train < 2 * N)  % *
  fprintf('# of examples (%d) is too different from N (%d).\n', num_ps, N);
  return
end

% Parse registries into an array.
Ps = cell(num_ps, 1);
for i = 1:num_ps
  end_idx = find(M(i, :) == 0, 1);
  if isempty(end_idx)
    Ps{i} = M(i, :);
  else
    Ps{i} = M(i, 1:end_idx-1);
  end
end
save(result_file, 'N', 'num_train', 'Ps', '-v7.3');

% Read in product names and save them in an array.
fid = fopen(names_file,'r');
names = textscan(fid,'%d %s', 'Delimiter', '\n');
fclose(fid);
save(result_file, 'names', '-append', '-v7.3');

% Set up for learning.
utils = opt_utils();
rutils = registries_utils();
eutils = eval_utils();

% Get learning method options.
em_opts = utils.get_em_opts(N);
ka_opts = utils.get_K_ascent_opts(N);

% Build train and test sets.
num_ps = numel(Ps);
num_train = floor(train_frac * num_ps);
assert(num_train > 0);
num_test = num_ps - num_train;
assert(num_test > 0);
train_ids = randsample(1:num_ps, num_train);
test_ids = setdiff(1:num_ps, train_ids);
T_train = rutils.build_training_set(Ps(train_ids), N);
T_train = utils.sort_training_set(T_train);
T_test = rutils.build_training_set(Ps(test_ids), N);
T_test = utils.sort_training_set(T_test);
% Note: We could also have sub-trials randomly sampling multiple train
%       sets here.  But our training sets are large enough that this
%       shouldn't have a large effect.

% Run the trials and collect a few basic statistics.
num_stats = 44;
results.stats_pg = zeros(1, num_stats, num_trials);
results.stats_eg = zeros(1, num_stats, num_trials);
results.stats_mm = zeros(1, num_stats, num_trials);
results.stats_diag = zeros(1, num_stats, num_trials);
for t = 1:num_trials
  % Handle the low-data setting.
  if max_numN_train > 0 && num_train > max_numN_train
    train_ids1 = randsample(train_ids, max_numN_train * N);
    T_train = rutils.build_training_set(Ps(train_ids1), N);
    T_train = utils.sort_training_set(T_train);
  end
  
  % Obtain an optimization starting point.
  max_eig = 0.999;
  max_L_eig = max_eig / (1 - max_eig);  % *
  if strcmp(start_type, 'mm')
     % Initialize K with data moments.
     K_start = utils.K_moments_init(T_train, max_eig);
  else
     if strcmp(start_type, 'wi')
       % Draw a random starting K from Wishart.
       K_start = utils.K_wishart_init(N, 1, max_L_eig);  % *
     else
       assert(false, 'Unknown start_type.  Try mm or wi.');
     end
  end
  results.Ks(1, t).start = K_start;
  
  % Run learning.
  tic;
  [K_em, obj_vals_em] = em(T_train, em_opts, utils, K_start);
  em_time = toc;
  results.Ks(1, t).em = K_em;
  results.obj_vals(1, t).em = obj_vals_em;
  [K_pg, obj_vals_pg, pg_time] = ...
      eutils.run_baseline(T_train, ka_opts, utils, K_start, 'pg');
  results.Ks(1, t).pg = K_pg;
  results.obj_vals(1, t).pg = obj_vals_pg;
  %[K_eg, obj_vals_eg, eg_time] = ...
  %    eutils.run_baseline(T_train, ka_opts, utils, K_start, 'eg');
  %results.Ks(1, t).eg = K_eg;
  %results.obj_vals(1, t).eg = obj_vals_eg;
  %[K_mm, obj_vals_mm, mm_time] = ...
  %    eutils.run_baseline(T_train, ka_opts, utils, K_start, 'mm');
  %results.Ks(1, t).mm = K_mm;
  %results.obj_vals(1, t).mm = obj_vals_mm;
  %[K_diag, obj_vals_diag, diag_time] = ...
  %    eutils.run_baseline(T_train, ka_opts, utils, K_start, 'diag');
  %results.Ks(1, t).diag = K_diag;
  %results.obj_vals(1, t).diag = obj_vals_diag;

  % Compare methods.
  % To re-use synth code we pass in K_baseline twice here;
  % be careful in the interpretation of the results!
  % Compare methods.
  results.stats_pg(1, :, t) = eutils.run_comparisons(K_pg, K_start, ...
      K_em, obj_vals_em, em_time, K_pg, obj_vals_pg, pg_time, T_train, ...
      T_test);
  %results.stats_eg(1, :, t) = eutils.run_comparisons(K_eg, K_start, ...
  %    K_em, obj_vals_em, em_time, K_eg, obj_vals_eg, eg_time, T_train, ...
  %    T_test);
  %results.stats_mm(1, :, t) = eutils.run_comparisons(K_mm, K_start, ...
  %    K_em, obj_vals_em, em_time, K_mm, obj_vals_mm, mm_time, T_train, ...
  %    T_test);
  %results.stats_diag(1, :, t) = eutils.run_comparisons(K_pg, K_start, ...
  %    K_em, obj_vals_em, em_time, K_diag, obj_vals_diag, diag_time, ...
  %    T_train, T_test);
end

% Compile summary statistics.
results.quantiles_pg = eutils.quantile_results(results.stats_pg);
%results.quantiles_eg = eutils.quantile_results(results.stats_eg);
%results.quantiles_mm = eutils.quantile_results(results.stats_mm);
%results.quantiles_diag = eutils.quantile_results(results.stats_diag);

save(result_file, 'results', '-append', '-v7.3');
