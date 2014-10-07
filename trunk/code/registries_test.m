function registries_test(sets_file, names_file, result_file, ...
			 train_frac, num_trials, start_type)
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
num_stats = 23;
results.stats = zeros(1, num_stats, num_trials);
for t = 1:num_trials
  % Draw a random starting K from Wishart.
  max_eig = 0.99;
  max_L_eig = max_eig / (1 - max_eig);  % *
  if strcmp(start_type, 'mm')
     K_start = utils.K_moments_init(T_train, max_eig);
  else
     if strcmp(start_type, 'wi')
       K_start = utils.K_wishart_init(N, 1, max_L_eig);  % *
     else
       assert(false, 'Unknown start_type.  Try mm or wi.');
     end
   end
  
  % Run learning.
  K_em = em(T_train, em_opts, utils, K_start);
  K_ka = K_ascent(T_train, ka_opts, utils, K_start, 'pg');

  % Compare learned Ks' likelihoods.
  % To re-use synth code we pass in K_em twice here;
  % be careful in the interpretation of the results!
  [lls, ll_diffs] = ...
    eutils.compare_learned_Ks_likelihoods(K_ka, K_start, K_em, K_ka, ...
      T_train, T_test);
   results.stats(1, :, t) = [lls, ll_diffs];
   
   % Save trial data.
   results.Ks(t).start = K_start;
   results.Ks(t).em = K_em;
   results.Ks(t).ka = K_ka;
end

% Compile summary statistics.
results.quantiles(1, :, :) = ...
    prctile(reshape(results.stats(1, :, :), num_stats, num_trials), ...
      [25, 50, 75], 2);

save(result_file, 'results', '-append', '-v7.3');
