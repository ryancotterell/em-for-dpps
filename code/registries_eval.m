function registries_eval(result_file)
% Prints information concerning the results from registries_test.m.

data = load(result_file);
if ~isfield(data, 'results')
  return;
end
r = data.results;
N = data.N;
num_train = data.num_train;
Ps = data.Ps;
names = data.names;

utils = opt_utils();
rutils = registries_utils();
eutils = eval_utils();

% Print basic EM-vs-baseline comparison statistics.
T = rutils.build_training_set(Ps(:), N);
k = median(T.Y_sizes);
num_stats = size(r.quantiles_pg, 2);

sum_str = sprintf('N = %d; k = %d; baseline = %s', N, k, 'pg');
r.quantiles = r.quantiles_pg;
eutils.eval_target_type(r, num_train, 1:num_stats, sum_str);

%sum_str = sprintf('N = %d; k = %d; baseline = %s', N, k, 'eg');
%r.quantiles = r.quantiles_eg;
%eutils.eval_target_type(r, num_train, 1:num_stats, sum_str);

%sum_str = sprintf('N = %d; k = %d; baseline = %s', N, k, 'mm');
%r.quantiles = r.quantiles_mm;
%eutils.eval_target_type(r, num_train, 1:num_stats, sum_str);

%sum_str = sprintf('N = %d; k = %d; baseline = %s', N, k, 'diag');
%r.quantiles = r.quantiles_diag;
%eutils.eval_target_type(r, num_train, 1:num_stats, sum_str);

% Check how important the off-diagonal is for this dataset.
% (Larger numbers here mean the off-diag is more important.)
Q = utils.get_low_moments(T);
diag_stren = (N / N^2) * norm(Q, 'fro') / norm(diag(Q), 2);
fprintf('Off-diagonal strength: %.5f\n', diag_stren);
