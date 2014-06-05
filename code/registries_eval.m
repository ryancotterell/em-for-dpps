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

% Print basic EM-KA comparison statistics.
T = rutils.build_training_set(Ps(:), N);
k = median(T.Y_sizes);
sum_str = sprintf('N = %d; k = %d', N, k);
eutils.eval_target_type(r, num_train, [9:20, 23:33], sum_str)

% Check how important the off-diagonal is for this dataset.
% (Smaller numbers here mean the off-diag is more important.)
Q = utils.get_low_moments(T);
Q_od = Q - diag(diag(Q));
diag_stren = (N / N^2) * norm(Q, 'fro') / norm(diag(Q), 2);
fprintf('Diagonal strength: %.5f\n', diag_stren);
