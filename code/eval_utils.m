function utils = eval_utils()
% Contains evaluation functionality common to synthetic- and
% real-data evaluations.  See synth_test.m, registries_test.m,
% synth_eval.m, and registries_eval.m for example uses.

utils.compare_learned_Ks = @compare_learned_Ks;
utils.compare_learned_Ks_likelihoods = @compare_learned_Ks_likelihoods;

utils.eval_target_type = @eval_target_type;


function [kls, kl_diffs] = compare_learned_Ks(K, K_start, K_em, K_ka)
N = size(K, 1);
num_kl_samples = 10 * N;
kl_div = @(K1, K2) dpp_kl_div(K1, K2, num_kl_samples);

em_true_Kkl = kl_div(K, K_em);
em_true_Kkl_diag = kl_div(K, decompose(diag(diag(K_em.M))));
em_start_Kkl = kl_div(K_start, K_em);
ka_true_Kkl = kl_div(K, K_ka);
ka_true_Kkl_diag = kl_div(K, decompose(diag(diag(K_ka.M))));
ka_start_Kkl = kl_div(K_start, K_ka);
start_true_Kkl = kl_div(K, K_start);
true_Kkl_diag = kl_div(K, decompose(diag(diag(K.M))));

kls = [em_true_Kkl, ...
  em_true_Kkl_diag, ...
  em_start_Kkl, ...
  ka_true_Kkl, ...
  ka_true_Kkl_diag, ...
  ka_start_Kkl, ...
  start_true_Kkl, ...
  true_Kkl_diag];

kl_diffs = [em_true_Kkl - ka_true_Kkl, ...
  em_true_Kkl - start_true_Kkl, ...
  em_true_Kkl - em_true_Kkl_diag, ...
  ka_true_Kkl - start_true_Kkl, ...
  ka_true_Kkl - ka_true_Kkl_diag];
  

function [lls, ll_diffs] = ...
  compare_learned_Ks_likelihoods(K, K_start, K_em, K_ka, T_train, T_test)
utils = opt_utils();

ll_em = utils.K_log_likelihood(T_train, K_em.M);
ll_test_em = utils.K_log_likelihood(T_test, K_em.M);
ll_test_em_diag = utils.K_log_likelihood(T_test, diag(diag(K_em.M)));
ll_ka = utils.K_log_likelihood(T_train, K_ka.M);
ll_test_ka = utils.K_log_likelihood(T_test, K_ka.M);
ll_test_ka_diag = utils.K_log_likelihood(T_test, diag(diag(K_ka.M)));
ll_start = utils.K_log_likelihood(T_train, K_start.M);
ll_test_start = utils.K_log_likelihood(T_test, K_start.M);
ll_test_start_diag = utils.K_log_likelihood(T_test, diag(diag(K_start.M)));
ll_true = utils.K_log_likelihood(T_train, K.M);
ll_test_true = utils.K_log_likelihood(T_test, K.M);
ll_test_diag = utils.K_log_likelihood(T_test, diag(diag(K.M)));

lls = [ll_em, ...
  ll_test_em, ...
  ll_test_em_diag, ... 
  ll_ka, ... 
  ll_test_ka, ...
  ll_test_ka_diag, ...
  ll_start, ...
  ll_test_start, ...
  ll_test_start_diag, ...
  ll_true, ...
  ll_test_true, ...
  ll_test_diag];

ll_diffs = [(ll_em - ll_ka) / abs(ll_true), ...
  (ll_em - ll_start) / abs(ll_true), ...
  (ll_ka - ll_start) / abs(ll_true), ...
  (ll_test_em - ll_test_ka) / abs(ll_test_true), ...
  (ll_test_em - ll_test_start) / abs(ll_test_true), ...
  (ll_test_em - ll_test_em_diag) / abs(ll_test_true), ...
  (ll_test_ka - ll_test_start) / abs(ll_test_true), ...
  (ll_test_ka - ll_test_ka_diag) / abs(ll_test_true), ...
  (ll_test_em - ll_test_true) / abs(ll_test_true), ...
  (ll_test_em - ll_test_diag) / abs(ll_test_true), ...
  (ll_test_true - ll_test_diag) / abs(ll_test_true)];

  
function eval_target_type(r, T_sizes, indices, sum_str)
T_num = 0;
for T_size = T_sizes
  T_num = T_num + 1;
  fstr = [' (%.4f) %.4f (%.4f) ---- ' sprintf('T_size = %d; ', T_size) ...
    sum_str '\n'];
  out_strs = {['em_true_Kkl' fstr], ... % 1
    ['em_true_Kkl_diag' fstr], ...
    ['em_start_Kkl' fstr], ...
    ['ka_true_Kkl' fstr], ...
    ['ka_true_Kkl_diag' fstr], ...
    ['ka_start_Kkl' fstr], ...
    ['start_true_Kkl' fstr], ...
    ['true_Kkl_diag' fstr], ...
    ['ll_em' fstr], ... % 9
    ['ll_test_em' fstr], ...
    ['ll_test_em_diag' fstr], ...
    ['ll_ka' fstr], ...
    ['ll_test_ka' fstr], ...
    ['ll_test_ka_diag' fstr], ...
    ['ll_start' fstr], ...
    ['ll_test_start' fstr], ...
    ['ll_test_start_diag' fstr], ...
    ['ll_true' fstr], ...
    ['ll_test_true' fstr], ...
    ['ll_test_diag' fstr], ...
    ['em_time' fstr], ... % 21
    ['ka_time' fstr], ...
    ['** (ll_em - ll_ka) / abs(ll_true)' fstr], ... % 23
    ['(ll_em - ll_start) / abs(ll_true)' fstr], ...
    ['(ll_ka - ll_start) / abs(ll_true)' fstr], ...
    ['** (ll_test_em - ll_test_ka) / abs(ll_test_true)' fstr], ...
    ['! (ll_test_em - ll_test_start) / abs(ll_test_true)' fstr], ...
    ['! (ll_test_em - ll_test_em_diag) / abs(ll_test_true)' fstr], ...
    ['! (ll_test_ka - ll_test_start) / abs(ll_test_true)' fstr], ...
    ['! (ll_test_ka - ll_test_ka_diag) / abs(ll_test_true)' fstr], ...
    ['(ll_test_em - ll_test_true) / abs(ll_test_true)' fstr], ...
    ['(ll_test_em - ll_test_diag) / abs(ll_test_true)' fstr], ...
    ['** (ll_test_true - ll_test_diag) / abs(ll_test_true)' fstr], ...
    ['** em_true_Kkl - ka_true_Kkl' fstr], ... % 34
    ['! em_true_Kkl - start_true_Kkl' fstr], ...
    ['! em_true_Kkl - em_true_Kkl_diag' fstr], ...
    ['! ka_true_Kkl - start_true_Kkl' fstr], ...
    ['! ka_true_Kkl - ka_true_Kkl_diag' fstr], ...
    ['(em_time - ka_time) / em_time' fstr]}; % 39

  fprintf([out_strs{indices} '\n'], ...
  reshape(r.quantiles(T_num, :, :), [], 3)');
end