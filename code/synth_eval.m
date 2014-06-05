function synth_eval(result_file)
% Prints information concerning the results from synth_test.m.

data = load(result_file);
if ~isfield(data, 'results')
  return;
end
r = data.results;
Ns = data.Ns;

for N_num = 1:numel(Ns)
  N = Ns(N_num);
  sum_str = sprintf('N = %d; ', N);
  
  rN = r(N_num);
  T_sizes = r(N_num).T_sizes;
  
  if isfield(rN, 'mm_start')
    sum_str = [sum_str sprintf('@ Moment matching start; ')];
    eval_all_target_types(rN.mm_start, T_sizes, sum_str);
  end
  
  if isfield(rN, 'mc_start')
    min_Y_sizes = rN.min_Y_sizes;
    base_sum_str = sum_str;
    for i = 1:numel(min_Y_sizes)
      sum_str = [base_sum_str ...
        sprintf('@ Truncated moments start with cutoff: %f; ', ...
        min_Y_sizes(i))];
      eval_all_target_types(rN.mc_start(i), T_sizes, sum_str);
    end
  end
  
  if isfield(rN, 'wi_start')
    L_scalings = rN.L_scalings;
    base_sum_str = sum_str;
    for i = 1:numel(L_scalings)
      sum_str = [base_sum_str ...
        sprintf('@ Wishart start with scaling: %f; ', ...
        L_scalings(i))];
      eval_all_target_types(rN.wi_start(i), T_sizes, sum_str);
    end
  end
  
  if isfield(r(N_num), 'Dshaken_start')
    shake_strengths = rN.shake_strengths;
    base_sum_str = sum_str;
    for i = 1:numel(shake_strengths)
      sum_str = [base_sum_str ...
        sprintf('@ Eigenvalue perturbation start: %f; ', ...
        shake_strengths(i))];
      eval_all_target_types(rN.Dshaken_start(i), T_sizes, sum_str);
    end    
  end
  
  if isfield(rN, 'Vrotated_start')
    rotation_strengths = rN.rotation_strengths;
    base_sum_str = sum_str;
    for i = 1:numel(rotation_strengths)
      sum_str = [base_sum_str ...
        sprintf('@ Eigenvector perturbation start: %f; ', ...
        rotation_strengths(i))];
      eval_all_target_types(rN.Vrotated_start(i), T_sizes, sum_str);
    end 
  end
end


function eval_all_target_types(r, T_sizes, sum_str)
eutils = eval_utils();
num_stats = size(r.wishart_tt.results(1).quantiles, 2);
indices = 1:num_stats;
base_sum_str = sum_str;

L_scalings = r.wishart_tt.L_scalings;
for i = 1:numel(L_scalings)
  sum_str = [base_sum_str ...
    sprintf('@ Wishart target with scaling: %f', L_scalings(i))];
  eutils.eval_target_type(r.wishart_tt.results(i), T_sizes, indices, ...
    sum_str);
end

cluster_closenesses = r.cluster_tt.cluster_closenesses;
kernel_closenesses = r.cluster_tt.kernel_closenesses;
for i = 1:numel(cluster_closenesses)
  for j = 1:numel(kernel_closenesses)
    sum_str = [base_sum_str, ...
      sprintf('@ Cluster target with parameters: %f, %f', ...
      cluster_closenesses(i), kernel_closenesses(j))];
    eutils.eval_target_type(r.cluster_tt.results(i, j), T_sizes, ...
      indices, sum_str);
  end
end
