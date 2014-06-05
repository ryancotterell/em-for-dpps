function utils = registries_utils()
% Contains functionality shared between registries_test.m and
% registries_eval.m.

utils.build_training_set = @build_training_set;


function T = build_training_set(Ps, N)
% Record sizes for convenience.
T.N = N;
T.n = numel(Ps);

% Use a map from sets to counts to build up a de-duplicated training set.
import EqualByValueIntArray;
Ys_map = java.util.HashMap;
for i = 1:T.n
  % Sample a set.
  Y = EqualByValueIntArray(sort(Ps{i}));

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
