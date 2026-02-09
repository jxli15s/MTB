function plot_band_path(g, k_nodes_frac, n_per_seg, labels)
%PLOT_BAND_PATH  Build path, solve E(k), plot with labels/ticks.

% cache hops
g = tbHFMF.prep_tb_cache(g);

% path
path = tbHFMF.build_kpath(g, k_nodes_frac, n_per_seg, labels);

% solve
out = tbHFMF.solve_band_path_E(g, path.klist_frac);

% plot
figure;
plot(path.kdist, out.E, 'LineWidth', 1.0);
grid on; box on;
xlabel('k-path distance (Å^{-1})');
ylabel('Energy (eV)');

% ticks + labels
if ~isempty(labels)
  xticks(path.kdist(path.tick_idx));
  xticklabels(labels);
end

% vertical separators at nodes
xline(path.kdist(path.tick_idx), ':');
set(gca,'TickLabelInterpreter','latex');
end
