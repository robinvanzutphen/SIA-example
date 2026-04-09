function plot_sia_fsd_reflectance_hists(Cases, R_full, R_sia)
% PLOT_SIA_FSD_REFLECTANCE_HISTS  Plot overlapping reflectance distributions
%                                   for the FSD and SIA approaches.
%
%  plot_sia_fsd_reflectance_hists(Cases, R_full, R_sia)
%
%  For each selected scenario, plots overlapping histograms of the reflectance
%  values from the full source-detector (FSD) and single-integral approximation
%  (SIA) Monte Carlo repetitions. Close overlap demonstrates equivalence.
%
%  INPUTS
%  ------
%  Cases   - (1 x n_cases) cell array of geometry structs (see build_sia_cases.m).
%              .display_name  - used as subplot title
%              .letter_group  - short identifier ('I'..'VI', 'extra')
%  R_full  - (n_cases x n_repetitions) FSD reflectance values [-]
%  R_sia   - (n_cases x n_repetitions) SIA reflectance values [-]

n_cases = numel(Cases);

% Layout: up to 4 columns, as many rows as needed.
n_cols = min(n_cases, 4);
n_rows = ceil(n_cases / n_cols);

% Full-screen named figure.
figure('Name', 'SIA vs FSD - Reflectance distributions', ...
       'NumberTitle', 'off', ...
       'WindowState', 'maximized');

% Histogram appearance.
n_bins_hist = 20;
alpha_val   = 0.55;
col_full    = [0.25 0.55 0.85];   % blue  = FSD
col_sia     = [0.90 0.35 0.20];   % red   = SIA

for c = 1:n_cases
    subplot(n_rows, n_cols, c);

    r_full_c = R_full(c, isfinite(R_full(c,:)));
    r_sia_c  = R_sia(c,  isfinite(R_sia(c,:)));

    if isempty(r_full_c) || isempty(r_sia_c)
        text(0.5, 0.5, 'No data', 'HorizontalAlignment', 'center', ...
             'Units', 'normalized', 'FontSize', 9);
        title(Cases{c}.display_name, 'Interpreter', 'none');
        xlabel('Reflectance R  [-]');
        ylabel('Count  [-]');
        continue;
    end

    % Common bin edges with a small margin so edge bars are fully visible.
    r_min     = min([r_full_c(:); r_sia_c(:)]);
    r_max     = max([r_full_c(:); r_sia_c(:)]);
    margin    = 0.05 * (r_max - r_min + eps);
    bin_edges = linspace(r_min - margin, r_max + margin, n_bins_hist + 1);

    % Overlapping histograms.
    histogram(r_full_c, bin_edges, 'FaceColor', col_full, 'FaceAlpha', alpha_val, ...
              'EdgeColor', 'none', 'DisplayName', 'FSD');
    hold on;
    histogram(r_sia_c,  bin_edges, 'FaceColor', col_sia,  'FaceAlpha', alpha_val, ...
              'EdgeColor', 'none', 'DisplayName', 'SIA');
    hold off;

    % Percentage difference between the two distribution means.
    delta_pct = 100 * (mean(r_sia_c) - mean(r_full_c)) / max(abs(mean(r_full_c)), eps);

    % Title shows the case label and the mean reflectance difference.
    title(sprintf('%s  |  \\DeltaR/R = %+.3f%%', Cases{c}.display_name, delta_pct), ...
          'Interpreter', 'tex', 'FontWeight', 'normal', 'FontSize', 9);

    % Axis labels on every panel.
    xlabel('Reflectance R  [-]');
    ylabel('Count  [-]');

    legend({'FSD', 'SIA'}, 'Location', 'best', 'FontSize', 8);
    box off;
end

% Add an overall figure title.
sgtitle('Reflectance distributions: FSD vs SIA', 'FontSize', 13, 'FontWeight', 'bold');

end % function plot_sia_fsd_reflectance_hists
