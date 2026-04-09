function plot_sia_fsd_pathlength_pdfs(Cases, L_mid, pL_full, pL_sia)
% PLOT_SIA_FSD_PATHLENGTH_PDFS  Plot detected pathlength PDFs for FSD and SIA.
%
%  plot_sia_fsd_pathlength_pdfs(Cases, L_mid, pL_full, pL_sia)
%
%  For each selected scenario, plots the absorption-weighted detected
%  pathlength probability density function (PDF) obtained from the full
%  source-detector (FSD) simulation and the single-integral approximation
%  (SIA). The curves are pooled over all repetitions before normalisation
%  (see the main demo script Section 9), yielding smoother estimates than
%  averaging per-repetition PDFs.
%
%  INPUTS
%  ------
%  Cases    - (1 x n_cases) cell array of geometry structs (see build_sia_cases.m).
%               .display_name  - used as subplot title
%               .letter_group  - short identifier ('I'..'VI', 'extra')
%  L_mid    - (1 x n_L_bins) pathlength bin centres [mm]
%  pL_full  - (n_cases x n_L_bins) pooled FSD pathlength PDF [1/mm]
%  pL_sia   - (n_cases x n_L_bins) pooled SIA pathlength PDF [1/mm]

n_cases = numel(Cases);

% Layout: up to 4 columns, as many rows as needed.
n_cols = min(n_cases, 4);
n_rows = ceil(n_cases / n_cols);

% Full-screen named figure.
figure('Name', 'SIA vs FSD - Detected pathlength distributions', ...
       'NumberTitle', 'off', ...
       'WindowState', 'maximized');

% Line style settings.
col_full = [0.25 0.55 0.85];   % blue  = FSD
col_sia  = [0.90 0.35 0.20];   % red   = SIA
lw_full  = 1.8;                % FSD line width [pt]
lw_sia   = 1.2;                % SIA line width [pt], slightly thinner so FSD shows through

for c = 1:n_cases
    subplot(n_rows, n_cols, c);

    pL_f = pL_full(c, :);
    pL_s = pL_sia(c, :);

    if all(pL_f == 0) && all(pL_s == 0)
        text(0.5, 0.5, 'No data', 'HorizontalAlignment', 'center', ...
             'Units', 'normalized', 'FontSize', 9);
        title(Cases{c}.display_name, 'Interpreter', 'none');
        xlabel('Pathlength  L  [mm]');
        ylabel('p(L)  [mm^{-1}]');
        continue;
    end

    % Plot FSD (solid) and SIA (dashed) on the same axes.
    plot(L_mid, pL_f, '-',  'Color', col_full, 'LineWidth', lw_full, 'DisplayName', 'FSD');
    hold on;
    plot(L_mid, pL_s, '--', 'Color', col_sia,  'LineWidth', lw_sia,  'DisplayName', 'SIA');
    hold off;

    % y-axis: 0 to 110% of peak value for a small top margin.
    y_max = max([pL_f(:); pL_s(:)]);
    if y_max > 0
        ylim([0, 1.1 * y_max]);
    end
    xlim([0, L_mid(end)]);

    % Title shows the case label.
    title(Cases{c}.display_name, 'Interpreter', 'none', ...
          'FontWeight', 'normal', 'FontSize', 9);

    % Axis labels on every panel.
    xlabel('Pathlength  L  [mm]');
    ylabel('p(L)  [mm^{-1}]');

    legend({'FSD', 'SIA'}, 'Location', 'best', 'FontSize', 8);
    box off;
end

% Overall figure title.
sgtitle('Detected pathlength distributions: FSD vs SIA', ...
        'FontSize', 13, 'FontWeight', 'bold');

end % function plot_sia_fsd_pathlength_pdfs
