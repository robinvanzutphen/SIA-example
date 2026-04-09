function plot_sia_fsd_prho_overview(Cases, rho_plot_store, pdf_plot_samp_store, rho_theo_store, pdf_theo_store, theo_label_store)
%PLOT_SIA_FSD_PRHO_OVERVIEW  Sampled-vs-analytic distance PDF overview.

nC = numel(Cases);

% Layout: up to 4 columns, as many rows as needed.
n_cols = min(nC, 4);
n_rows = ceil(nC / n_cols);

figure('Name', 'SIA vs FSD - Source-detector distance distributions', ...
       'NumberTitle', 'off', ...
       'Color', 'w', ...
       'WindowState', 'maximized');

for c = 1:nC
    ax = subplot(n_rows, n_cols, c);
    hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on');

    plot(ax, rho_plot_store{c}, pdf_plot_samp_store{c}, '-', ...
         'LineWidth', 2.2, 'DisplayName', 'Sampled');

    Ptheo = pdf_theo_store{c};
    if ~isempty(Ptheo)
        plot(ax, rho_theo_store{c}, Ptheo, '--', 'LineWidth', 2.2, ...
             'DisplayName', sprintf('Analytical (%s)', theo_label_store{c}));
    end

    xlabel(ax, '\rho  [mm]');
    ylabel(ax, 'p(\rho)  [mm^{-1}]');
    title(ax, get_case_title(Cases{c}), 'Interpreter', 'none', ...
          'FontWeight', 'normal', 'FontSize', 9);

    legend(ax, 'Location', 'best', 'FontSize', 8);
    box off;
end

sgtitle('Source–detector distance distributions: sampled vs analytical', ...
        'FontSize', 13, 'FontWeight', 'bold');
end

% -------------------------------------------------------------------------
function title_str = get_case_title(C)
if     isfield(C, 'display_name'),  title_str = C.display_name;
elseif isfield(C, 'name'),          title_str = C.name;
else,                                title_str = 'Case';
end
end
