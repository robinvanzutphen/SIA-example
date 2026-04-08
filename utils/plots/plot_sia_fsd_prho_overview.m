function plot_sia_fsd_prho_overview(Cases, rho_plot_store, pdf_plot_samp_store, rho_theo_store, pdf_theo_store, theo_label_store)
%PLOT_SIA_FSD_PRHO_OVERVIEW  Sampled-vs-analytic distance PDF overview.

nC = numel(Cases);
figure('Color','w');
tiledlayout(nC,1,'TileSpacing','compact','Padding','compact');

for c = 1:nC
    ax = nexttile; hold(ax,'on'); box(ax,'on'); grid(ax,'on');

    rho_plot = rho_plot_store{c};
    pdf_plot = pdf_plot_samp_store{c};

    plot(ax, rho_plot, pdf_plot, '-', 'LineWidth', 2.2, 'DisplayName','Sampled');

    Ptheo = pdf_theo_store{c};
    if ~isempty(Ptheo)
        plot(ax, rho_theo_store{c}, Ptheo, '--', 'LineWidth', 2.2, ...
             'DisplayName', sprintf('Theoretical (%s)', theo_label_store{c}));
    end

    xlabel(ax,'\rho [mm]');
    ylabel(ax,'p(\rho)');
    title(ax, get_case_title(Cases{c}), 'Interpreter','none');

    if c == 1
        legend(ax,'Location','best');
    end
end
end

function title_str = get_case_title(C)
if isfield(C,'display_name')
    title_str = C.display_name;
elseif isfield(C,'name')
    title_str = C.name;
else
    title_str = 'Case';
end
end
