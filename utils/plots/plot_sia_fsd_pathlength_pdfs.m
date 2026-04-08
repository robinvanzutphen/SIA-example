function plot_sia_fsd_pathlength_pdfs(Cases, L_mid, pL_full_mean, pL_sia_mean)
%PLOT_SIA_FSD_PATHLENGTH_PDFS  Display detected pathlength PDFs p(L).

nC = numel(Cases);
figure('Color','w');
tiledlayout(nC,1,'TileSpacing','compact','Padding','compact');

for c = 1:nC
    ax = nexttile; hold(ax,'on'); box(ax,'on'); grid(ax,'on');

    plot(ax, L_mid, pL_full_mean(c,:), '-',  'LineWidth', 2.0, 'DisplayName','Full S-D');
    plot(ax, L_mid, pL_sia_mean (c,:), '--', 'LineWidth', 2.0, 'DisplayName','SIA');

    xlabel(ax,'Pathlength L [mm]');
    ylabel(ax,'p(L) (absorption-weighted, normalized)');
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
