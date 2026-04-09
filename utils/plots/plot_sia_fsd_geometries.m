function plot_sia_fsd_geometries(Cases, geom_src_pts, geom_det_pts)
%PLOT_SIA_FSD_GEOMETRIES  Geometry overview with sampled source/detector points.

nC = numel(Cases);

n_cols = min(nC, 4);
n_rows = ceil(nC / n_cols);

figure('Name', 'SIA vs FSD - Source/detector geometries', ...
       'NumberTitle', 'off', ...
       'Color', 'w', ...
       'WindowState', 'maximized');

for c = 1:nC
    ax = subplot(n_rows, n_cols, c);
    hold(ax, 'on'); box(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');

    Ps = geom_src_pts{c};
    Pd = geom_det_pts{c};

    Ns = min(500, size(Ps,1));
    Nd = min(500, size(Pd,1));

    s1 = scatter(ax, Ps(1:Ns,1), Ps(1:Ns,2), 6, 'filled', 'DisplayName', 'Source');
    s2 = scatter(ax, Pd(1:Nd,1), Pd(1:Nd,2), 6, 'filled', 'DisplayName', 'Detector');

    draw_case_outline(ax, Cases{c});

    xlabel(ax, 'x  [mm]');
    ylabel(ax, 'y  [mm]');
    title(ax, get_case_title(Cases{c}), 'Interpreter', 'none', ...
          'FontWeight', 'normal', 'FontSize', 9);

    lim = case_plot_limit(Cases{c}, Ps, Pd);
    xlim(ax, [-lim lim]);
    ylim(ax, [-lim lim]);

    legend(ax, [s1 s2], {'Source', 'Detector'}, 'Location', 'best', 'FontSize', 8);
end

sgtitle('Source–detector geometry overview', 'FontSize', 13, 'FontWeight', 'bold');
end

% -------------------------------------------------------------------------
function draw_case_outline(ax, C)
t = linspace(0, 2*pi, 512);
C = normalize_case_fields(C);

switch C.type
    case 'uniform_shapes'
        draw_uniform_shape(ax, [0 0],           C.src_shape, C.src_rin, C.src_rout, C.src_side, 'k--');
        draw_uniform_shape(ax, C.det_offset_mm, C.det_shape, C.det_rin, C.det_rout, C.det_side, 'k-');

    case 'square_square'
        draw_square_outline(ax, 0, 0, C.src_side, 'k--');
        draw_square_outline(ax, C.det_offset_mm(1), C.det_offset_mm(2), C.det_side, 'k-');

    case 'gaussian_overlap'
        Rsrc_cap = C.R_det_mm + C.Nw_cap * C.w_mm;
        plot(ax, Rsrc_cap*cos(t), Rsrc_cap*sin(t), 'k--', 'LineWidth', 1.0);
        plot(ax, C.R_det_mm*cos(t), C.R_det_mm*sin(t), 'k-', 'LineWidth', 1.0);

    otherwise
        error('Unknown case type.');
end
end

% -------------------------------------------------------------------------
function draw_uniform_shape(ax, ctr, shape, rin, rout, side, ls)
t = linspace(0, 2*pi, 512);
switch shape
    case 'disk'
        plot(ax, ctr(1)+rout*cos(t), ctr(2)+rout*sin(t), ls, 'LineWidth', 1.0);
    case 'annulus'
        plot(ax, ctr(1)+rout*cos(t), ctr(2)+rout*sin(t), ls, 'LineWidth', 1.0);
        plot(ax, ctr(1)+ rin*cos(t), ctr(2)+ rin*sin(t), ls, 'LineWidth', 1.0);
    case 'square'
        draw_square_outline(ax, ctr(1), ctr(2), side, ls);
    otherwise
        error('Unknown shape.');
end
end

% -------------------------------------------------------------------------
function draw_square_outline(ax, cx, cy, side_mm, ls)
h = side_mm / 2;
x = cx + [-h, +h, +h, -h, -h];
y = cy + [-h, -h, +h, +h, -h];
plot(ax, x, y, ls, 'LineWidth', 1.0);
end

% -------------------------------------------------------------------------
function lim = case_plot_limit(C, Ps, Pd)
C = normalize_case_fields(C);

switch C.type
    case 'uniform_shapes'
        switch C.src_shape
            case {'disk', 'annulus'},  rs = C.src_rout;
            otherwise,                 rs = 0.5 * C.src_side;
        end
        switch C.det_shape
            case {'disk', 'annulus'},  rd = norm(C.det_offset_mm) + C.det_rout;
            otherwise,                 rd = norm(C.det_offset_mm) + 0.5 * C.det_side;
        end
        lim = 1.10 * max([rs, rd, max(hypot(Ps(:,1),Ps(:,2))), max(hypot(Pd(:,1),Pd(:,2))), 1e-3]);

    case 'square_square'
        xmin = min(-C.src_side/2,  C.det_offset_mm(1) - C.det_side/2);
        xmax = max( C.src_side/2,  C.det_offset_mm(1) + C.det_side/2);
        ymin = min(-C.src_side/2,  C.det_offset_mm(2) - C.det_side/2);
        ymax = max( C.src_side/2,  C.det_offset_mm(2) + C.det_side/2);
        lim  = max(abs([xmin xmax ymin ymax])) * 1.15;

    case 'gaussian_overlap'
        lim = 1.0;   % fixed ±1 mm window for the Gaussian case

    otherwise
        lim = 1.1 * max([max(abs(Ps(:))), max(abs(Pd(:))), 1]);
end
end

% -------------------------------------------------------------------------
function title_str = get_case_title(C)
if     isfield(C, 'display_name'),  title_str = C.display_name;
elseif isfield(C, 'name'),          title_str = C.name;
else,                                title_str = 'Case';
end
end

% -------------------------------------------------------------------------
function C = normalize_case_fields(C)
if isfield(C,'source_shape')       && ~isfield(C,'src_shape'),      C.src_shape     = C.source_shape; end
if isfield(C,'detector_shape')     && ~isfield(C,'det_shape'),      C.det_shape     = C.detector_shape; end
if isfield(C,'source_rin')         && ~isfield(C,'src_rin'),        C.src_rin       = C.source_rin; end
if isfield(C,'source_rout')        && ~isfield(C,'src_rout'),       C.src_rout      = C.source_rout; end
if isfield(C,'source_side')        && ~isfield(C,'src_side'),       C.src_side      = C.source_side; end
if isfield(C,'detector_rin')       && ~isfield(C,'det_rin'),        C.det_rin       = C.detector_rin; end
if isfield(C,'detector_rout')      && ~isfield(C,'det_rout'),       C.det_rout      = C.detector_rout; end
if isfield(C,'detector_side')      && ~isfield(C,'det_side'),       C.det_side      = C.detector_side; end
if isfield(C,'detector_offset_mm') && ~isfield(C,'det_offset_mm'), C.det_offset_mm = C.detector_offset_mm; end
end
