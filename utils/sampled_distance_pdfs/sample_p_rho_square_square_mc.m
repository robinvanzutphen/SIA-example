function [rho_cent, pdf_samp, src_plot, det_plot, edges, binw] = sample_p_rho_square_square_mc(s_src, s_det, dx, dy, N_pairs_total, N_plot, Nbins)
%SAMPLE_P_RHO_SQUARE_SQUARE_MC  Monte-Carlo sampled p(rho) for two
% axis-aligned squares with center offset [dx,dy].

rho_max  = hypot(dx + (s_src/2 + s_det/2), dy + (s_src/2 + s_det/2));
edges    = linspace(0, rho_max, Nbins+1);
binw     = edges(2) - edges(1);
rho_cent = (edges(1:end-1) + edges(2:end))/2;

Ps = sample_square_local(N_pairs_total, s_src);
Pd = sample_square_local(N_pairs_total, s_det) + [dx, dy];

rho = hypot(Ps(:,1) - Pd(:,1), Ps(:,2) - Pd(:,2));

counts   = histcounts(rho, edges);
pdf_samp = counts ./ (sum(counts) * binw + eps);

src_plot = Ps(1:N_plot, :);
det_plot = Pd(1:N_plot, :);
end

function P = sample_square_local(N, side_mm)
half = side_mm/2;
P = (rand(N,2) - 0.5) * (2*half);
end
