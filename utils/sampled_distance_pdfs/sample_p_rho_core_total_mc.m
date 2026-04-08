function [rho_cent, pdf_samp, src_plot, det_plot, edges, binw] = sample_p_rho_core_total_mc(r_core, R_total, N_pairs_total, N_plot, Nbins)
%SAMPLE_P_RHO_CORE_TOTAL_MC  Monte-Carlo sampled p(rho) for one point in a
% core disk and one point in a larger concentric total disk.

rho_max  = R_total + r_core;
edges    = linspace(0, rho_max, Nbins+1);
binw     = edges(2) - edges(1);
rho_cent = (edges(1:end-1) + edges(2:end))/2;

th_s = 2*pi*rand(N_pairs_total,1);
r_s  = r_core*sqrt(rand(N_pairs_total,1));
xs   = r_s .* cos(th_s);
ys   = r_s .* sin(th_s);

th_d = 2*pi*rand(N_pairs_total,1);
r_d  = R_total*sqrt(rand(N_pairs_total,1));
xd   = r_d .* cos(th_d);
yd   = r_d .* sin(th_d);

rho = hypot(xs - xd, ys - yd);

counts   = histcounts(rho, edges);
pdf_samp = counts ./ (sum(counts) * binw + eps);

src_plot = [xs(1:N_plot), ys(1:N_plot)];
det_plot = [xd(1:N_plot), yd(1:N_plot)];
end
