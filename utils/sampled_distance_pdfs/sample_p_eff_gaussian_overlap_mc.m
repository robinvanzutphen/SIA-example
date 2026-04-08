function [rho_mid, p_samp, src_plot, det_plot, rho_edges, drho] = sample_p_eff_gaussian_overlap_mc(w_mm, R_s, R_d, N_pairs_total, N_plot, Nbins)
%SAMPLE_P_EFF_GAUSSIAN_OVERLAP_MC  Weighted Monte-Carlo sampled p_eff(rho)
% for overlapping, truncated Gaussian source/detector weights.

if any([w_mm, R_s, R_d] <= 0)
    error('w_mm, R_s, R_d must be > 0');
end

rho_max   = R_s + R_d;
rho_edges = linspace(0, rho_max, Nbins+1);
drho      = rho_edges(2) - rho_edges(1);
rho_mid   = (rho_edges(1:end-1) + rho_edges(2:end))/2;
rho_mid   = rho_mid(:);

inv_r_gauss_trunc = @(u, Rtr) (w_mm/sqrt(2)) * sqrt( ...
    -log( 1 - u*(1 - exp(-2*Rtr^2/w_mm^2)) ) );

u1  = rand(N_pairs_total,1);  th1 = 2*pi*rand(N_pairs_total,1);
u2  = rand(N_pairs_total,1);  th2 = 2*pi*rand(N_pairs_total,1);

rs  = inv_r_gauss_trunc(u1, R_s);
rd  = inv_r_gauss_trunc(u2, R_d);

xs  = rs .* cos(th1);  ys = rs .* sin(th1);
xd  = rd .* cos(th2);  yd = rd .* sin(th2);

rho_samp = hypot(xs - xd, ys - yd);

cnt    = histcounts(rho_samp, rho_edges);
p_samp = cnt ./ (sum(cnt)*drho + eps);
p_samp = p_samp(:);

src_plot = [xs(1:N_plot), ys(1:N_plot)];
det_plot = [xd(1:N_plot), yd(1:N_plot)];
end
