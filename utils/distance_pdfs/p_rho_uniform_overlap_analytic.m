function p = p_rho_uniform_overlap_analytic(rho, diam)
%P_RHO_UNIFORM_OVERLAP_ANALYTIC  Distance PDF p(rho) for two independent
% uniform points in the same disk of diameter "diam".
%
%   p = p_rho_uniform_overlap_analytic(rho, diam)
%
% Inputs
%   rho  : vector of distances at which to evaluate p(rho)
%   diam : disk diameter (same units as rho)
%
% Output
%   p    : p(rho) evaluated at rho (same size as rho)
%
% Closed-form result (support 0 <= rho <= diam):
%   p(rho) = (16*rho)/(pi*diam^2) * acos(rho/diam)
%            - (16*rho^2)/(pi*diam^3) * sqrt(1 - (rho/diam)^2)

rho_in = rho;
rho    = double(rho(:));
d      = double(diam);

p = zeros(size(rho));

in = (rho >= 0) & (rho <= d);
if any(in)
    x = rho(in) ./ d;
    x = min(max(x, 0), 1);

    p(in) = (16.*rho(in))./(pi*d.^2).*acos(x) ...
          - (16.*rho(in).^2)./(pi*d.^3).*sqrt(max(0, 1 - x.^2));
end

p(~isfinite(p)) = 0;
p(p < 0) = 0;

p = reshape(p, size(rho_in));
end
