function p = p_rho_core_total_analytic(rho, r_core, R_total)
%P_RHO_CORE_TOTAL_ANALYTIC  Analytical distance PDF p(rho) for concentric disks:
%   - point A uniform in disk radius r_core
%   - point B uniform in disk radius R_total (R_total >= r_core)

rho_in = rho;
rho    = double(rho(:));
a      = double(r_core);
R      = double(R_total);

p = zeros(size(rho));

in = (rho >= 0) & (rho <= (R + a));
if ~any(in)
    p = reshape(p, size(rho_in));
    return;
end

r = rho(in);

r1 = (r > 0) & (r <= (R - a));
p_in = zeros(size(r));

p_in(r1) = (2 .* r(r1)) ./ (R.^2);

r2 = (r > (R - a)) & (r <= (R + a)) & (r > 0);
if any(r2)
    rr = r(r2);

    c1 = (a^2 + rr.^2 - R^2) ./ (2*a.*rr);
    c2 = (R^2 + rr.^2 - a^2) ./ (2*R.*rr);

    c1 = min(max(c1, -1), 1);
    c2 = min(max(c2, -1), 1);

    alpha = acos(c1);
    beta  = acos(c2);

    p_in(r2) = (rr/pi) .* ( ...
        (2*alpha - sin(2*alpha))./(R.^2) + ...
        (2*beta  - sin(2*beta ))./(a.^2) );
end

p(in) = p_in;

p(~isfinite(p)) = 0;
p(p < 0) = 0;

p = reshape(p, size(rho_in));
end
