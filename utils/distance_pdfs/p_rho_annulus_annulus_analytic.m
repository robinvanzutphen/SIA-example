function p = p_rho_annulus_annulus_analytic(rho, Rs_in, Rs_out, Rd_in, Rd_out)
%P_RHO_ANNULUS_ANNULUS_ANALYTIC  Analytic p(rho) for annulus -> annulus (concentric)

rho_in = rho;
rho    = double(rho(:));
p      = zeros(size(rho));

A_S = pi * (Rs_out^2 - Rs_in^2);
A_D = pi * (Rd_out^2 - Rd_in^2);

if A_S <= 0 || A_D <= 0
    error('Annulus areas must be positive. Check inner/outer radii.');
end

Acap = circle_overlap_area_local(rho, Rs_out, Rd_out) ...
     - circle_overlap_area_local(rho, Rs_out, Rd_in ) ...
     - circle_overlap_area_local(rho, Rs_in , Rd_out) ...
     + circle_overlap_area_local(rho, Rs_in , Rd_in );

p = (2*pi .* rho ./ (A_S * A_D)) .* Acap;
p(p < 0 & p > -1e-14) = 0;
p(~isfinite(p)) = 0;

p = reshape(p, size(rho_in));
end

function A = circle_overlap_area_local(d, R1, R2)
d_in = d;
d  = double(d(:));
R1 = double(R1);
R2 = double(R2);

A = zeros(size(d));

if R1 < 0 || R2 < 0
    error('Radii must be non-negative.');
end
if R1 == 0 || R2 == 0
    A = reshape(A, size(d_in));
    return;
end

noOverlap = d >= (R1 + R2);
contained = d <= abs(R1 - R2);
partial   = ~(noOverlap | contained);

A(contained) = pi * min(R1, R2)^2;

dp = d(partial);
dp(dp == 0) = eps;

arg1 = (dp.^2 + R1^2 - R2^2) ./ (2*dp*R1);
arg2 = (dp.^2 + R2^2 - R1^2) ./ (2*dp*R2);
arg1 = max(-1, min(1, arg1));
arg2 = max(-1, min(1, arg2));

term1 = R1^2 .* acos(arg1);
term2 = R2^2 .* acos(arg2);

rootTerm = (-dp + R1 + R2) .* (dp + R1 - R2) .* (dp - R1 + R2) .* (dp + R1 + R2);
rootTerm = max(0, rootTerm);

term3 = 0.5 .* sqrt(rootTerm);

A(partial) = term1 + term2 - term3;
A = reshape(A, size(d_in));
end
