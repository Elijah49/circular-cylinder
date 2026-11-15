function unc = cpUncertainty(Bmm, uTap, uStat, uTot, Cp, thDeg, thetaSep, nFine, K)
%CPUNCERTAINTY  GUM uncertainty of Cp, of CD and of the separation angle.
%
%   UNC = cpUncertainty(BMM, UTAP, USTAT, UTOT, CP, THDEG, THETASEP, NFINE, K)
%
%   Standalone module: it knows nothing about the run structure, the config
%   struct or the raw readings. It takes the already-computed pressure
%   coefficients plus the uncertainty of each reading, and returns a small
%   summary. The covariance algebra, the interpolation grid and the slope
%   threshold never leave this file, so the calling script only consumes UNC.
%
%   INPUTS
%     BMM       dynamic head as a column height difference [mm]
%               (= |h_total - h_static| of the Pitot pair)
%     UTAP      std uncertainty of the tapping readings [mm]; either a
%               scalar (same for every tapping) or a 1xN vector
%     USTAT     std uncertainty of the STATIC Pitot reading [mm]
%     UTOT      std uncertainty of the TOTAL Pitot reading [mm]
%     CP        1xN pressure coefficients of the tappings, N >= 3
%     THDEG     1xN matching tap angles in degrees, increasing
%     THETASEP  separation angle from the separation criterion [deg];
%               pass NaN if it is unknown or was set by hand
%     NFINE     grid step for the slope/curvature estimate [deg]
%     K         coverage factor for the expanded uncertainty (usually 2)
%
%   OUTPUT  struct with
%     uCp       1xN standard uncertainty of every Cp            [-]
%     U95Cp     1xN expanded uncertainty, U = K*uCp             [-]
%     uCpPeak   standard uncertainty at the suction peak        [-]
%     uCD       combined standard uncertainty of CD'            [-]
%     uCDrel    uCD as a fraction of this run's CD'             [-]
%     uThetaSep standard uncertainty of the separation angle    [deg]
%     K         coverage factor used
%     uTap      1xN tapping uncertainty as used                 [mm]
%     uStat, uTot  Pitot uncertainties as used                  [mm]
%     nFine     grid step used                                  [deg]
%
%   WHY THESE EXPRESSIONS
%     Cp_i = (h_i - h_s)/(h_t - h_s) = A_i/B, so
%       dCp_i/dh_i =  1/B      dCp_i/dh_s = (Cp_i - 1)/B
%       dCp_i/dh_t = -Cp_i/B
%     With mutually independent readings:
%       Var(Cp_i)      = (1/B)^2 [ uTap_i^2 + (Cp_i uTot)^2 + ((Cp_i-1) uStat)^2 ]
%       Cov(Cp_i,Cp_j) = (1/B)^2 [ uTot^2 Cp_i Cp_j
%                                + uStat^2 (Cp_i-1)(Cp_j-1) ]      i ~= j
%     The uTap term is absent off the diagonal (each tapping is its own
%     reading), while h_s and h_t sit in EVERY Cp and are therefore
%     COHERENT: a Pitot error is not averaged down by adding more holes.
%     A single scalar u reproduces the familiar
%       u(Cp) = (u/B) sqrt(1 + Cp^2 + (Cp-1)^2).
%
%     CD is the quadratic form g'*Cov*g, g_i = w_i*cos(theta_i), w the
%     trapezoidal weights in radians (the same ones TRAPZ uses). Off-diagonal
%     sums are written as (sum)^2 - sum(squares), which collapses to
%       u(CD) = sqrt( A + B1 + B2 )
%       A  = (1/B)^2 sum( g.^2 .* (uTap.^2 + (Cp*uTot).^2 + ((Cp-1)*uStat).^2) )
%       B1 = (uTot/B)^2  [ (sum(g.*Cp))^2        - sum((g.*Cp).^2) ]
%       B2 = (uStat/B)^2 [ (sum(g.*(Cp-1)))^2    - sum((g.*(Cp-1)).^2) ]
%     A shrinks with more holes; B1 and B2 do not.
%
%     The separation angle inherits the SLOPE noise; a coherent scaling
%     leaves the shape unchanged and is ignored:
%       u(slope) = sqrt(2)*(median(uTap)/B)/spacing     [Cp per degree]
%       u(theta) = u(slope)/|d2Cp/dtheta^2|             at thetaSep
%     Floored at half the tap spacing (the data cannot resolve the knee
%     better) and capped at 45 deg (no knee was located at all).
%
%   See also CYLINDER_LAB, DRAGINTEGRAL, SEPDETECT.

cp = Cp(:); th = thDeg(:);

if numel(cp) ~= numel(th)
    error('cpUncertainty:sizeMismatch', ...
          'Cp and thDeg must have the same number of elements.');
end
if ~(Bmm > 0)
    error('cpUncertainty:badB', 'Bmm must be positive.');
end

n = numel(cp);
uTap = uTap(:);
if isscalar(uTap), uTap = repmat(uTap, n, 1); end
if numel(uTap) ~= n
    error('cpUncertainty:sizeMismatch', ...
          'uTap must be a scalar or have one element per tapping.');
end

% ---- per-tapping standard uncertainty
uCp = sqrt( uTap.^2 + (cp*uTot).^2 + ((cp-1)*uStat).^2 ) / Bmm;

% ---- trapezoidal weights in radians, exactly as TRAPZ builds them
w = zeros(n,1);
for i = 1:n
    if i == 1
        w(i) = (th(2) - th(1))/2;
    elseif i == n
        w(i) = (th(n) - th(n-1))/2;
    else
        w(i) = (th(i+1) - th(i-1))/2;
    end
end
g = (w*pi/180).*cosd(th);

% ---- drag coefficient: per-tapping term + the two coherent terms
A  = sum( g.^2 .* (uTap.^2 + (cp*uTot).^2 + ((cp-1)*uStat).^2) ) / Bmm^2;
s1 = sum(g.*cp);        q1 = sum((g.*cp).^2);
s2 = sum(g.*(cp-1));    q2 = sum((g.*(cp-1)).^2);
B1 = (uTot/Bmm)^2  * (s1^2 - q1);
B2 = (uStat/Bmm)^2 * (s2^2 - q2);
uCD = sqrt(max(A + B1 + B2, 0));

CDd = sum(g.*cp);                              % CD' from the same weights
if abs(CDd) > 1e-9
    uRel = uCD/abs(CDd);
else
    uRel = NaN;
end

% ---- separation angle
thF = (0:nFine:180)';
cpF = interp1(th, cp, thF, 'pchip');
d1  = gradient(cpF, nFine);
d2  = gradient(d1,  nFine);
dTh = median(diff(th));                        % tap spacing [deg]
uSlope = sqrt(2)*(median(uTap)/Bmm)/dTh;       % first-difference slope noise

if isfinite(thetaSep)
    [~, iS]   = min(abs(thF - thetaSep));
    curv      = max(abs(d2(iS)), 1e-6);
    uThetaSep = uSlope/curv;
    uThetaSep = min(max(uThetaSep, 0.5*dTh), 45);   % floor / cap
else
    uThetaSep = NaN;
end

% ---- summary (no matrices, no grids, nothing else escapes)
[~, iPk]      = min(cp);
unc.uCp       = uCp;
unc.U95Cp     = K*uCp;
unc.uCpPeak   = uCp(iPk);
unc.uCD       = uCD;
unc.uCDrel    = uRel;
unc.uThetaSep = uThetaSep;
unc.K         = K;
unc.uTap      = uTap;
unc.uStat     = uStat;
unc.uTot      = uTot;
unc.nFine     = nFine;
end
