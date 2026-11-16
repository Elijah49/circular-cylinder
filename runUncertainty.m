function u = runUncertainty(inp)
%RUNUNCERTAINTY  Uncertainty of the quantities that sit OUTSIDE Cp.
%
%   U = RUNUNCERTAINTY(INP), INP is a struct with these fields:
%     .U        uncorrected free-stream speed [m/s]
%     .Re       uncorrected Reynolds number
%     .ReCorr   blockage-corrected Reynolds number
%     .rho      air density actually used [kg/m3]
%     .mu       air viscosity actually used [Pa s]
%     .T_K      air temperature [K]
%     .p        barometric pressure [Pa]
%     .S        Sutherland constant [K]
%     .Bmm      Pitot column difference [mm]  (= |h_t - h_s|)
%     .uStat    uncertainty of the static Pitot reading [mm]
%     .uTot     uncertainty of the total Pitot reading [mm]
%     .uTairK   uncertainty of the air temperature [K]
%     .uPbaro   uncertainty of the barometric pressure [Pa]
%     .r        blockage ratio d/b
%     .CDp      uncorrected drag coefficient
%     .uCD      its standard uncertainty
%     .thDeg, .Cp   tapping angles [deg] and Cp values, same length
%     .uTheta   std uncertainty of ONE tapping's angular position [deg]
%     .dThetaAlign  coherent rotation of the whole ring implied by the
%               -30/+30 mismatch [deg]; NaN if it could not be estimated
%
%   OUTPUT struct
%     uQRel, uRhoRel, uMuRel     relative uncertainty of q, rho, mu
%     uURel                      relative uncertainty of U
%     uReRel, uRe, uReCorrRel, uReCorr
%     uCDtheta                   u(CD) that comes from the tap POSITIONS
%     uThetaSep                  angular term to add to u(thetaSep)
%
%   HOW THE CHAIN IS BUILT
%     q  = rho_e g |dh|                      -> u(q)/q = u(dh)/dh
%     U  = sqrt(2q/rho)                      -> dlnU  = 0.5 dlnq - 0.5 dln rho
%     Re = rho U d / mu                      -> dlnRe = 0.5 dln rho + 0.5 dln q
%                                                        - dln mu
%     rho appears with +0.5 in Re, not -0.5, because it enters twice: once
%     as the density in rho U d and once through U = sqrt(2q/rho). The three
%     inputs are independent, so they combine in quadrature.
%       rho = p/(R T)                -> dln rho = dln p - dln T
%       mu  = mu0 (T/T0)^1.5 (T0+S)/(T+S)
%                                    -> dln mu/dT = 1.5/T - 1/(T+S)
%     The blockage factor V/V' = 1 + 0.25 CD'(d/b) + 0.82 (d/b)^2 depends
%     on CD', so the corrected Re carries one extra term,
%       u(V/V') = 0.25 (d/b) u(CD')     applied relative to V/V' itself.
%
%   THE TAP ANGLES
%     Each tapping sits at theta_i + e_i with independent e_i of std uTheta.
%     The drag integral is evaluated at the nominal angles, so to first order
%       u^2(CD) = uTheta^2 * sum_i (dCD/dtheta_i)^2
%     evaluated here by a central difference on each angle separately. The
%     difference is taken on the whole integral because moving one tap also
%     changes the two trapezoidal weights around it, which an analytical
%     dC_p/dtheta term would miss.
%     The separation angle lives on the same angular axis, so uTheta adds
%     directly: u(thetaSep) = sqrt(u(thetaSep)_slope^2 + uTheta^2).
%
%   THE -30/+30 PAIR IS A DIFFERENT ANIMAL
%     Holes 1 and 7 are mirror images about hole 2, so a perfectly aligned
%     cylinder must give Cp(-30) = Cp(+30). A mismatch means the flow datum
%     is rotated by dThetaAlign, which shifts EVERY tap by the same angle.
%     That is a coherent bias, not a random one:
%       - it does not shrink like 1/sqrt(N) when averaged over 27 taps;
%       - its effect on CD is the outright difference CD(theta+dTheta) -
%         CD(theta), taken here as a bias of that magnitude;
%       - it moves thetaSep by dTheta itself, so the two enter linearly in
%         the reported uncertainty rather than being averaged away.
%     (dThetaAlign is estimated in the main file from dCp30 and the local
%     slope g = dCp/dtheta at +30 deg, as dTheta = -dCp30/(2g).)

uQRel   = hypot(inp.uStat, inp.uTot)/inp.Bmm;
uRhoRel = hypot(inp.uPbaro/inp.p, inp.uTairK/inp.T_K);
uMuRel  = abs(1.5/inp.T_K - 1/(inp.T_K + inp.S))*inp.uTairK;

uURel  = hypot(0.5*uQRel, 0.5*uRhoRel);
uReRel = hypot(hypot(0.5*uQRel, 0.5*uRhoRel), uMuRel);
uRe    = inp.Re*uReRel;

VV      = 1 + 0.25*inp.CDp*inp.r + 0.82*inp.r^2;
uVVRel  = 0.25*inp.r*inp.uCD/VV;
uReCorrRel = hypot(uReRel, uVVRel);
uReCorr    = inp.ReCorr*uReCorrRel;

% ---- drag uncertainty that comes from the tap positions -----------------
uCDtheta = 0;
if isfinite(inp.uTheta) && inp.uTheta > 0 && numel(inp.thDeg) >= 3
    th = inp.thDeg(:)';
    cp = inp.Cp(:)';
    s  = 0;
    for i = 1:numel(th)
        e = zeros(size(th));
        e(i) = inp.uTheta;                       % one tap moved, others fixed
        s = s + ((cdOf(th + e, cp) - cdOf(th - e, cp))/2)^2;
    end
    uCDtheta = sqrt(s);
end

uCDtap = uCDtheta;                    % random part, kept for reporting
dCDcoh = 0;
dTheta = inp.dThetaAlign;
if isfinite(dTheta) && dTheta ~= 0
    th = inp.thDeg(:)';
    cp = inp.Cp(:)';
    dCDcoh = abs(cdOf(th + dTheta, cp) - cdOf(th, cp));   % bias, not averaged
end
uCDtheta   = hypot(uCDtap, dCDcoh);
uThetaSep  = hypot(inp.uTheta, dTheta);                   % bias + random

u = struct('uQRel',uQRel,'uRhoRel',uRhoRel,'uMuRel',uMuRel, ...
           'uURel',uURel,'uReRel',uReRel,'uRe',uRe, ...
           'uReCorrRel',uReCorrRel,'uReCorr',uReCorr, ...
           'uCDtheta',uCDtheta,'uThetaSep',uThetaSep, ...
           'uCDtap',uCDtap,'dCDcoh',dCDcoh,'dThetaAlign',dTheta);
end

% -----------------------------------------------------------------------
function CD = cdOf(thDeg, Cp)
% Same integral the main script uses: CD = int_0^pi Cp cos(theta) dtheta.
CD = trapz(thDeg(:)*pi/180, Cp(:).*cos(thDeg(:)*pi/180));
end
