%% =========================================================================
%  cylinder_lab.m
%  Pressure distribution on a circular cylinder -> Cp, CL, CD, separation
%  point and blockage correction (Allen & Vincenti 1944, as used by
%  Roshko 1961, JFM 10, 345-356, eqs 1-3)
%
%  Engineering Practice 2 - Aerodynamics 2, circular cylinder experiment
%  Open-circuit wind tunnel RODH 115 (ex L50), d = 102 mm, b = 460 mm
%
%  HOW TO USE
%    1. Dry run: leave cfg.demo = true and press Run. Synthetic readings are
%       generated so the whole pipeline, the plots and the table can be
%       checked before the lab.
%    2. In the lab: set cfg.demo = false. The 28 readings, the three Pitot /
%       room-condition repeats and the board scale division are typed into
%       readings.csv (one row per hole, repeats side by side); the script
%       picks that file up automatically and then reports BOTH a conservative
%       uncertainty (cfg.uRead) and the one actually measured from the
%       repeats. Without the file the manual block in section 2 is used and
%       only the conservative figure is available; a run left as [] is
%       skipped instead of crashing the script.
%    3. Outputs go to ./results : CSV table + PNG figures + .fig files.
%
%  CONVENTIONS
%    theta = 0 deg at the forward stagnation point (hole 2), positive
%    towards hole 28 (theta = 180 deg, rear of the cylinder); Table 1 of
%    the handout gives the 28 tapping positions from -30 to 180 deg.
%    Cp is referenced to the free-stream static pressure of the Pitot tube,
%    so only DIFFERENCES of manometer readings enter and a common datum
%    offset of the multi-tube manometer cancels.
%
%    Reading convention: a top-fed board (tubing to the top of each tube,
%    lower end in a common open ethanol reservoir) gives a column that
%    FALLS as the applied pressure rises; a bottom-fed board gives the
%    opposite. The script does not care. Because the liquid density, g and
%    the sign of the board are all common to every tube, Cp reduces to a
%    ratio of two column differences, Cp = (h_i - h_s)/(h_t - h_s), and the
%    sign cancels. The only place the sign could bite is "which Pitot tube
%    is total", and that is settled physically, not by magnitude: at hole 2
%    (theta = 0) the surface pressure IS the total pressure, so the total
%    tube is the one whose reading matches h(hole 2). The script prints this
%    comparison for every run - check it, it should be within a millimetre
%    or two.
%
%  Reference equations (Roshko 1961, eqs 1-3), V' = measured tunnel speed,
%  V = corrected speed, h = tunnel breadth (= the handout's b):
%    V /V'    = 1 + (1/4) CD'  (d/h) + 0.82 (d/h)^2
%    CD /CD'  = 1 - (1/2) CD'  (d/h) - 2.5  (d/h)^2
%    Cp - 1   = (V'/V)^2 (Cp' - 1)
%% =========================================================================

clear; clc; close all;

%% ------------------------------ 1. CONFIG -------------------------------
cfg.demo        = true;        % <<< SET TO false IN THE LAB
cfg.dia         = 0.102;       % cylinder diameter [m]
cfg.tunnelB     = 0.46;        % tunnel test-section breadth [m]  (= handout b)
cfg.rhoEth      = 789;         % ethanol density [kg/m^3] (handout)
cfg.g           = 9.80665;     % gravity [m/s^2]
cfg.Rair        = 287.05;      % gas constant for air [J/kg/K]
cfg.mu0         = 1.716e-5;    % Sutherland reference viscosity [Pa s]
cfg.T0          = 273.15;      % [K]
cfg.Suth        = 110.4;       % Sutherland constant [K]
cfg.units       = 'mm';        % manometer reading unit ('mm' or 'm')
cfg.pitotOrder  = 'auto';      % 'auto' | 'tube35static' | 'tube36static'
cfg.nFine       = 0.25;        % refined integration grid [deg]
cfg.thetaSepManual = [];       % per-run override, e.g. [80 88 118 124 NaN]
cfg.uRead       = 0.5;         % CONSERVATIVE std uncertainty of one reading [mm]
                               % merged type A + type B; used when no repeat
                               % readings are available (see cfg.dataFile)
cfg.uScale      = 1.0;         % board scale resolution [mm] -> type B term
                               % u_B = (uScale/2)/sqrt(3) = 0.289 mm
cfg.dataFile    = 'readings.csv';   % repeat-reading log. If this file exists
                               % and cfg.demo is false, both the readings AND
                               % the measured type A uncertainty come from it
cfg.sepUncK     = 2;           % coverage factor for the reported U95
cfg.outDir      = fullfile(pwd,'results');
cfg.TairDefault  = 21.0;       % <<< update to the tunnel room temperature
cfg.pBaroDefault = 101300;     % <<< update to the day's barometric pressure

% Reference CD(Re) curve: APPROXIMATE hand-digitisation of Anderson Fig 3.44
% and Roshko Fig 2. Used only for the "equivalent Reynolds number" step.
% Replace with your own readings of the figures if your copy differs.
cfg.refRe = [1e4   1e5   2e5   2.5e5 3e5   5e5   7e5   1e6   2e6   3.5e6 1e7];
cfg.refCd = [1.20  1.20  1.05  0.60  0.40  0.35  0.35  0.38  0.50  0.70  0.70];

% Tapping positions, Table 1 of the handout (deg)
thetaHole = [-30  0   5  10  15  20  30  40  50  55  60  65  70  75  80  85 ...
              90  95 100 105 110 115 120 130 140 150 165 180];

%% ------------------- 2. RAW DATA  (manometer readings) ------------------
% Three sources, in this order:
%   cfg.demo = true          -> synthetic runs, nothing to fill in
%   readings.csv present     -> readings AND repeat-reading uncertainty
%                               (cfg.demo must be false)
%   otherwise                -> the manual block below, uncertainty then
%                               falls back to the conservative cfg.uRead
uMeas = [];                    % one entry per run, aligned with runs
if cfg.demo
    runs = demoRuns(thetaHole, cfg);
elseif exist(cfg.dataFile,'file') == 2
    [runs, uMeas] = loadReadings(cfg.dataFile, cfg);
    fprintf('[data] readings and repeat-reading uncertainty taken from %s\n', ...
            cfg.dataFile);
else
    runs = struct('label',{},'tripped',{},'rpm',{},'Tair',{},'pBaro',{}, ...
                  'h',{},'hPitot',{});

    runs(1).label   = 'Smooth - low speed';
    runs(1).tripped = false;
    runs(1).rpm     = 700;      % ~15-16 m/s
    runs(1).Tair    = [ ];      % <<< air temperature [degC]
    runs(1).pBaro   = [ ];      % <<< barometric pressure [Pa], e.g. 101300
    runs(1).hPitot  = [ ];      % <<< [ tube35 , tube36 ]
    runs(1).h       = [ ];      % <<< 28 readings, holes 1..28

    runs(2).label   = 'Smooth - high speed';
    runs(2).tripped = false;
    runs(2).rpm     = 1350;     % ~30 m/s
    runs(2).Tair    = [ ];
    runs(2).pBaro   = [ ];
    runs(2).hPitot  = [ ];
    runs(2).h       = [ ];

    runs(3).label   = 'Trip - low speed';
    runs(3).tripped = true;
    runs(3).rpm     = 700;
    runs(3).Tair    = [ ];
    runs(3).pBaro   = [ ];
    runs(3).hPitot  = [ ];
    runs(3).h       = [ ];

    runs(4).label   = 'Trip - high speed';
    runs(4).tripped = true;
    runs(4).rpm     = 1350;
    runs(4).Tair    = [ ];
    runs(4).pBaro   = [ ];
    runs(4).hPitot  = [ ];
    runs(4).h       = [ ];

    runs(5).label   = 'Trip - extra speed (choose with demonstrator)';
    runs(5).tripped = true;
    runs(5).rpm     = [ ];
    runs(5).Tair    = [ ];
    runs(5).pBaro   = [ ];
    runs(5).hPitot  = [ ];
    runs(5).h       = [ ];
end

%% --------------------------- 3. PROCESS --------------------------------
if ~exist(cfg.outDir,'dir'), mkdir(cfg.outDir); end
pre = ''; if cfg.demo, pre = 'DEMO_'; end

r = cfg.dia/cfg.tunnelB;                       % blockage ratio d/h
fprintf('\n================ CYLINDER LAB ================\n');
fprintf('d = %.1f mm, tunnel breadth b = %.1f mm, d/b = %.4f\n', ...
        cfg.dia*1e3, cfg.tunnelB*1e3, r);
if cfg.demo
    fprintf(['[NOTE] cfg.demo = true: readings below are SYNTHETIC, for\n' ...
             '       testing the code only. Do NOT report them.\n']);
end

thTheory = (0:0.5:180)';
cpTheory = 1 - 4*sind(thTheory).^2;

fields = {'label','tripped','rpm','Tair','pBaro','rho','mu','U','Ucorr', ...
          'Re','Recorr','Cp','CpCorr','CDprime','CDfine','CLsym','CL', ...
          'dCp30','CDcorrInt','CDcorrEq2','VV','thetaSep','thetaCpmin', ...
          'Cpmin','hStat','hTot','iStat','ReEq','unc','uncA'};
res = struct();
for fi = 1:numel(fields), res.(fields{fi}) = []; end
res = repmat(res, 1, 0);               % 0x0 struct array, all fields present

for k = 1:numel(runs)
    run = runs(k);

    if isempty(run.h) || numel(run.h) ~= 28
        fprintf('\n[skip] %-40s : no complete set of 28 readings\n', run.label);
        continue
    end
    if isempty(run.hPitot) || numel(run.hPitot) ~= 2 || any(isnan(run.hPitot))
        fprintf('\n[skip] %-40s : Pitot readings (tubes 35/36) missing\n', run.label);
        continue
    end

    Tair  = pick(run.Tair , cfg.TairDefault);
    pBaro = pick(run.pBaro, cfg.pBaroDefault);
    [rhoAir, muAir] = airProps(Tair, pBaro, cfg);

    % ---- Pitot-static tube -> tunnel speed
    hStag = run.h(2);                      % hole 2 sits at theta = 0
    [hStat, hTot, iStat, iTot] = assignPitot(run.hPitot, hStag, cfg);
    dq = toMetres(hTot - hStat, cfg);                % signed, metres of ethanol
    q  = cfg.rhoEth*cfg.g*abs(dq);                   % dynamic head [Pa], > 0
    if ~isfinite(q) || q <= 0
        fprintf('\n[skip] %-40s : the two Pitot tubes read the same\n', run.label);
        continue
    end
    U  = sqrt(2*q/rhoAir);
    Re = rhoAir*U*cfg.dia/muAir;

    % ---- pressure coefficient
    % Cp is a ratio of two column differences taken on the same bank, so
    % rho_e, g and the sign convention of the board all cancel: whether the
    % column rises or falls with pressure, Cp comes out identical.
    Cp = toMetres(run.h - hStat, cfg) / dq;

    % ---- stagnation cross-check: the pressure at hole 2 IS the total pressure
    fprintf('  stagnation check: h(hole 2) = %.1f mm, h(total tube) = %.1f mm, diff = %+.1f mm\n', ...
            hStag, hTot, hStag - hTot);

    % ---- drag and lift on the measured data
    CDp    = dragIntegral(thetaHole(2:end), Cp(2:end), cfg);
    CDfine = dragIntegralFine(thetaHole(2:end), Cp(2:end), cfg);
    [CLsym, CLasym] = liftIntegral(thetaHole, Cp, cfg);

    % ---- blockage correction (Roshko 1961, eqs 1-3)
    VV        = 1 + 0.25*CDp*r + 0.82*r^2;                       % V/V'
    CDcorrEq2 = (1 - 0.5*CDp*r - 2.5*r^2)*CDp;                   % eq 2
    CpCorr    = 1 + (1/VV)^2*(Cp - 1);                           % eq 3
    CDcorrInt = dragIntegral(thetaHole(2:end), CpCorr(2:end), cfg);
    Ucorr     = U*VV;
    Recorr    = rhoAir*Ucorr*cfg.dia/muAir;

    % ---- separation point
    if ~isempty(cfg.thetaSepManual) && numel(cfg.thetaSepManual) >= k ...
            && ~isnan(cfg.thetaSepManual(k))
        thetaSep   = cfg.thetaSepManual(k);
        thetaCpmin = NaN;
    else
        [thetaSep, thetaCpmin] = sepDetect(thetaHole(2:end), Cp(2:end), cfg);
    end

    % ---- uncertainty of Cp, CD and the separation angle (GUM)
    % The maths lives in cpUncertainty.m; only its summary comes back. The
    % 1..28 angles are passed WITHOUT the extra -30 deg point so that the
    % trapezoidal weights match dragIntegral exactly.
    %   unc  - conservative: every reading taken as +-cfg.uRead
    %   uncA - measured: type A from the repeat readings in readings.csv,
    %          combined with the type B scale term; NaN if not available
    if exist('cpUncertainty','file') == 2
        Bmm = abs(hTot - hStat);
        unc = cpUncertainty(Bmm, cfg.uRead, cfg.uRead, cfg.uRead, ...
                            Cp(2:end), thetaHole(2:end), thetaSep, ...
                            cfg.nFine, cfg.sepUncK);
        if ~isempty(uMeas)
            % which Pitot tube is static was decided above (iStat)
            if iStat == 35
                uSt = uMeas(k).u35;  uTt = uMeas(k).u36;
            else
                uSt = uMeas(k).u36;  uTt = uMeas(k).u35;
            end
            % uTap is indexed by hole 1..28, Cp by hole 2..28 -> drop the first
            uncA = cpUncertainty(Bmm, uMeas(k).uTap(2:end), uSt, uTt, ...
                                 Cp(2:end), thetaHole(2:end), ...
                                 thetaSep, cfg.nFine, cfg.sepUncK);
        else
            uncA = nanUnc(cfg);
        end
    else
        % Module absent (e.g. only the main file was copied across). Keep
        % running and report NaN rather than pretending it was quantified.
        fprintf('  [note] cpUncertainty.m not found - uncertainty skipped\n');
        unc  = nanUnc(cfg);
        uncA = nanUnc(cfg);
    end

    % ---- equivalent Reynolds number of the tripped case
    ReEq = equivalentRe(CDcorrEq2, cfg);

    % ---- symmetry diagnostic across the -30 / +30 pair
    dCp30 = Cp(thetaHole == -30) - Cp(thetaHole == 30);

    res(end+1) = struct('label',run.label,'tripped',run.tripped,'rpm',run.rpm, ...
        'Tair',Tair,'pBaro',pBaro,'rho',rhoAir,'mu',muAir,'U',U,'Ucorr',Ucorr, ...
        'Re',Re,'Recorr',Recorr,'Cp',Cp,'CpCorr',CpCorr,'CDprime',CDp, ...
        'CDfine',CDfine,'CLsym',CLsym,'CL',CLasym,'dCp30',dCp30, ...
        'CDcorrInt',CDcorrInt,'CDcorrEq2',CDcorrEq2,'VV',VV, ...
        'thetaSep',thetaSep,'thetaCpmin',thetaCpmin,'Cpmin',min(Cp), ...
        'hStat',hStat,'hTot',hTot,'iStat',iStat,'ReEq',ReEq,'unc',unc, ...
        'uncA',uncA);

    fprintf('\n--- %s  (tripped = %d, %.0f rpm) ---\n', run.label, run.tripped, run.rpm);
    fprintf('  rho  = %.3f kg/m3 , mu = %.3e Pa s , T = %.1f C , p = %.0f Pa\n', ...
            rhoAir, muAir, Tair, pBaro);
    fprintf('  Pitot: tube %d = static (%.1f mm), tube %d = total (%.1f mm), q = %.2f Pa\n', ...
            iStat, hStat, iTot, hTot, q);
    fprintf('  U'' = %.2f m/s  ->  U = %.2f m/s  ,  Re'' = %.3e  ->  Re = %.3e\n', ...
            U, Ucorr, Re, Recorr);
    fprintf('  CD'' = %.4f (trapz on tappings) | %.4f (refined)    CL'' = %.4f\n', ...
            CDp, CDfine, CLasym);
    fprintf('  CD eq(2) = %.4f   CD from integrated corrected Cp = %.4f   (diff %.1f %%)\n', ...
            CDcorrEq2, CDcorrInt, 100*(CDcorrInt - CDcorrEq2)/CDcorrEq2);
    fprintf('  separation = %.1f deg   (suction peak at %.1f deg, Cp_min = %.2f)\n', ...
            thetaSep, thetaCpmin, min(Cp));
    if isfinite(unc.uCD)
        fprintf('  u (cons., u=%.1f mm): u(Cp_peak)=%.3f  u(CD'')=%.3f (%.1f %%)  u(sep)=%.1f deg\n', ...
                cfg.uRead, unc.uCpPeak, unc.uCD, 100*unc.uCDrel, unc.uThetaSep);
    end
    if isfinite(uncA.uCD)
        fprintf('  u (measured, n>=2):   u(Cp_peak)=%.3f  u(CD'')=%.3f (%.1f %%)  u(sep)=%.1f deg\n', ...
                uncA.uCpPeak, uncA.uCD, 100*uncA.uCDrel, uncA.uThetaSep);
    end
    fprintf('  +-30 deg mismatch: dCp = %+.3f  (CL by symmetry = %.4f, from data = %+.4f)\n', ...
            dCp30, CLsym, CLasym);
end

if isempty(res)
    fprintf('\nNothing processed. Set cfg.demo = false and fill at least one run.\n');
    return
end

%% --------------------------- 4. SUMMARY TABLE --------------------------
printTable(res);

csvFile = fullfile(cfg.outDir,[pre 'cylinder_results.csv']);
writeCsv(res, csvFile);
fprintf('Table written to: %s\n', csvFile);

%% ------------------------------ 5. PLOTS -------------------------------
fig1 = plotCp(res, thetaHole, thTheory, cpTheory, false);
fig2 = plotCp(res, thetaHole, thTheory, cpTheory, true);
fig3 = plotCdRe(res, cfg);
fig4 = plotPolar(res, thetaHole);

% .fig first: exportgraphics leaves the handles in a state savefig rejects
savefig(fig1, fullfile(cfg.outDir,[pre 'Cp_vs_theta.fig']));
savefig(fig2, fullfile(cfg.outDir,[pre 'Cp_corrected.fig']));
for fk = [fig1 fig2 fig3 fig4]
    for a = findall(fk,'Type','axes')'
        if isprop(a,'Toolbar') && ~isempty(a.Toolbar)
            a.Toolbar.Visible = 'off';   % keeps the export warning quiet
        end
    end
end
exportgraphics(fig1, fullfile(cfg.outDir,[pre 'Cp_vs_theta.png']), 'Resolution', 150);
exportgraphics(fig2, fullfile(cfg.outDir,[pre 'Cp_corrected.png']), 'Resolution', 150);
exportgraphics(fig3, fullfile(cfg.outDir,[pre 'CD_vs_Re.png']), 'Resolution', 150);
exportgraphics(fig4, fullfile(cfg.outDir,[pre 'pressure_polar.png']), 'Resolution', 150);
fprintf('Figures written to: %s\n', cfg.outDir);

fprintf('\n================ END ================\n');


%% ====================== LOCAL FUNCTIONS ================================
function v = pick(v, default)
if isempty(v) || any(isnan(v)), v = default; end
end

% -----------------------------------------------------------------------
function h = toMetres(h, cfg)
if strcmpi(cfg.units,'mm'), h = h*1e-3; end
end

% -----------------------------------------------------------------------
function [rho, mu] = airProps(TdegC, pPa, cfg)
T = TdegC + 273.15;
rho = pPa/(cfg.Rair*T);                                            % ideal gas
mu  = cfg.mu0*(T/cfg.T0)^1.5*(cfg.T0 + cfg.Suth)/(T + cfg.Suth);   % Sutherland
end

% -----------------------------------------------------------------------
function [hStat, hTot, iStat, iTot] = assignPitot(hp, hStag, cfg)
% Which of the two Pitot tubes is the total port?
% Do NOT decide this by "the larger reading". On a top-fed board (flexible
% tubing to the TOP of each tube, lower end dipping into a common ethanol
% reservoir) the column FALLS as the applied pressure rises, so the larger
% reading belongs to the STATIC port; on a bottom-fed board it is the
% opposite. Use physics instead: hole 2 sits at theta = 0, the stagnation
% point, where the surface pressure IS the total pressure. So h(hole 2) must
% equal the reading of the total tube, in either convention. The resulting
% speed is then still checked against 700 rpm -> 15-16 m/s and
% 1300-1400 rpm -> ~30 m/s; if it is grossly off, set cfg.pitotOrder by hand.
switch lower(cfg.pitotOrder)
    case 'tube35static', iStat = 1; iTot = 2;
    case 'tube36static', iStat = 2; iTot = 1;
    otherwise
        if isfinite(hStag)
            d = abs(hp - hStag);
            [~, iTot] = min(d);  iStat = 3 - iTot;
            if abs(d(1) - d(2)) < 0.5        % both within noise: nothing to gain
                [~, iTot] = max(hp);  iStat = 3 - iTot;
            end
        else
            [~, iTot] = max(hp);  iStat = 3 - iTot;
        end
end
hStat = hp(iStat);  hTot = hp(iTot);
end

% -----------------------------------------------------------------------
function CD = dragIntegral(thDeg, Cp, cfg)   %#ok<INUSD>
th = thDeg(:)*pi/180;  cp = Cp(:);
CD = trapz(th, cp.*cos(th));                  % CD = int_0^pi Cp cos(theta) dtheta
end

% -----------------------------------------------------------------------
function CDfine = dragIntegralFine(thDeg, Cp, cfg)
thF = (0:cfg.nFine:180)';
cpF = interp1(thDeg(:), Cp(:), thF, 'pchip');
CDfine = trapz(thF*pi/180, cpF.*cos(thF*pi/180));
end

% -----------------------------------------------------------------------
function [CLsym, CLasym] = liftIntegral(thDeg, Cp, cfg)
% CL = (1/2) int_0^{2pi} Cp sin(theta) dtheta.
% Only -30 .. 180 deg is instrumented, so the missing arc is mirrored.
% CLsym (pure mirror) is identically zero; CLasym keeps the single measured
% point at -30 deg, so a non-zero value means the flow was not symmetric
% (cylinder misalignment, or a genuine lift).
tf = (-180:cfg.nFine:180)';
pos = thDeg >= 0;
[thP, ip] = sort(thDeg(pos));
cpP = Cp(pos);  cpP = cpP(ip);
cpMirror = interp1(thP, cpP, abs(tf), 'pchip');
CLsym = 0.5*trapz(tf*pi/180, cpMirror.*sin(tf*pi/180));

[~, iM30] = min(abs(thDeg + 30));
cpM30 = Cp(iM30);
cpAsym = cpMirror;
in = tf >= -30 & tf < 0;
w = (tf(in) + 30)/30;                     % 0 at -30 deg, 1 at 0 deg
cpAsym(in) = (1-w).*cpM30 + w.*cpMirror(in);
CLasym = 0.5*trapz(tf*pi/180, cpAsym.*sin(tf*pi/180));
end

% -----------------------------------------------------------------------
function [thetaSep, thetaCpmin] = sepDetect(thDeg, Cp, cfg)
% Separation from the pressure distribution. Two numbers are returned:
%   thetaCpmin - angle of the suction peak, i.e. where the attached flow
%                reaches its lowest pressure and recovery starts;
%   thetaSep   - angle where the pressure recovery is finished and the
%                constant base-pressure plateau begins. This is the
%                operational definition used in the report: behind
%                separation the surface pressure no longer responds to
%                the outer flow, so dCp/dtheta falls to almost zero.
% The search deliberately starts AFTER the recovery gradient has been
% established, because dCp/dtheta is also zero at the suction peak itself.
thF = (0:cfg.nFine:180)';
cpF = interp1(thDeg(:), Cp(:), thF, 'pchip');
d   = gradient(cpF, cfg.nFine);           % dCp/dtheta, per degree
ds  = movmean(d, 9);

[~, iMin]   = min(cpF);
thetaCpmin  = thF(iMin);

post  = ds(iMin:end);
dmax  = max(post);
iRec  = find(post > 0.4*dmax, 1, 'first');
if isempty(iRec) || dmax < 1e-3
    thetaSep = thetaCpmin;                % no clear recovery: report the peak
    return
end
iRec = iMin + iRec - 1;

nWin = max(3, round(12/cfg.nFine));       % require a 12 deg flat window
thetaSep = thF(end);
for i = iRec:(numel(ds) - nWin)
    w = ds(i:i+nWin);
    if abs(ds(i)) < 0.1*dmax && mean(abs(w)) < 0.12*dmax
        thetaSep = thF(i);
        break
    end
end
end

% -----------------------------------------------------------------------
function ReEq = equivalentRe(CD, cfg)
% Equivalent Reynolds number from CD(Re). CD(Re) is double valued around
% the critical point, so only the post-minimum (rising) branch is used,
% where a tripped, transitioned boundary layer belongs.
refRe = cfg.refRe(:);  refCd = cfg.refCd(:);
[~, iMin] = min(refCd);
x = log10(refRe(iMin:end));  y = refCd(iMin:end);
[y, ia] = unique(y, 'stable');   x = x(ia);   % CD repeats above the minimum
if CD <= y(1) || CD >= y(end)
    ReEq = NaN;                            % outside the tabulated range
else
    ReEq = 10^interp1(y, x, CD, 'linear');
end
end

% -----------------------------------------------------------------------
function printTable(res)
fprintf('\n===================== RESULTS TABLE =====================\n');
fprintf('%-24s %4s %7s %7s %9s %7s %7s %7s %7s %7s\n', ...
        'Case','trip','U''','U','Re''','CD''','CD eq2','CD int','CL''','sep');
fprintf('%s\n', repmat('-',1,80));
for k = 1:numel(res)
    q = res(k);
    fprintf('%-24s %4d %7.2f %7.2f %9.2e %7.3f %7.3f %7.3f %7.3f %7.1f\n', ...
        shorten(q.label,24), q.tripped, q.U, q.Ucorr, q.Re, q.CDprime, ...
        q.CDcorrEq2, q.CDcorrInt, q.CL, q.thetaSep);
end
fprintf('%s\n', repmat('-',1,80));
fprintf('CD eq2 = Roshko equation (2), direct correction of CD''\n');
fprintf('CD int = integration of the corrected Cp distribution (eq 3)\n');
fprintf('CL''     = symmetry assumed; the printed value above uses the\n');
fprintf('         measured -30 deg point, i.e. it is the asymmetry term.\n');
for k = 1:numel(res)
    if res(k).tripped
        fprintf('  %-22s CD = %.3f  ->  equivalent Re ~ %.2e\n', ...
                shorten(res(k).label,22), res(k).CDcorrEq2, res(k).ReEq);
    end
end

% ---- uncertainty summary (random part only; systematic terms are listed
%      separately in the report and are NOT added here)
if isfinite(res(1).unc.uCD)
    fprintf('\nUncertainty, CONSERVATIVE u(h) = %.1f mm, k = %d\n', ...
            res(1).unc.uStat, res(1).unc.K);
    for k = 1:numel(res)
        q = res(k);
        fprintf('  %-24s u(Cp_peak) = %.3f   u(CD'') = %.3f (%.1f %%)   u(sep) = %.1f deg\n', ...
                shorten(q.label,24), q.unc.uCpPeak, q.unc.uCD, ...
                100*q.unc.uCDrel, q.unc.uThetaSep);
    end
end
if isfinite(res(1).uncA.uCD)
    fprintf('\nUncertainty, MEASURED from repeat readings (type A + type B), k = %d\n', ...
            res(1).uncA.K);
    for k = 1:numel(res)
        q = res(k);
        fprintf('  %-24s u(Cp_peak) = %.3f   u(CD'') = %.3f (%.1f %%)   u(sep) = %.1f deg\n', ...
                shorten(q.label,24), q.uncA.uCpPeak, q.uncA.uCD, ...
                100*q.uncA.uCDrel, q.uncA.uThetaSep);
    end
end
end

function s = shorten(s, n)
if numel(s) > n, s = [s(1:n-2) '..']; end
end

% -----------------------------------------------------------------------
function writeCsv(res, file)
hdr = {'label','tripped','rpm','Tair_C','pBaro_Pa','rho_kgm3','mu_Pas', ...
       'U_prime_ms','U_corr_ms','Re_prime','Re_corr','CD_prime','CD_refined', ...
       'CL_prime','CL_symmetry','dCp_30deg','CD_corr_eq2','CD_corr_integral', ...
       'V_over_Vprime','theta_sep_deg','theta_Cpmin_deg','Cp_min','Re_equiv', ...
       'u_Cp_peak','u_CD','u_thetaSep_deg', ...
       'u_Cp_peak_meas','u_CD_meas','u_thetaSep_meas'};
fid = fopen(file,'w');
fprintf(fid,'%s,',hdr{1:end-1});  fprintf(fid,'%s\n',hdr{end});
for k = 1:numel(res)
    q = res(k);
    fprintf(fid,'%s,%d,%.0f,%.2f,%.0f,%.4f,%.4e,%.3f,%.3f,%.4e,%.4e,', ...
        q.label, q.tripped, q.rpm, q.Tair, q.pBaro, q.rho, q.mu, ...
        q.U, q.Ucorr, q.Re, q.Recorr);
    fprintf(fid,'%.4f,%.4f,%.4f,%.4f,%.5f,%.4f,%.4f,%.4f,%.1f,%.1f,%.3f,%.3e,%.3f,%.3f,%.1f,%.3f,%.3f,%.1f\n', ...
        q.CDprime, q.CDfine, q.CL, q.CLsym, q.dCp30, q.CDcorrEq2, ...
        q.CDcorrInt, q.VV, q.thetaSep, q.thetaCpmin, q.Cpmin, q.ReEq, ...
        q.unc.uCpPeak, q.unc.uCD, q.unc.uThetaSep, ...
        q.uncA.uCpPeak, q.uncA.uCD, q.uncA.uThetaSep);
end
fclose(fid);
end

% -----------------------------------------------------------------------
function f = plotCp(res, thetaHole, thTheory, cpTheory, corrected)
% Figure 1: measured Cp of every case against the potential-flow curve.
% Figure 2: for every case BOTH the uncorrected (dashed) and the
%           blocked-corrected (solid) distribution, as the handout asks.
if corrected
    ttl = 'Blockage-corrected (solid) and uncorrected (dashed) C_p';
else
    ttl = 'C_p against \theta';
end
f = figure('Color','w','Position',[80 80 1000 560]); hold on; box on; grid on;
n = numel(res);
cols = lines(n);
ht = plot(thTheory, cpTheory, 'k--', 'LineWidth', 1.4);

if corrected
    h = gobjects(1, 2*n);
    lab = cell(1, 2*n + 1);
    for k = 1:n
        h(2*k-1) = plot(thetaHole, res(k).Cp, '--', 'Color', cols(k,:), 'LineWidth', 1.0);
        h(2*k)   = plot(thetaHole, res(k).CpCorr, '-o', 'Color', cols(k,:), ...
                        'LineWidth', 1.3, 'MarkerSize', 5, 'MarkerFaceColor','w');
        lab{2*k-1} = [res(k).label ' - uncorrected'];
        lab{2*k}   = [res(k).label ' - corrected'];
    end
else
    h = gobjects(1, n);
    lab = cell(1, n + 1);
    for k = 1:n
        h(k) = plot(thetaHole, res(k).Cp, '-o', 'Color', cols(k,:), ...
                    'LineWidth', 1.2, 'MarkerSize', 5, 'MarkerFaceColor','w');
        cpSep = interp1(thetaHole, res(k).Cp, res(k).thetaSep, 'pchip');
        plot(res(k).thetaSep, cpSep, 'x', 'Color', cols(k,:), ...
             'MarkerSize', 10, 'LineWidth', 1.6);
        lab{k} = sprintf('%s  (sep %.0f%c)', res(k).label, res(k).thetaSep, 176);
    end
end
lab{end} = 'Potential flow:  C_p = 1 - 4 sin^2\theta';

xlabel('\theta  [deg]'); ylabel('C_p');
xlim([0 180]); set(gca,'XTick',0:20:180);
title(ttl);
legend([h ht], lab, 'Location','eastoutside','FontSize',8);
hold off
end

% -----------------------------------------------------------------------
function f = plotCdRe(res, cfg)
f = figure('Color','w','Position',[120 120 940 540]); hold on; box on; grid on;
colsc = lines(numel(res));
n = numel(res);
href = semilogx(cfg.refRe, cfg.refCd, 'k--', 'LineWidth', 1.4);
h1 = gobjects(1,n);  h2 = gobjects(1,n);
hdum = plot(nan, nan, 'p', 'Color','k', 'MarkerFaceColor',[1 0.85 0], 'MarkerSize', 12);
for k = 1:n
    q = res(k);
    h1(k) = semilogx(q.Re, q.CDprime, 'o', 'Color', colsc(k,:), ...
                     'MarkerFaceColor','w', 'MarkerSize', 9);
    h2(k) = semilogx(q.Recorr, q.CDcorrEq2, 'o', 'Color', colsc(k,:), ...
                     'MarkerFaceColor', colsc(k,:), 'MarkerSize', 9);
    if q.tripped && ~isnan(q.ReEq)
        semilogx(q.ReEq, q.CDcorrEq2, 'p', 'Color','k', ...
                 'MarkerFaceColor',[1 0.85 0], 'MarkerSize', 13);
    end
end
set(gca, 'XScale','log');          % the dummy plot handle resets it to linear
xlim([5e4 2e7]);
xlabel('Re = \rho U d / \mu'); ylabel('C_D');
title('C_D against Re   (hollow = uncorrected, filled = corrected)');
lab = cell(1, 2*n + 2);
lab{1} = 'Reference curve (approx., Anderson Fig 3.44 / Roshko Fig 2)';
for k = 1:n, lab{1+k}     = [res(k).label ' - uncorrected']; end
for k = 1:n, lab{1+n+k}   = [res(k).label ' - corrected'];  end
lab{end} = 'Equivalent Re of the tripped case';
legend([href, h1, h2, hdum], lab, 'Location','eastoutside','FontSize',8);
hold off
end

% -----------------------------------------------------------------------
function f = plotPolar(res, thetaHole)
% Cross-section plane, r = -Cp, mirrored about the flow axis.
f = figure('Color','w','Position',[160 160 640 640]); hold on; box on; grid on;
cols = lines(numel(res));
thF = (0:1:180)';
h = gobjects(1, numel(res));
for k = 1:numel(res)
    cp = interp1(thetaHole, res(k).Cp, thF, 'pchip');
    x = -cp.*cosd(thF);  y = cp.*sind(thF);
    h(k) = plot(x, y, '-', 'Color', cols(k,:), 'LineWidth', 1.2);
    plot(x, -y, '-', 'Color', cols(k,:), 'LineWidth', 1.2, 'HandleVisibility','off');
end
axis equal; xlabel('-C_p cos\theta'); ylabel('C_p sin\theta');
title('Pressure distribution in the cross-section plane (r = -C_p)');
lab = cell(1,numel(res));
for k = 1:numel(res), lab{k} = res(k).label; end
legend(h, lab, 'Location','best','FontSize',8);
hold off
end

% -----------------------------------------------------------------------
function runs = demoRuns(thetaHole, cfg)
% SYNTHETIC readings so the whole pipeline can be tested before the lab.
% The Cp shapes are hand-built anchor tables with literature-typical
% features (suction peak, separation, base pressure) - they are NOT a
% prediction of what your rig will give.
rng(7);
% anchors: [theta(deg) Cp]
A = { ...
 [0 1; 15 0.62; 30 0.05; 45 -0.65; 60 -1.42; 70 -2.05; 76 -2.32; ...
  80 -2.22; 85 -1.62; 90 -1.40; 110 -1.31; 180 -1.28], ...           % smooth low
 [0 1; 15 0.62; 30 0.10; 45 -0.62; 60 -1.38; 72 -1.95; 79 -2.10; ...
  83 -1.98; 87 -1.62; 91 -1.32; 96 -1.16; 110 -1.08; 180 -1.06], ... % smooth high
 [0 1; 15 0.50; 30 -0.15; 45 -0.85; 60 -1.30; 80 -1.55; 92 -1.58; ...
  102 -1.45; 112 -1.20; 120 -0.80; 128 -0.62; 180 -0.60], ...        % trip low
 [0 1; 15 0.50; 30 -0.20; 45 -0.85; 60 -1.28; 80 -1.46; 95 -1.50; ...
  105 -1.35; 115 -1.10; 124 -0.85; 132 -0.76; 180 -0.74]};           % trip high
labs = {'Smooth - low speed','Smooth - high speed','Trip - low speed','Trip - high speed'};
rpmv = [700 1350 700 1350];
trip = [false false true true];
Uv   = [15.5 30.5 15.5 30.5];

runs = struct('label',{},'tripped',{},'rpm',{},'Tair',{},'pBaro',{},'h',{},'hPitot',{});
for k = 1:4
    Tair  = 21.0 + 0.5*randn;
    pBaro = 101300 + 200*randn;
    [rho, ~] = airProps(Tair, pBaro, cfg);
    dd = cfg.rhoEth*cfg.g;                             % N/m^3 per metre
    q = 0.5*rho*Uv(k)^2;
    cp = synthCp(thetaHole, A{k});
    hStat = 250;                                       % arbitrary datum, mm
    h = hStat + cp*q/dd*1e3 + 0.4*randn(1,28);         % 0.4 mm reading noise
    runs(end+1) = struct('label',['DEMO - ' labs{k}],'tripped',trip(k),'rpm',rpmv(k), ...
        'Tair',Tair,'pBaro',pBaro,'h',h, ...
        'hPitot',[hStat, hStat + q/dd*1e3]);
end
% 5th run left empty on purpose: demonstrates the graceful skip
runs(end+1) = struct('label','Trip - extra speed (choose with demonstrator)', ...
    'tripped',true,'rpm',NaN,'Tair',[],'pBaro',[],'hPitot',[],'h',[]);
end

% -----------------------------------------------------------------------
function cp = synthCp(thDeg, anchors)
th = thDeg(:);
% negative angles (hole 1) are mirrored, as a symmetric flow would give
cp = interp1(anchors(:,1), anchors(:,2), abs(th), 'pchip')';
end

% -----------------------------------------------------------------------
function u = nanUnc(cfg)
% Same fields as cpUncertainty, all empty/NaN. Keeps every print and CSV
% column valid when the module or the repeat readings are unavailable.
u = struct('uCp',[],'U95Cp',[],'uCpPeak',NaN,'uCD',NaN,'uCDrel',NaN, ...
           'uThetaSep',NaN,'K',cfg.sepUncK,'uTap',[],'uStat',NaN, ...
           'uTot',NaN,'nFine',cfg.nFine);
end

% -----------------------------------------------------------------------
function [runs, uM] = loadReadings(file, cfg)
%LOUDREADINGS  Build runs and repeat-reading uncertainties from readings.csv.
%
% The log is WIDE: one row per hole, the repeat readings side by side, plus
% the room conditions (also read three times) and the board scale division:
%
%   run,tripped,rpm,scale_mm,Tair_C_1,Tair_C_2,Tair_C_3, ...
%       pBaro_Pa_1,pBaro_Pa_2,pBaro_Pa_3,hole,r_1,r_2,r_3
%   1,0,700,1.0,21.0,21.1,20.9,101300,101299,101301,1,268.4,268.5,
%   1,0,700,1.0,,,,,,,,2,268.9,268.8,
%
%   hole      1..28 = tapping, 35/36 = the two Pitot tubes
%   r_1..r_3  repeats of the same reading; leave the unused ones blank
%   scale_mm  one small division of the board [mm] -> type B term
%             u_B = (scale_mm/2)/sqrt(3); blank falls back to cfg.uScale
%   Tair_C_i / pBaro_Pa_i may be filled once or three times per run; all
%   non-blank values are averaged. Blank -> the cfg defaults.
%
% Every (run, hole) row is averaged over its non-blank repeats. With n >= 2
% the standard uncertainty of the mean, s/sqrt(n), is the type A term and is
% combined in quadrature with u_B. A single reading carries no type A
% information, so it falls back to the conservative cfg.uRead.
%
% A run is skipped (and reported) unless all 28 holes and both Pitot tubes
% have at least one reading.

T = readtable(file, 'Delimiter', ',', 'VariableNamingRule', 'preserve');

names = {'Smooth - low speed','Smooth - high speed','Trip - low speed', ...
         'Trip - high speed','Trip - extra speed (choose with demonstrator)'};

uB   = (cfg.uScale/2)/sqrt(3);          % type B, uniform distribution
runs = struct('label',{},'tripped',{},'rpm',{},'Tair',{},'pBaro',{}, ...
              'h',{},'hPitot',{});
uM   = struct('uTap',{},'u35',{},'u36',{});

ids = unique(T.run(~isnan(T.run)))';
for k = 1:numel(ids)
    id  = ids(k);
    sel = (T.run == id);

    % the scale division can be given per run; blank -> cfg.uScale
    scale = meanValid(T.scale_mm(sel), cfg.uScale);
    uB    = (scale/2)/sqrt(3);

    hv  = cell(1,28);  hv(:) = {[]};    % readings per tapping
    v35 = [];  v36 = [];

    for i = find(sel)'
        hh = T.hole(i);
        if isnan(hh), continue; end
        rr = [T.r_1(i) T.r_2(i) T.r_3(i)];
        rr = rr(~isnan(rr));            % unused repeat columns stay blank
        if isempty(rr), continue; end
        if hh >= 1 && hh <= 28
            hv{hh} = [hv{hh}; rr(:)];
        elseif hh == 35
            v35 = [v35; rr(:)];
        elseif hh == 36
            v36 = [v36; rr(:)];
        end
    end

    if isempty(v35) || isempty(v36)
        fprintf('[skip] run %d: no Pitot readings in %s\n', id, file);
        continue
    end

    h  = NaN(1,28);
    uT = NaN(1,28);
    ok = true;
    for j = 1:28
        v = hv{j};
        if isempty(v), ok = false; break; end
        h(j) = mean(v);
        if numel(v) >= 2
            uT(j) = sqrt((std(v)/sqrt(numel(v)))^2 + uB^2);
        else
            % a single reading carries no type A information; fall back to
            % the conservative value rather than quoting the type B alone
            uT(j) = cfg.uRead;
        end
    end
    if ~ok
        fprintf('[skip] run %d: at least one tapping has no reading in %s\n', ...
                id, file);
        continue
    end

    % metadata: average every non-blank entry, blanks fall back to cfg
    ta  = [T.Tair_C_1(sel);   T.Tair_C_2(sel);   T.Tair_C_3(sel)];
    pa  = [T.pBaro_Pa_1(sel); T.pBaro_Pa_2(sel); T.pBaro_Pa_3(sel)];
    rp  = meanValid(T.rpm(sel), NaN);
    tr  = meanValid(T.tripped(sel), 0);

    if id <= numel(names)
        lbl = names{id};
    else
        lbl = sprintf('Run %d', id);
    end

    runs(end+1) = struct('label',lbl,'tripped',logical(tr),'rpm',rp, ...
        'Tair',meanValid(ta,NaN),'pBaro',meanValid(pa,NaN), ...
        'h',h,'hPitot',[mean(v35) mean(v36)]);

    uM(end+1) = struct('uTap',uT, ...
        'u35',pitotU(v35,uB,cfg.uRead),'u36',pitotU(v36,uB,cfg.uRead));
end
end

% -----------------------------------------------------------------------
function u = pitotU(v, uB, uFallback)
% Type A + type B uncertainty of one Pitot reading; a single reading gets
% the conservative fallback (see loadReadings).
if numel(v) >= 2
    u = sqrt((std(v)/sqrt(numel(v)))^2 + uB^2);
else
    u = uFallback;
end
end

% -----------------------------------------------------------------------
function x = meanValid(v, dflt)
% Mean of the non-blank entries; dflt when every entry is blank.
v = v(~isnan(v));
if isempty(v)
    x = dflt;
else
    x = mean(v);
end
end
