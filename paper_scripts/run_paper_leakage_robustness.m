function T = run_paper_leakage_robustness()
%RUN_PAPER_LEAKAGE_ROBUSTNESS Does the indoor architecture ranking survive an UNMODELLED leakage loss?
% Post-processing of the intrinsic production curves (results/paper_indoor_production, CIE 6500 K, 100 and 1000 lx): a parallel leakage
% path is subtracted from each curve, J_t(V) = J(V) - J_leak(V), for
%   (a) ohmic shunt: J_leak = V/Rsh, Rsh = 1e5 ... 1e7 Ohm cm2;
%   (b) diode + series-resistor branch fitted at 1000 lx to the excess current of the measured Klitzke homojunction (klitzke_parasitic_path.m,
%       re-digitized HJ curves, 2026-10-29): V = Vd + Rp*J, J = J0p*(exp(Vd/(np*kT/q)) - 1), J0p = 4.97e-10 A/cm2 (4.97e-4 uA/cm2), np = 3.72,
%       Rp = 8.82e3 Ohm cm2 (and J0p x0.5 / x2). Previous version used a threshold-linear branch (Vth = 0.72 V, Rp = 1.9e4), kept in _eski_esikli/.
% This is a SENSITIVITY BOUND, not a prediction: the leakage is not derived for GaAsPN. Curves are pchip-interpolated to 2 mV before the
% subtraction. 100 lx curves are diagnostic (SCREEN class) in the paper.
assertCanonicalPath();
root = fileparts(fileparts(mfilename('fullpath'))); inDir = fullfile(root,'results','paper_indoor_production');
outDir = fullfile(root,'results','paper_leakage_robustness'); if ~isfolder(outDir), mkdir(outDir); end
S = readtable(fullfile(inDir,'indoor_production_summary.csv'),VariableNamingRule='preserve');
codes = ["HJ-0","F-HJ","RBSF-GaP","RHJ"]; luxs = [1000 100];
cases = {'intrinsic',0,0; 'Rsh 1e7',1e7,0; 'Rsh 3e6',3e6,0; 'Rsh 1e6',1e6,0; 'Rsh 3e5',3e5,0; 'Rsh 1e5',1e5,0; ...
         'branch J0 x0.5',0,0.5; 'branch J0 x1',0,1; 'branch J0 x2',0,2};
rows = {};
for lux = luxs
  for ia = 1:numel(codes)
    q = S(strcmp(S.Architecture,codes(ia)) & strcmp(S.Spectrum,'6500K') & S.Lux==lux,:); Pin = q.Pin_Wm2; cls = q.CurveClass{1};
    c = readtable(fullfile(inDir,sprintf('%s_6500K_%dlux_jv.csv',codes(ia),lux)),VariableNamingRule='preserve');
    keep = isfinite(c.V_V) & isfinite(c.J_Acm2); V0 = c.V_V(keep); J0 = c.J_Acm2(keep);
    V = V0; J = J0;   % unified: native 0.02 V grid, no interpolation (same grid as the main tables)
    for ic = 1:size(cases,1)
      Jl = zeros(size(V));
      if cases{ic,2} > 0, Jl = V/cases{ic,2}; elseif cases{ic,3} > 0, Jl = leakBranch(V, 4.97e-10*cases{ic,3}, 3.72, 8.82e3, 0.025852); end
      m = pvMetrics([V, J-Jl], Pin);
      rows(end+1,:) = {char(codes(ia)), lux, cls, cases{ic,1}, m.Jsc*1e6, m.Voc*1e3, m.FF*100, m.PCE}; %#ok<AGROW>
    end
  end
end
T = cell2table(rows,'VariableNames',{'Architecture','Lux','CurveClass','Leakage','Jsc_uAcm2','Voc_mV','FF_pct','PCE_pct'});
writetable(T,fullfile(outDir,'leakage_robustness.csv'));
for lux = luxs
  fprintf('\n=== 6500 K, %d lx: PCE (%%) per leakage case; rank order and gaps ===\n',lux);
  fprintf('%-18s %8s %8s %9s %8s | %-24s | RBSF-RHJ  RHJ-HJ0  FHJ-HJ0\n','case','HJ-0','F-HJ','RBSF-GaP','RHJ','order');
  for ic = 1:size(cases,1)
    p = zeros(1,4); f = zeros(1,4);
    for ia = 1:4
      r = T(strcmp(T.Architecture,codes(ia)) & T.Lux==lux & strcmp(T.Leakage,cases{ic,1}),:); p(ia) = r.PCE_pct; f(ia) = r.FF_pct;
    end
    [~,o] = sort(p,'descend'); ord = strjoin(cellstr(codes(o)),'>');
    fprintf('%-18s %8.2f %8.2f %9.2f %8.2f | %-24s | %7.2f  %7.2f  %7.2f   (FF %.1f %.1f %.1f %.1f)\n',cases{ic,1},p(1),p(2),p(3),p(4),ord,p(3)-p(4),p(4)-p(1),p(2)-p(1),f);
  end
end
end

function J = leakBranch(V, J0p, np, Rp, VT)
%LEAKBRANCH diode (J0p, np) in series with Rp: V = Vd + Rp*J, J = J0p*(exp(Vd/(np*VT))-1) (A/cm2).
J = zeros(size(V));
for i = 1:numel(V)
  if V(i) <= 0, continue; end
  f = @(Vd) Vd + Rp*J0p*(exp(min(Vd/(np*VT),700))-1) - V(i);
  lo = 0; hi = V(i);
  for k = 1:80, mid = (lo+hi)/2; if f(mid) > 0, hi = mid; else, lo = mid; end, end
  J(i) = J0p*(exp(min(lo/(np*VT),700))-1);
end
end
