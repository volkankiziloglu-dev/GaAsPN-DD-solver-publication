function T = run_paper_leakage_archspecific()
%RUN_PAPER_LEAKAGE_ARCHSPECIFIC Architecture-specific leakage test (referee check).
% The earlier test (run_paper_leakage_robustness.m) subtracted a branch fitted to the measured excess current of the Klitzke HOMOjunction from all four
% architectures, although the excess current of the Klitzke RHJ behaves differently (sublinear in light, supplement, Table "additional current of the
% lower-doped RHJ"). Here the HJ-type architectures (HJ-0, F-HJ, RBSF-GaP) receive the homojunction branch (J0p x0.5/x1/x2) and the RHJ receives the
% measured RHJ excess current dJ(V) (tabulated, log-linear interpolation/extrapolation in V, same x0.5/x1/x2 scaling), for 1000 and 100 lx.
% Post-processing of the intrinsic production curves, as in the earlier test; sensitivity bound, not a prediction.
% Output: results/paper_leakage_robustness/leakage_archspecific.csv.
assertCanonicalPath();
root = fileparts(fileparts(mfilename('fullpath'))); inDir = fullfile(root,'results','paper_indoor_production');
outDir = fullfile(root,'results','paper_leakage_robustness');
S = readtable(fullfile(inDir,'indoor_production_summary.csv'),VariableNamingRule='preserve');
codes = ["HJ-0","F-HJ","RBSF-GaP","RHJ"]; luxs = [1000 100];
% measured excess current of the lower-doped Klitzke RHJ (uA/cm2), supplement Table (additional current)
Vt = [0.80 0.90 1.00 1.10 1.20];
dJ = struct('lux1000',[0.99 1.32 1.69 2.93 5.91], 'lux100',[0.12 0.23 0.38 0.61 2.00]);
rows = {};
scales = [0 0.5 1 2];
for lux = luxs
  tab = dJ.(sprintf('lux%d',lux));
  for ia = 1:numel(codes)
    q = S(strcmp(S.Architecture,codes(ia)) & strcmp(S.Spectrum,'6500K') & S.Lux==lux,:); Pin = q.Pin_Wm2; cls = q.CurveClass{1};
    c = readtable(fullfile(inDir,sprintf('%s_6500K_%dlux_jv.csv',codes(ia),lux)),VariableNamingRule='preserve');
    keep = isfinite(c.V_V) & isfinite(c.J_Acm2); V0 = c.V_V(keep); J0 = c.J_Acm2(keep);
    V = V0; J = J0;   % unified: native 0.02 V grid
    for sc = scales
      for mode = ["archspecific","rhj_branch_all","hj_branch_all"]
        if sc == 0 && mode ~= "archspecific", continue; end
        useRHJ = (mode == "rhj_branch_all") || (mode == "archspecific" && codes(ia) == "RHJ");
        if sc == 0
          Jl = zeros(size(V));
        elseif useRHJ
          lnT = log(tab*1e-6);
          p1 = polyfit(Vt(1:2),lnT(1:2),1); p2 = polyfit(Vt(4:5),lnT(4:5),1);
          lnJ = interp1(Vt,lnT,V,'linear','extrap');
          lnJ(V<Vt(1)) = polyval(p1,V(V<Vt(1))); lnJ(V>Vt(5)) = polyval(p2,V(V>Vt(5)));
          Jl = sc*exp(lnJ); Jl(V<=0) = 0;
        else
          Jl = leakBranch(V,4.97e-10*sc,3.72,8.82e3,0.025852);
        end
        m = pvMetrics([V, J-Jl],Pin);
        rows(end+1,:) = {char(codes(ia)), lux, cls, char(mode), sc, m.Jsc*1e6, m.Voc*1e3, m.FF*100, m.PCE}; %#ok<AGROW>
      end
    end
  end
end
T = cell2table(rows,'VariableNames',{'Architecture','Lux','CurveClass','Mode','Scale','Jsc_uAcm2','Voc_mV','FF_pct','PCE_pct'});
writetable(T,fullfile(outDir,'leakage_archspecific.csv'));
for lux = luxs
  fprintf('\n=== 6500 K, %d lx ===\n',lux);
  for mode = ["archspecific","rhj_branch_all","hj_branch_all"]
    for sc = scales
      if sc == 0 && mode ~= "archspecific", continue; end
      p = zeros(1,4);
      for ia = 1:4
        r = T(strcmp(T.Architecture,codes(ia)) & T.Lux==lux & strcmp(T.Mode,mode) & T.Scale==sc,:); p(ia) = r.PCE_pct;
      end
      [~,o] = sort(p,'descend');
      fprintf('%-15s x%-4g HJ-0 %6.2f F-HJ %6.2f RBSF %6.2f RHJ %6.2f | %s | RBSF-RHJ %+.2f  RHJ-HJ0 %+.2f\n',mode,sc,p(1),p(2),p(3),p(4),strjoin(cellstr(codes(o)),'>'),p(3)-p(4),p(4)-p(1));
    end
  end
end
end

function J = leakBranch(V, J0p, np, Rp, VT)
J = zeros(size(V));
for i = 1:numel(V)
  if V(i) <= 0, continue; end
  f = @(Vd) Vd + Rp*J0p*(exp(min(Vd/(np*VT),700))-1) - V(i);
  lo = 0; hi = V(i);
  for k = 1:80, mid = (lo+hi)/2; if f(mid) > 0, hi = mid; else, lo = mid; end, end
  J(i) = J0p*(exp(min(lo/(np*VT),700))-1);
end
end
