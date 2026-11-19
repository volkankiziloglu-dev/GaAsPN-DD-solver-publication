function P = fit_klitzke_hj_branch_v112()
%FIT_KLITZKE_HJ_BRANCH_V112 Refit of the light-independent leakage branch of the Klitzke homojunction with this project's HJ model.
% Branch: J_leak(V) = J0p*(exp(Vd/(n*VT))-1), Vd = V - J_leak*Rp (diode plus series resistor), subtracted from the intrinsic model:
% J(V) = J_model(V) - J_leak(V). Fitted at 1000 lx (as the earlier fit, Klitzke_Validasyon/cases/klitzke_branch_vs_rsrsh_eval.js:
% J0p = 4.97e-10 A/cm2, n = 3.72, Rp = 8.82e3 Ohm cm2 on the Klitzke_Validasyon HJ model) to the measured curve (V >= 0.3 V,
% J > 0.02 Jsc, ~80 points, rms error in % of the measured Jsc), then evaluated at 100 lx.
% Model curves: results/paper_klitzke_one_solver/HJ_<lux>lx_mesh1_qfp.csv (dd_solver v1.0.0 = v1.1.2 production).
% Output: results/paper_klitzke_rhj_excess_v112/hj_branch_fit.csv
root = fileparts(fileparts(mfilename('fullpath')));
Tm = readtable(resolveDataFile(fullfile('experimental','klitzke2025','klitzke2025_fig4_iv_curves.csv')));
VT = 0.025852; luxs = [1000 100]; D = struct();
for l = luxs
  c = readtable(fullfile(root,'results','paper_klitzke_one_solver',sprintf('HJ_%dlx_mesh1_qfp.csv',l)));
  c = c(c.Finite==1,:); mo = [c.V_V, c.J_Acm2*1e6];
  sel = strcmp(Tm.device,'Homojunction') & Tm.lux == l; me = sortrows([Tm.V_V(sel), Tm.J_uAcm2(sel)]);
  Jsc = mean(me(me(:,1)>0.06 & me(:,1)<0.3,2)); s = me(me(:,1)>=0.3 & me(:,2)>0.02*Jsc,:); st = max(1,floor(size(s,1)/80)); s = s(1:st:end,:);
  D.(sprintf('L%d',l)) = struct('mo',mo,'p',s,'Jsc',Jsc);
end
leak = @(V,J0,n,Rp) arrayfun(@(v) leak1(v,J0,n,Rp,VT), V);
% Fit on the excess recombination current dJ = Jrec,meas - Jrec,model (Jrec = Jsc - J), so that the Jsc offset between model and
% digitized curve does not enter (as for the RHJ table); n bounded to [1, 6] by a logistic map.
for l = luxs
  d = D.(sprintf('L%d',l)); JscD = d.mo(1,2);
  d.dJ = (d.Jsc - d.p(:,2)) - (JscD - interp1(d.mo(:,1),d.mo(:,2),d.p(:,1),'linear','extrap')); D.(sprintf('L%d',l)) = d;
end
nmap = @(u) 1 + 5./(1+exp(-u)); umap = @(n) -log(5./(n-1) - 1);
cost = @(x,d) 100*sqrt(mean((leak(d.p(:,1),10^x(1),nmap(x(2)),10^x(3)) - d.dJ).^2))/d.Jsc;
x0 = [log10(4.97e-4) umap(3.72) log10(8.82e3)];   % J0p in uA/cm2
c0 = cost(x0, D.L1000);
% joint fit at 100 and 1000 lx (the branch is light-independent by construction); a 1000 lx-only fit runs to n = 1
x = fminsearch(@(x) hypot(cost(x,D.L1000), cost(x,D.L100)), x0, optimset('MaxFunEvals',6000,'MaxIter',6000,'TolX',1e-7,'TolFun',1e-9));
P = table(10^x(1)*1e-6, nmap(x(2)), 10^x(3), cost(x,D.L1000), cost(x,D.L100), c0, cost(x0,D.L100), ...
  'VariableNames',{'J0p_Acm2','n','Rp_Ohmcm2','rms1000_pct','rms100_pct','rms1000_oldParams_pct','rms100_oldParams_pct'});
V = (0.6:0.1:1.2)'; E = table(V, interp1(D.L1000.p(:,1),D.L1000.dJ,V,'linear',NaN), interp1(D.L100.p(:,1),D.L100.dJ,V,'linear',NaN), ...
  leak(V,10^x(1),nmap(x(2)),10^x(3)), leak(V,4.97e-4,3.72,8.82e3), 'VariableNames',{'V_V','dJ1000_uA','dJ100_uA','newBranch_uA','oldBranch_uA'});
disp(E); writetable(E, fullfile(root,'results','paper_klitzke_rhj_excess_v112','hj_excess_current.csv'));
disp(P); writetable(P, fullfile(root,'results','paper_klitzke_rhj_excess_v112','hj_branch_fit.csv'));
end
function J = leak1(V,J0,n,Rp,VT)
if V <= 0, J = 0; return; end
lo = 0; hi = V;
for k = 1:80, m = (lo+hi)/2; f = m + Rp*1e-6*J0*(exp(min(m/(n*VT),700))-1) - V; if f > 0, hi = m; else, lo = m; end, end
J = J0*(exp(min(lo/(n*VT),700))-1);
end
