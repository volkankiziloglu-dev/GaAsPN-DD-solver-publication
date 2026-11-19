function analyze_sukrittanon_scaledSRH()
%ANALYZE_SUKRITTANON_SCALEDSRH Fit k*SRH + free-n diode to ALL dark points
% (Rs-corrected), then predict the light curve: J = J_light - (k-1) J_SRH - J_leak.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'src'),'-begin');
rd = fullfile(root,'results','sukrittanon_ganp');
Vt = 0.025852; Rs = 22; Pin = spectrumPower(loadRealAM15G());
names = {'GaNP-800','GaP-800','GaNP-2000'};
dark = {[0.60 -6.45; 0.80 -5.0; 0.90 -3.7; 1.00 -2.6; 1.10 -1.2; 1.20 -0.23; 1.30 1.35; 1.34 1.95], ...
        [0.60 -8.3; 0.70 -7.0; 0.80 -5.9; 0.90 -4.9; 1.00 -3.0; 1.10 -1.5; 1.20 0.3; 1.29 2.0], ...
        [0.90 -7.3; 1.00 -6.0; 1.05 -5.3; 1.10 -4.4; 1.15 -3.2; 1.20 -2.0; 1.25 -1.2; 1.30 -0.55; 1.35 -0.3]};
ex = {[5.52 1260 73 5.1],[2.69 1190 77 2.5],[7.80 1310 71 7.3]};
lt = {'GaNP-800_tau5e-09_S0_sh0.10_mu1.00_Na0_AR0_tB1e-09_tENaN_tG5e-09_final.csv', ...
      'GaP-800_tau5e-09_S0_sh0.10_mu1.00_Na0_AR0_tB1e-09_tENaN_tG5e-09_final.csv', ...
      'GaNP-2000_tau5e-09_S0_sh0.10_mu1.00_Na0_AR0_tB1e-09_tENaN_tG5e-09_final.csv'};
fprintf('cell | k (=5ns/tau_eff) | tau_eff ns | J0p, n_p | rms ln J (all dark) | Jsc Voc FF PCE (model) | exp\n');
for c = 1:3
  D = readtable(fullfile(rd, ['dark_' names{c} '_tau5e-09.csv']));
  JS = @(V) exp(interp1(D.V, log(D.Jdark_Acm2), V, 'linear', 'extrap'));
  L = readtable(fullfile(rd, lt{c}));
  e = dark{c};
  lk = @(p,V) exp(p(2)).*(exp(V./(p(3)*Vt))-1);
  tot = @(p,V) exp(p(1))*JS(V) + lk(p,V);
  obj = @(p) sum((lnD(tot,p,e(:,1),Rs) + log(1e3) - e(:,2)).^2);
  best = []; bf = Inf;
  for k0 = [-2 0 2]
    for n0 = [2.5 3.5]
      p = fminsearch(obj,[k0 log(1e-9) n0],optimset('MaxFunEvals',6000,'TolX',1e-8,'TolFun',1e-10));
      if p(3) > 1 && p(3) < 8 && obj(p) < bf, bf = obj(p); best = p; end
    end
  end
  p = best; k = exp(p(1));
  Jt = L.J_Acm2 - (k-1)*JS(L.V) - lk(p,L.V);
  Vterm = L.V - Jt*Rs; [Vterm,ia] = unique(Vterm); Jt = Jt(ia);
  m = pvMetrics([Vterm Jt], Pin);
  fprintf('%s | k=%.3g | tau_eff=%.2f ns | J0p=%.2e n_p=%.2f | rms=%.2f | %.2f %.0f %.1f %.2f | exp %.2f %.0f %.0f %.1f\n', ...
    names{c}, k, 5/k, exp(p(2)), p(3), sqrt(bf/size(e,1)), m.Jsc*1e3, m.Voc*1e3, m.FF*100, m.PCE, ex{c});
end
end

function lnJ = lnD(tot, p, Vterm, Rs)
Vint = Vterm(:);
for it = 1:60
  J = tot(p,Vint); Vint = 0.5*Vint + 0.5*(Vterm(:) - J*Rs);
end
lnJ = log(tot(p,Vint));
end
