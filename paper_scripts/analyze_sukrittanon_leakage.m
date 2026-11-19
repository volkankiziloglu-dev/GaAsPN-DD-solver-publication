function analyze_sukrittanon_leakage()
%ANALYZE_SUKRITTANON_LEAKAGE Add a parallel leakage branch to the intrinsic
% (tau_i = 5 ns) model of the Sukrittanon cells, fit it to the DARK J-V
% only (V <= 1.1 V, Rs-corrected), then test it on held-out data: dark J-V
% above 1.1 V and the light Voc / FF / PCE.  No light data enter the fit.
%
% Branch forms:  ohmic shunt  J = V/Rsh;  diode with n_p = 2;  diode with free n_p.
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'src'),'-begin');
rd = fullfile(root,'results','sukrittanon_ganp');
Vt = 0.025852; Rs = 22;
Pin = spectrumPower(loadRealAM15G());

% digitised dark data, [V, ln(J/mA cm-2)] (thesis Fig. 4.7a, +-0.3)
cells = struct( ...
 'name',{'GaNP-800','GaP-800','GaNP-2000'}, ...
 'dark',{[0.60 -6.45; 0.80 -5.0; 0.90 -3.7; 1.00 -2.6; 1.10 -1.2; 1.20 -0.23; 1.30 1.35; 1.34 1.95], ...
         [0.60 -8.3; 0.70 -7.0; 0.80 -5.9; 0.90 -4.9; 1.00 -3.0; 1.10 -1.5; 1.20 0.3; 1.29 2.0], ...
         [0.90 -7.3; 1.00 -6.0; 1.05 -5.3; 1.10 -4.4; 1.15 -3.2; 1.20 -2.0; 1.25 -1.2; 1.30 -0.55; 1.35 -0.3]}, ...
 'exp',{[5.52 1260 73 5.1],[2.69 1190 77 2.5],[7.80 1310 71 7.3]}, ...
 'light',{'GaNP-800_tau5e-09_S0_sh0.10_mu1.00_Na0_AR0_tB1e-09_tENaN_tG5e-09_final.csv', ...
          'GaP-800_tau5e-09_S0_sh0.10_mu1.00_Na0_AR0_tB1e-09_tENaN_tG5e-09_final.csv', ...
          'GaNP-2000_tau5e-09_S0_sh0.10_mu1.00_Na0_AR0_tB1e-09_tENaN_tG5e-09_final.csv'});
Vfit = 1.10;
out = {};
for c = 1:numel(cells)
  C = cells(c);
  D = readtable(fullfile(rd, ['dark_' C.name '_tau5e-09.csv']));
  lnSRH = @(V) interp1(D.V, log(D.Jdark_Acm2), V, 'linear', 'extrap');   % ln A/cm2
  L = readtable(fullfile(rd, C.light));                                    % V, J (A/cm2), intrinsic
  forms = {'none','ohmic shunt','diode n=2','diode free n'};
  for f = 1:numel(forms)
    form = forms{f};
    if strcmp(form,'none')
      p = [];  rms = NaN;
      leak = @(V) 0*V;
    else
      sel = C.dark(:,1) <= Vfit;
      Ve = C.dark(sel,1); lnJe = C.dark(sel,2);
      switch form
        case 'ohmic shunt', mk = @(p) @(V) V./exp(p(1));            % p = ln Rsh
        case 'diode n=2',   mk = @(p) @(V) exp(p(1)).*(exp(V./(2*Vt))-1);
        case 'diode free n',mk = @(p) @(V) exp(p(1)).*(exp(V./(p(2)*Vt))-1);
      end
      if strcmp(form,'ohmic shunt'), p0 = log(1e5); else, p0 = [log(1e-9) 3]; end
      if strcmp(form,'diode n=2'), p0 = log(1e-9); end
      obj = @(p) sum((lnDark(mk(p), lnSRH, Ve, Rs) + log(1e3) - lnJe).^2);
      p = fminsearch(obj, p0, optimset('TolX',1e-8,'TolFun',1e-10,'MaxFunEvals',4000));
      leak = mk(p);
      rms = sqrt(obj(p)/numel(Ve));
    end
    % held-out dark (V > Vfit)
    hv = C.dark(C.dark(:,1) > Vfit, :);
    dres = NaN;
    if ~isempty(hv)
      dres = sqrt(mean((lnDark(leak, lnSRH, hv(:,1), Rs) + log(1e3) - hv(:,2)).^2));
    end
    % light: J = Jlight_int(Vint) - leak(Vint); terminal V = Vint - J*Rs
    Jt = L.J_Acm2 - leak(L.V);
    Vterm = L.V - Jt*Rs;
    [Vterm, ia] = unique(Vterm); Jt = Jt(ia);
    m = pvMetrics([Vterm Jt], Pin);
    out(end+1,:) = {C.name, form, rms, dres, m.Jsc*1e3, m.Voc*1e3, m.FF*100, m.PCE, C.exp(1), C.exp(2), C.exp(3), C.exp(4), p}; %#ok<AGROW>
  end
end
T = cell2table(out(:,1:12),'VariableNames',{'cell','branch','rmsLowV_lnJ','rmsHeldOutDark_lnJ','Jsc','Voc_mV','FF','PCE','Jsc_exp','Voc_exp','FF_exp','PCE_exp'});
disp(T);
writetable(T, fullfile(rd,'leakage_test.csv'));
for k = 1:size(out,1)
  if ~isempty(out{k,13}), fprintf('%s | %s | p = %s\n', out{k,1}, out{k,2}, mat2str(out{k,13},4)); end
end
end

function lnJ = lnDark(leak, lnSRH, Vterm, Rs)
% ln(J_dark [A/cm2]) at terminal voltage Vterm, with Rs: V_int = V - J*Rs.
Vint = Vterm(:);
for it = 1:60
  J = exp(lnSRH(Vint)) + leak(Vint);
  VintNew = Vterm(:) - J*Rs;
  Vint = 0.5*Vint + 0.5*VintNew;
end
J = exp(lnSRH(Vint)) + leak(Vint);
lnJ = log(J);
end
