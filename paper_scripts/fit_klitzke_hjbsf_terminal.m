function T = fit_klitzke_hjbsf_terminal()
%FIT_KLITZKE_HJBSF_TERMINAL Can one light-independent loss turn the Klitzke homojunction WITH back-surface field
% (run_klitzke_hj_bsf.m, p-AlInP BSF) into the measured device at 100 lx, 1000 lx and AM1.5D?
% Works on the saved intrinsic curves (no re-solving). Measured metrics: Klitzke et al. Table II.
% Losses tested (all light-independent, applied identically at every condition):
%   RsRsh   series + shunt resistance
%   Rsh     shunt only
%   branch  parallel diode with series resistor, J = J0p*(exp(Vd/(n VT))-1), Vd = V - J*Rp
%           'old'  parameters of the paper (J0p 4.97e-10 A/cm2, n 3.72, Rp 8.82e3 Ohm cm2), not refitted
%           'fit1000' fitted at 1000 lx only, 100 lx and AM1.5D predicted
%           'fit2'    fitted jointly at 100 and 1000 lx, AM1.5D predicted
% Objective per condition: ((Voc-Voc_e)/10 mV)^2 + ((FF-FF_e)/1 point)^2 + ((PCE-PCE_e)/0.5 point)^2 (J_sc not fitted).
% Output: results/paper_klitzke_hj_bsf/terminal_fit.csv
root = fileparts(fileparts(mfilename('fullpath'))); d = fullfile(root,'results','paper_klitzke_hj_bsf');
cn = {'100','1000','AM15D'}; ex = [1.040 64.2 24.7; 1.240 75.5 34.5; 1.440 86.7 17.6];   % Voc (V), FF (%), PCE (%)
Tspd = readtable(resolveDataFile('Klitzke2025_CWLED_SPD_digitized.csv')); led = [Tspd.lambda_nm, Tspd.relative_irradiance];
Pin = [0.31 3.09 spectrumPower(loadRealAM15D())]; Pin(1:2) = arrayfun(@(p) spectrumPower([led(:,1) led(:,2)*p/spectrumPower(led)]), Pin(1:2));
for k = 1:3
  c = readtable(fullfile(d,sprintf('HJBSF_AlInP_%s_curve.csv',cn{k}))); c = c(c.Finite==1,:); C{k} = sortrows([c.V_V c.J_Acm2]); %#ok<AGROW>
end
VT = 0.025852; rows = {};
mods = { 'none', [], [];  'RsRsh', [1 2], [log10(10) 5.5];  'Rsh', [1 2], 5.5; ...
         'branch old', [], [];  'branch fit1000', 2, [log10(4.97e-10) 3.72 log10(8.82e3)];  'branch fit2', [1 2], [log10(4.97e-10) 3.72 log10(8.82e3)];
         'RsRsh all3', [1 2 3], [log10(10) 5.5];  'branch all3', [1 2 3], [log10(4.97e-10) 3.72 log10(8.82e3)] };
opt = optimset('Display','off','TolX',1e-6,'TolFun',1e-8,'MaxFunEvals',4000,'MaxIter',4000);
for i = 1:size(mods,1)
  name = mods{i,1}; use = mods{i,2}; x = mods{i,3};
  if ~isempty(use), x = fminsearch(@(x) obj(x,name,use), x, opt); end
  for k = 1:3
    m = metr(x,name,k);
    rows(end+1,:) = {name, mat2str(use), cn{k}, pars(x,name), m.Voc*1e3, m.FF*100, m.PCE, ex(k,1)*1e3, ex(k,2), ex(k,3), any(use==k)}; %#ok<AGROW>
  end
end
T = cell2table(rows,'VariableNames',{'Loss','FittedOn','Condition','Parameters','Voc_mV','FF_pct','PCE_pct','ExpVoc_mV','ExpFF','ExpPCE','Fitted'});
writetable(T,fullfile(d,'terminal_fit.csv')); disp(T)

  function e = obj(x,name,use)
    e = 0;
    for kk = use
      m = metr(x,name,kk); if ~isfinite(m.PCE), e = e + 1e6; continue; end
      e = e + ((m.Voc-ex(kk,1))/0.010)^2 + ((m.FF*100-ex(kk,2))/1)^2 + ((m.PCE-ex(kk,3))/0.5)^2;
    end
  end
  function m = metr(x,name,k)
    p = C{k};
    if startsWith(name,'RsRsh'), p = applyTerminalResistance(p,10^x(1),10^x(2));
    elseif strcmp(name,'Rsh'), p = applyTerminalResistance(p,0,10^x(1));
    elseif startsWith(name,'branch')
      if strcmp(name,'branch old'), J0 = 4.97e-10; n = 3.72; Rp = 8.82e3; else, J0 = 10^x(1); n = min(max(x(2),1),8); Rp = 10^x(3); end
      p(:,2) = p(:,2) - arrayfun(@(v) leak1(v,J0,n,Rp,VT), p(:,1));
    end
    try, m = pvMetrics(p,Pin(k)); catch, m = struct('Voc',NaN,'FF',NaN,'PCE',NaN); end
  end
  function s = pars(x,name)
    if startsWith(name,'RsRsh'), s = sprintf('Rs=%.3g Rsh=%.3g',10^x(1),10^x(2));
    elseif strcmp(name,'Rsh'), s = sprintf('Rsh=%.3g',10^x(1));
    elseif strcmp(name,'branch old'), s = 'J0p=4.97e-10 n=3.72 Rp=8.82e3';
    elseif startsWith(name,'branch'), s = sprintf('J0p=%.3g n=%.3g Rp=%.3g',10^x(1),min(max(x(2),1),8),10^x(3));
    else, s = '-'; end
  end
end

function J = leak1(V,J0,n,Rp,VT)
if V <= 0, J = 0; return; end
lo = 0; hi = V;
for it = 1:80
  mid = 0.5*(lo+hi);
  if mid + Rp*J0*(exp(min(mid/(n*VT),700))-1) - V > 0, hi = mid; else, lo = mid; end
end
J = J0*(exp(min(lo/(n*VT),700))-1);
end
