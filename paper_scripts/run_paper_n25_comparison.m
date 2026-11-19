function T = run_paper_n25_comparison(Npct, condSet)
%RUN_PAPER_N25_COMPARISON Composition check: architecture comparison with a 2.5% N absorber.
% The Si-matched GaAsPN absorber is changed from 3% N (Eg 1.857 eV) to 2.5% N (Eg 1.921 eV, BAC, close to the
% 1.90 eV of GaInP). Consistently changed: electrical bandgap, Zou-Goodnick n,k table (nk_25N), radiative
% coefficient (van Roosbroeck-Shockley from that table and gap, Nc = Nv = 1e19) and electron affinity (gap change
% taken in the conduction band, chi_abs = 4.0 - dEg). Everything else nominal; the ARC is kept at the Stage-A
% design (not re-optimized). Npct = 2.5 (check) or 3.0 (control, must reproduce production).
% condSet 'main': AM1.5G, CIE 6500 K at 1000 lx and 100 lx; 'spectra': CIE 2700 K, 4000 K, FL11 at 1000 lx.
% Output: results/paper_n25_comparison/n<Npct*10>_<condSet>.csv
assertCanonicalPath(); warning('off','all');
root = fileparts(fileparts(mfilename('fullpath'))); outDir = fullfile(root,'results','paper_n25_comparison');
if ~isfolder(outDir), mkdir(outDir); end
arc = load(fullfile(root,'results','paper_stageA','stageA_comparison.mat'),'d1nm','d2nm');
arcFn = @(lam) tmmReflectance(lam,[1.0,1.38,2.35,gapRefractiveIndexReal(lam)],[arc.d1nm,arc.d2nm]);
Bmat = paperRadiativeCoefficients(300);
absBand = gaasnpBandgap(Npct/100);
optAbs = loadGaAsNPOptical(Npct);
Babs = computeBradVRS(optAbs,absBand.Eg,1e19,1e19,300);
switch condSet
  case 'main'
    conds = {'AM1.5G','6500K_1000lx','6500K_100lx'};
    specs = {loadRealAM15G(), buildCIEIndoorSpectrum('6500K',1000), buildCIEIndoorSpectrum('6500K',100)};
    Vgr = {(0:0.02:1.76)', (0:0.02:1.60)', (0:0.01:1.60)'};
  case 'spectra'
    conds = {'2700K_1000lx','4000K_1000lx','4000F_1000lx'};
    specs = {buildCIEIndoorSpectrum('2700K',1000), buildCIEIndoorSpectrum('4000K',1000), buildCIEIndoorSpectrum('4000F',1000)};
    Vgr = {(0:0.02:1.60)', (0:0.02:1.60)', (0:0.02:1.60)'};
  case 'cwled'   % Klitzke et al. CW-LED, scaled to their 1000 lx and 100 lx input powers
    Tk = readtable(resolveDataFile('Klitzke2025_CWLED_SPD_digitized.csv'));
    s0 = [Tk.lambda_nm, Tk.relative_irradiance];
    conds = {'CWLED_3.09W','CWLED_0.31W'};
    specs = {[s0(:,1), s0(:,2)*(3.09/spectrumPower(s0))], [s0(:,1), s0(:,2)*(0.31/spectrumPower(s0))]};
    Vgr = {(0:0.02:1.56)', (0:0.01:1.56)'};
end
Pins = cellfun(@spectrumPower,specs);
common = {'AbsorberNFraction',Npct/100,'AbsorberOpticalN',Npct,'BradAbs',Babs,'BradGaP',Bmat.GaP_cm3ps,'BradGaNP',Bmat.GaPN_Ege18_cm3ps, ...
  'CnAugerAbs',3e-30,'CpAugerAbs',3e-30,'CnAugerGaP',3e-30,'CpAugerGaP',3e-30,'CnAugerGaNP',3e-30,'CpAugerGaNP',3e-30};
codes = {'HJ-0','F-HJ','RBSF-GaP','RHJ'};
fprintf('N = %.1f%%: Eg = %.4f eV, chi_abs = %.4f eV, B_rad = %.4g cm3/s\n',Npct,absBand.Eg,4.0-(absBand.Eg-gaasnpBandgap(0.03).Eg),Babs);
rows = {}; ck = fullfile(outDir,sprintf('n%02d_%s.csv',round(Npct*10),condSet));
for icd = 1:numel(conds)
  for ic = 1:4
    code = codes{ic}; t0 = tic;
    try
      dev = buildPaperArchitecture(code,common{:});
      G = deviceGenAndTransmitML(dev,specs{icd},arcFn);
      if strcmp(conds{icd},'AM1.5G')
        [~,~,m] = runAcceptedJVSweep(dev,G,Vgr{icd},Pins(icd),true); cls = 'AM15G';
        Jsc = m.Jsc*1e3; Voc = m.Voc*1e3; FF = m.FF*100; PCE = m.PCE;   % Jsc in mA/cm2
      else
        gate = indoorSweepLocal(dev,G,Vgr{icd},Pins(icd)); cls = char(gate.curveClass);
        Jsc = gate.Jsc*1e6; Voc = gate.Voc*1e3; FF = gate.FF*100; PCE = gate.PCE;   % Jsc in uA/cm2
      end
    catch ME
      Jsc = NaN; Voc = NaN; FF = NaN; PCE = NaN; cls = ['ERR ' strrep(ME.message,',',';')];
    end
    rows(end+1,:) = {code, Npct, absBand.Eg, Babs, conds{icd}, Pins(icd), Jsc, Voc, FF, PCE, cls}; %#ok<AGROW>
    T = cell2table(rows,'VariableNames',{'Architecture','N_pct','Eg_eV','BradAbs','Condition','Pin_Wm2','Jsc','Voc_mV','FF_pct','PCE_pct','Class'});
    writetable(T,ck);
    fprintf('N=%.1f %s %s: Jsc=%.4f Voc=%.1f FF=%.2f PCE=%.3f %s (%.0f s)\n',Npct,code,conds{icd},Jsc,Voc,FF,PCE,cls,toc(t0));
  end
end
end

function gate = indoorSweepLocal(dev,G,Vgrid,Pin)
% Same sweep, gate and unified extractor as run_paper_eg_optics_shift.
N = dev.N; idxStart = find(abs(Vgrid-0.02)<1e-12,1);
r002 = gummelSolveQFPAnalytic(dev,0.02,G,struct('psiInit',dev.psiEq,'phinInit',zeros(N,1),'phipInit',zeros(N,1)));
prev = r002;
for v = [0.01,0.005,0], r0 = gummelSolveQFPAnalytic(dev,v,G,struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip)); prev = r0; end
states = repmat(r002,numel(Vgrid),1); states(1) = r0; states(idxStart) = r002; prev = r002;
for k = idxStart+1:numel(Vgrid)
  rk = gummelSolveQFPAnalytic(dev,Vgrid(k),G,struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip));
  states(k) = rk; if rk.isFinite, prev = rk; end
end
gate = indoorAcceptanceGate(dev,G,Vgrid,states,Pin);
fU = gate.isFinite(:); mU = pvMetrics([gate.Vs(fU), gate.Jterm(fU)],Pin);
gate.Jsc = mU.Jsc; gate.Voc = mU.Voc; gate.FF = mU.FF; gate.PCE = mU.PCE;
end
