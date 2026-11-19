function T = run_paper_eps_dos_check(variant)
%RUN_PAPER_EPS_DOS_CHECK Absorber permittivity and density-of-states check (round-5 review R2): nominal devices (HJ-0, F-HJ, RBSF-GaP, RHJ),
% AM1.5G and CIE 6500 K/1000 lx, production sweeps and unified extractor, with
%   variant 'nominal' : eps = 12.9, Nc = Nv = 1e19 (production)
%   variant 'eps'     : eps = 11.1 (Vegard interpolation for GaN0.03P0.924As0.046 with GaP 11.1, GaAs 12.9, GaN 9.7)
%   variant 'dos'     : Nc = 1.8e19, Nv = 1.9e19 (GaP values), B_rad recomputed by van Roosbroeck-Shockley for these Nc, Nv
% Output: results/paper_eps_dos_check/<variant>.csv
assertCanonicalPath(); warning('off','all');
root = fileparts(fileparts(mfilename('fullpath'))); outDir = fullfile(root,'results','paper_eps_dos_check'); if ~isfolder(outDir), mkdir(outDir); end
arc = load(fullfile(root,'results','paper_stageA','stageA_comparison.mat'),'d1nm','d2nm');
arcFn = @(lam) tmmReflectance(lam,[1.0,1.38,2.35,gapRefractiveIndexReal(lam)],[arc.d1nm,arc.d2nm]);
Bmat = paperRadiativeCoefficients(300); Babs = Bmat.GaAsPN_cm3ps; extra = {};
switch variant
  case 'nominal'
  case 'eps', extra = {'AbsorberEps',11.1};
  case 'dos'
    Babs = computeBradVRS(loadGaAsNPOptical(3.0),gaasnpBandgap(0.03).Eg,1.8e19,1.9e19,300);
    extra = {'AbsorberNc',1.8e19,'AbsorberNv',1.9e19};
  otherwise, error('unknown variant');
end
specs = {loadRealAM15G(), buildCIEIndoorSpectrum('6500K',1000)}; names = {'AM1.5G','6500K_1000lx'};
Vg = {(0:0.02:1.70)', (0:0.02:1.56)'}; rows = {};
for code = {'HJ-0','F-HJ','RBSF-GaP','RHJ'}
  dev = buildPaperArchitecture(code{1},'BradAbs',Babs,'BradGaP',Bmat.GaP_cm3ps,'BradGaNP',Bmat.GaPN_Ege18_cm3ps, ...
    'CnAugerAbs',3e-30,'CpAugerAbs',3e-30,'CnAugerGaP',3e-30,'CpAugerGaP',3e-30,'CnAugerGaNP',3e-30,'CpAugerGaNP',3e-30, extra{:});
  for ic = 1:2
    G = deviceGenAndTransmitML(dev,specs{ic},arcFn); Pin = spectrumPower(specs{ic});
    if ic == 1
      [~,~,m] = runAcceptedJVSweep(dev,G,Vg{ic},Pin,true);
    else
      gate = sweepIndoor(dev,G,Vg{ic},Pin); fU = gate.isFinite(:); m = pvMetrics([gate.Vs(fU), gate.Jterm(fU)],Pin);
    end
    rows(end+1,:) = {variant, code{1}, names{ic}, m.Jsc, m.Voc*1e3, m.FF*100, m.PCE}; %#ok<AGROW>
    fprintf('%s %s %s PCE=%.4f\n', variant, code{1}, names{ic}, m.PCE);
  end
end
T = cell2table(rows,'VariableNames',{'Variant','Architecture','Condition','Jsc','Voc_mV','FF_pct','PCE_pct'});
writetable(T,fullfile(outDir,[variant '.csv']));
end

function gate = sweepIndoor(dev,G,Vgrid,Pin)
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
end
