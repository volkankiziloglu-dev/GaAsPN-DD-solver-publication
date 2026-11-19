function run_dasilva_pin_qfp()
%RUN_DASILVA_PIN_QFP da Silva et al. 2015's real PIN structure
% (n+GaP/i-GaAsPN/p+GaP, two real measured samples: 1um and 0.3um
% intrinsic-region thickness), rebuilt on the current QFP production
% stack (gummelSolveQFP + homotopy fallback + accepted-gate) for
% consistency with N.1/#49-D1/D2/#58 -- run_gaasnp_pin_dasilva.m (still
% present, untouched) used the old sweepJVML pipeline. Realistic
% interface recombination Sn_int=Sp_int=8910cm/s (Ilahi 2020 measured
% IRV, same value used for MD421/#58) at both n+/i and i/p+ junctions --
% the ideal-interface case was already shown (old solver, see
% run_gaasnp_pin_dasilva.m's own printed history) to hit the same
% Voc/FF=NaN model-validity wall #58's ideal-interface RHJ case did.
% ARC already used here (@arcReflectanceGaP), same as #58/N.1/#49-D1/D2.
assertCanonicalPath();
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:singularMatrix');
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'src'),'-begin');
addpath(fullfile(root,'paper_scripts'),'-begin');
% (publication package: the project tests/ folder is not needed)
resultsDir = fullfile(root,'results');

am15 = loadRealAM15G();
Pin = spectrumPower(am15);
IRV = 8910;

samples = { ...
  struct('tag','1um',  't_i',1000e-7,'N_i',220,'ref','Voc=1180mV Jsc=2.8mA/cm2 FF=50%% eta=1.6%% (Numune 1)'), ...
  struct('tag','0p3um','t_i',300e-7, 'N_i',80, 'ref','Voc=880mV  Jsc=3.8mA/cm2 FF=71%% eta=2.3%% (Numune 2)') ...
};

for si = 1:numel(samples)
  s = samples{si};
  fprintf('da Silva PIN %s -- gummelSolveQFP, accepted-gate + homotopy fallback\n', s.tag);
  dev = buildDaSilvaDev(s.t_i, s.N_i, IRV);
  N = dev.N;
  G = deviceGenAndTransmitML(dev, am15, @arcReflectanceGaP);
  % Vmax capped at a fixed, pre-bracketed value rather than the dynamic
  % max(1.4,-Vbi*1.05) formula: this device's Vbi=2.73V (undoped
  % intrinsic region, much larger depletion swing than the doped-
  % homojunction devices elsewhere in this project) would give
  % Vmax~2.87V, scanning ~1.3V past the real Voc crossing (bracketed at
  % V=1.30-1.50 for both samples via a quick direct-solve-only check)
  % into an all-rejected tail -- the same wasteful-Vmax mistake already
  % made once for #58 earlier this session. 1.55V comfortably covers
  % both samples' Voc with minimal wasted tail.
  Vmax = 1.55;
  Vs = 0:0.02:Vmax;
  prev = struct('psi',dev.psiEq,'phin',zeros(N,1),'phip',zeros(N,1));
  pts = [];
  nRejected = 0;
  tAll = tic;
  for k = 1:numel(Vs)
    V = Vs(k);
    r = gummelSolveQFP(dev, V, G, struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip));
    if ~r.isFinite, fprintf('NONFINITE at V=%.2f\n', V); break; end
    usedHomotopy = false;
    if ~qfpResultAccepted(r)
      rH = gummelSolveQFPHomotopy(dev, V, G, struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip));
      if rH.isFinite && rH.accepted
        r = rH; usedHomotopy = true;
      end
    end
    Jterm = mean(r.Jseg(N-20:N-1));
    if qfpResultAccepted(r)
      pts(end+1,:) = [V, Jterm]; %#ok<AGROW>
      if usedHomotopy, fprintf('  V=%.2f accepted via homotopy\n', V); end
    else
      nRejected = nRejected + 1;
      fprintf('  REJECTED V=%.2f\n', V);
    end
    prev = r;
  end
  fprintf('Accepted %d / %d points (%d rejected), total %.1f min\n', size(pts,1), numel(Vs), nRejected, toc(tAll)/60);
  m = pvMetrics(pts, Pin);
  fprintf('da Silva %s: Jsc=%.4f mA/cm2 Voc=%.1f mV FF=%.2f%% eta=%.2f%%\n', s.tag, m.Jsc*1e3, m.Voc*1e3, m.FF*100, m.PCE);
  fprintf('  Referans (real measured): %s\n\n', s.ref);

  T = array2table(pts,'VariableNames',{'V','J_Acm2'});
  writetable(T, fullfile(resultsDir, sprintf('dasilva_pin_%s_qfp_jv.csv', s.tag)));
  save(fullfile(resultsDir, sprintf('dasilva_pin_%s_qfp.mat', s.tag)), 'pts', 'm');
end
end

function dev = buildDaSilvaDev(t_i, N_i, Sint)
gap = gapParams();
Nprct = 3.0;
bac = gaasnpBandgap(Nprct/100);
opt = loadGaAsNPOptical(Nprct);
tauAbs = 0.6e-6; muP_Abs = 112;

layers(1) = struct('t',30e-7,'N',30,'Eg',gap.Eg,'chi',gap.chi,'eps',gap.eps, ...
  'Nd',1e19,'Na',0,'muN',gap.muN_doped,'muP',gap.muP_doped,'tau',gap.tau, ...
  'NcNv',gap.NcNv,'Brad',0,'alphaMax',0,'edgeWidth',0.05,'useOpticalTable',false,'opticalTable',[]);
layers(2) = struct('t',t_i,'N',N_i,'Eg',bac.Eg,'chi',4.0,'eps',12.9, ...
  'Nd',0,'Na',0,'muN',500,'muP',muP_Abs,'tau',tauAbs,'NcNv',1e19, ...
  'Brad',0,'alphaMax',0,'edgeWidth',0.05,'useOpticalTable',true,'opticalTable',opt);
% N=20 (was 60): thin (100nm), heavily-doped (Na=1e19) rear layers
% amplify the Scharfetter-Gummel flux sensitivity when meshed too finely
% (N.16.10) -- a controlled A/B/C test on #50's analogous BSF isolated
% mesh coarseness, not tau/chi/mobility, as the cause of the continuous
% V>=0.5 rejection "wall" this device hit at both Sint=0 and Sint=8910;
% N=20 (matching #58's own clean rear-layer mesh) fully clears it here
% (8/8 accepted at every tested V, both interface treatments).
layers(3) = struct('t',100e-7,'N',20,'Eg',gap.Eg,'chi',4.3,'eps',gap.eps, ...
  'Nd',0,'Na',1e19,'muN',gap.muN_doped,'muP',gap.muP_doped,'tau',gap.tau, ...
  'NcNv',gap.NcNv,'Brad',0,'alphaMax',0,'edgeWidth',0.05,'useOpticalTable',false,'opticalTable',[]);

junctions = struct('Sn_int',{Sint,Sint},'Sp_int',{Sint,Sint});
dev = buildDeviceML(layers, 2, junctions, struct('Sback',1e5));
dev.psiEq = solveEquilibrium(dev);
dev.alphaTable = opt;
end
