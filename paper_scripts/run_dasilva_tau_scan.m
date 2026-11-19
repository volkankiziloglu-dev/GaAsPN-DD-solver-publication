function T = run_dasilva_tau_scan()
%RUN_DASILVA_TAU_SCAN Can a shorter absorber lifetime alone explain the low J_sc of the measured da Silva GaAsPN p-i-n cells (2.8 and 3.8 mA/cm2)? (referee check R2-5)
% The p-i-n stack of run_dasilva_pin_qfp.m is solved at V = 0 under AM1.5G for i-region lifetimes from 1e-11 to 6e-7 s (single midgap centre, nominal mobilities).
% Output: results/paper_calibrated_regime/dasilva_tau_scan.csv
assertCanonicalPath(); warning('off','all');
root = fileparts(fileparts(mfilename('fullpath'))); outDir = fullfile(root,'results','paper_calibrated_regime'); if ~isfolder(outDir), mkdir(outDir); end
am15 = loadRealAM15G(); Pin = spectrumPower(am15);
samples = {struct('tag','1um','t_i',1000e-7,'N_i',220,'Jmeas',2.8), struct('tag','0p3um','t_i',300e-7,'N_i',80,'Jmeas',3.8)};
rows = {};
for si = 1:2
  s = samples{si};
  for tau = [1e-11 1e-10 1e-9 1e-8 1e-7 6e-7]
    dev = buildDaSilvaDev(s.t_i,s.N_i,8910,tau); N = dev.N;
    G = deviceGenAndTransmitML(dev,am15,@arcReflectanceGaP);
    r = gummelSolveQFPAnalytic(dev,0.02,G,struct('psiInit',dev.psiEq,'phinInit',zeros(N,1),'phipInit',zeros(N,1))); prev = r;
    for v = [0.01 0.005 0], r0 = gummelSolveQFPAnalytic(dev,v,G,struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip)); prev = r0; end
    Jsc = mean(r0.Jseg(N-20:N-1))*1e3; Jgen = 1.602176634e-19*trapz(dev.x,G)*1e3;
    rows(end+1,:) = {s.tag, tau, Jsc, s.Jmeas, Jgen, r0.isFinite}; %#ok<AGROW>
    disp([s.tag ' tau=' num2str(tau) ' Jsc=' num2str(Jsc,'%.3f') ' mA/cm2 (generated ' num2str(Jgen,'%.2f') ')']);
  end
end
T = cell2table(rows,'VariableNames',{'Sample','tau_s','Jsc_model_mAcm2','Jsc_measured_mAcm2','Jgen_mAcm2','Finite'});
writetable(T,fullfile(outDir,'dasilva_tau_scan.csv'));
end

function dev = buildDaSilvaDev(t_i, N_i, Sint, tauAbs)
gap = gapParams();
Nprct = 3.0;
bac = gaasnpBandgap(Nprct/100);
opt = loadGaAsNPOptical(Nprct);
muP_Abs = 112;

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
