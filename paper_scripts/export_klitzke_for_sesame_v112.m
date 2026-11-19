function export_klitzke_for_sesame_v112()
%EXPORT_KLITZKE_FOR_SESAME_V112 (2026-11-13) As export_klitzke_for_sesame.m, but with the optics of the paper devices
% (AlInP window optical table, 85/45 nm MgF2/ZnS ARC; generation as run_klitzke_table_v112.m) and the RHJ of the main text only
% (N_D = 2e16, rear N_A = 2e17, 30/190/30 nodes = production device). Window 30 nodes also in the HJ. Tags *_v112.
% Original header: Export Klitzke-type GaInP homojunction (HJ) and rear-heterojunction (RHJ) devices under the digitized CW-LED SPD
% (100 and 1000 lx) for an independent Sesame cross-check (external_validation/sesame_klitzke_crosscheck.py).
% Same equations in both codes: ONE SRH centre per layer at E_i (tau_n, tau_p may differ), radiative B and Auger Cn, Cp where set, ideal heterointerfaces
% (Sn_int = Sp_int = 0), minority velocity Sback at the rear, ideal front contact. The multi-centre SRH of the HJ base of the paper is replaced by
% its third centre (tau_n = 39.5 ns, tau_p = 89.6 ns), because Sesame has one bulk SRH centre per material; these are therefore solver-check devices,
% not the calibrated devices of the validation table.
% Devices: KHJ (SRH only), KHJrad (same stack with B = 1e-10 in the GaInP layers and Auger 3e-30 everywhere), KRHJhi (N_D = 5e17, rear N_A = 1e19),
% KRHJlo (N_D = 2e16, rear N_A = 2e17); the RHJ stacks and parameters are those of run_klitzke2025_rhj_cwled_doping_sweep.m.
% Output: external_validation/sesame_input/<tag>_nodes.csv (x, layer, Eg, chi, eps, Nc, Nv, Nd, Na, muN, muP, tauN, tauP, Brad, Cn, Cp, G), <tag>_mlJV.csv, <tag>_meta.json.
assertCanonicalPath(); warning('off','all');
root = fileparts(fileparts(mfilename('fullpath'))); outDir = fullfile(root,'external_validation','sesame_input'); if ~isfolder(outDir), mkdir(outDir); end
T = readtable(resolveDataFile('Klitzke2025_CWLED_SPD_digitized.csv'));
specRaw = [T.lambda_nm,T.relative_irradiance];
devs = {'KHJ','KHJrad','KRHJlo'};
alOpt = loadAlInPOptical(); arcFn = @(lam) tmmReflectance(lam,[1,1.38,2.35,cidx(alOpt,lam)],[85,45]);
lux = [100 1000]; targetPin = [0.31 3.09];
Vg = (0:0.02:1.60)';
for id = 1:numel(devs)
  dev = buildDev(devs{id},alOpt);
  for il = 1:2
    tag = sprintf('%s_v112_%dlx',devs{id},lux(il));
    spec = specRaw; spec(:,2) = spec(:,2)*(targetPin(il)/spectrumPower(spec)); Pin = spectrumPower(spec);
    G = deviceGenAndTransmitML(dev,spec,arcFn);
    N = dev.N;
    Tn = table(dev.x(:), dev.layerIdxOfNode(:), dev.Eg(:), dev.chi(:), dev.eps(:), dev.Nc(:), dev.Nv(:), dev.Nd(:), dev.Na(:), dev.muN(:), dev.muP(:), ...
      dev.srhTauN(:,1), dev.srhTauP(:,1), dev.Brad(:), dev.CnAuger(:), dev.CpAuger(:), G(:), ...
      'VariableNames',{'x_cm','layer','Eg','chi','eps','Nc','Nv','Nd','Na','muN','muP','tauN','tauP','Brad','Cn','Cp','G'});
    writetable(Tn,fullfile(outDir,sprintf('%s_nodes.csv',tag)));
    [cv,m] = indoorCurve(dev,G,Vg,Pin);
    writetable(cv,fullfile(outDir,sprintf('%s_mlJV.csv',tag)));
    meta = struct('code',devs{id},'tag',tag,'N',N,'Sback',dev.Sback,'Pin_Wm2',Pin,'Jsc_mAcm2',m.Jsc*1e3,'Voc_mV',m.Voc*1e3,'FF_pct',m.FF*100,'PCE_pct',m.PCE);
    fid = fopen(fullfile(outDir,sprintf('%s_meta.json',tag)),'w'); fprintf(fid,'%s',jsonencode(meta)); fclose(fid);
    fprintf('%s: Jsc=%.4g uA/cm2 Voc=%.1f mV FF=%.2f PCE=%.3f (N=%d, accepted %d/%d)\n',tag,m.Jsc*1e6,m.Voc*1e3,m.FF*100,m.PCE,N,sum(cv.Accepted),height(cv));
  end
end
end

function dev = buildDev(name,alOpt)
gp = struct('Eg',1.90,'chi',4.10,'eps',11.8,'Nc',6.5e17,'Nv',1.5e19,'muN',500,'muP',30,'EgB',2.35,'chiB',3.78,'NcB',2.5e18,'NvB',7e18,'muNB',100,'muPB',10);
ctr3 = struct('tauN',39.5e-9,'tauP',89.6e-9,'EtRelEi_eV',0);
switch name
  case {'KHJ','KHJrad'}
    rad = strcmp(name,'KHJrad'); B = 1e-10*rad; CA = 3e-30*rad;
    L(1) = lay(25e-7,30,gp.EgB,gp.chiB,gp.NcB,gp.NvB,8e18,0,gp.muNB,gp.muPB,1e-9,0,CA,[]);
    L(2) = lay(150e-7,60,gp.Eg,gp.chi,gp.Nc,gp.Nv,2e18,0,gp.muN,gp.muP,248e-9,B,CA,struct('tauN',337e-9,'tauP',248e-9,'EtRelEi_eV',0));
    L(3) = lay(700e-7,160,gp.Eg,gp.chi,gp.Nc,gp.Nv,0,2e16,gp.muN,gp.muP,248e-9,B,CA,ctr3);
    jn = struct('Sn_int',{0,0},'Sp_int',{0,0});
  case {'KRHJhi','KRHJlo'}
    if strcmp(name,'KRHJhi'), ND = 5e17; NAr = 1e19; else, ND = 2e16; NAr = 2e17; end
    L(1) = lay(25e-7,30,gp.EgB,gp.chiB,gp.NcB,gp.NvB,8e18,0,gp.muNB,gp.muPB,1e-9,0,3e-30,[]);
    L(2) = lay(850e-7,190,gp.Eg,gp.chi,gp.Nc,gp.Nv,ND,0,gp.muN,gp.muP,248e-9,1e-10,3e-30,struct('tauN',337e-9,'tauP',248e-9,'EtRelEi_eV',0));
    L(3) = lay(100e-7,30,gp.EgB,gp.chiB,gp.NcB,gp.NvB,0,NAr,gp.muNB,gp.muPB,1e-9,0,3e-30,[]);
    jn = struct('Sn_int',{0,0},'Sp_int',{0,0});
end
L(1).useOpticalTable = true; L(1).opticalTable = alOpt;
dev = buildDeviceML(L,2,jn,struct('Sback',1e5));
dev.psiEq = solveEquilibrium(dev);
end

function L = lay(t,N,Eg,chi,Nc,Nv,Nd,Na,muN,muP,tau,B,CA,ctr)
L = struct('t',t,'N',N,'Eg',Eg,'chi',chi,'eps',11.8,'Nd',Nd,'Na',Na,'muN',muN,'muP',muP,'tau',tau,'Nc',Nc,'Nv',Nv,'NcNv',sqrt(Nc*Nv), ...
  'Brad',B,'CnAuger',CA,'CpAuger',CA,'alphaMax',3e4,'edgeWidth',0.05,'useOpticalTable',false,'opticalTable',[],'srhCenters',ctr);
if Eg > 2.2, L.alphaMax = 0; end    % AlInP layers (window/rear) are transparent to the CW-LED in the sweep scripts (alphaMax = 0)
end

function [cv,m] = indoorCurve(dev,G,Vg,Pin)
N = dev.N; idx = find(abs(Vg-0.02)<1e-12,1);
r02 = gummelSolveQFPAnalytic(dev,0.02,G,struct('psiInit',dev.psiEq,'phinInit',zeros(N,1),'phipInit',zeros(N,1))); prev = r02;
for v = [0.01 0.005 0], r0 = gummelSolveQFPAnalytic(dev,v,G,struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip)); prev = r0; end
st = repmat(r02,numel(Vg),1); st(1) = r0; st(idx) = r02; prev = r02;
for k = idx+1:numel(Vg)
  rk = gummelSolveQFPAnalytic(dev,Vg(k),G,struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip)); st(k) = rk; if rk.isFinite, prev = rk; end
end
J = nan(numel(Vg),1); acc = false(numel(Vg),1);
for k = 1:numel(Vg), if st(k).isFinite, J(k) = mean(st(k).Jseg(N-20:N-1)); acc(k) = qfpResultAccepted(st(k)); end, end
cv = table(Vg,J,double(acc),'VariableNames',{'V_V','J_Acm2','Accepted'});
m = pvMetrics([Vg(acc),J(acc)],Pin);
end

function z=cidx(opt,lam)
E=1239.842/lam; n=interp1(opt.E,opt.n,E,'linear','extrap'); a=interp1(opt.E,opt.alpha,E,'linear',0); k=max(a,0)*(lam*1e-7)/(4*pi); z=n-1i*k;
end
