function out = run_sukrittanon_ganp_validation(mode, varargin)
%RUN_SUKRITTANON_GANP_VALIDATION Sukrittanon 2015 GaP / GaNP p-i-n cells.
%
% Experimental reference: S. Sukrittanon, PhD thesis, UC San Diego (2015),
% Chapter 4 (published as Sukrittanon et al., APL 107, 153901 (2015),
% doi 10.1063/1.4933317) -- GaP-800, GaNP-800,
% GaNP-2000 on p-GaP(001), AM1.5G, with/without Si3N4(50)/SiO2(20) AR.
%
% Layers (thesis Fig. 4.1/Sec. 4.2): n+ GaP 15 nm (8e18) / n GaP 100 nm
% (4e18) / i-(Ga)(N)P t_i / p-GaP substrate (6.5e17, 225 um; truncated
% here to a few diffusion lengths).  No window layer.  [N]~1.8 %, Eg 2.05 eV.
%
% mode:
%   'jsc'    quick optical/collection check: Jsc and EQE, no J-V sweep
%   'jv'     full AM1.5G J-V sweep (args: tauI, Sint)
%
% Reference values (thesis Table 4.1, with AR):
%   GaP-800    Jsc 3.17 Voc 1.19 FF 76 eta 2.9
%   GaNP-800   Jsc 6.08 Voc 1.26 FF 71 eta 5.5
%   GaNP-2000  Jsc 8.53 Voc 1.33 FF 69 eta 7.9
%   (without AR: 2.69/1.19/77/2.5, 5.52/1.26/73/5.1, 7.80/1.31/71/7.3)
assertCanonicalPath();
warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:singularMatrix');
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'src'),'-begin');
addpath(fullfile(root,'paper_scripts'),'-begin');
% (publication package: the project tests/ folder is not needed)
if nargin < 1, mode = 'jsc'; end
p = inputParser;
p.addParameter('sample','GaNP-2000');
p.addParameter('tauI',30e-9);
p.addParameter('Sint',0);
p.addParameter('shift',0.10);       % eV, 2.1%N table -> 1.8%N
p.addParameter('muScale',1.0);      % i-layer mobility scale
p.addParameter('NaI',0);            % background p-type doping in i-layer
p.addParameter('AR',true);
p.addParameter('Rs',0);             % ohm cm2, external series resistance (post-processed)
p.addParameter('tauBase',NaN);      % s, p-GaP base SRH lifetime (NaN: thesis Lb)
p.addParameter('tauEmit',NaN);      % s, n-GaP emitter SRH lifetime (NaN: thesis Le)
p.addParameter('tauGaPi',NaN);      % s, i-GaP lifetime in GaP-800 control (NaN: = tauI)
p.addParameter('tag','');
p.parse(varargin{:});
o = p.Results;

s = sampleSpec(o.sample);
dev = buildSukrittanonDev(s, o);
N = dev.N;
am15 = loadRealAM15G();
Pin = spectrumPower(am15);
spec = am15(am15(:,1)>=280 & am15(:,1)<=800,:);
if o.AR, frontR = @arcSukrittanon; else, frontR = @arcBareGaP; end

switch lower(mode)
  case 'jsc'
    out = jscAndEqe(dev, s, spec, frontR, Pin);
  case 'dev'
    out = struct('dev',dev,'spec',{spec},'frontR',{frontR},'Pin',Pin);
  case 'dark'
    out = darkSweep(dev, s, o);
  case 'jv'
    out = jvSweep(dev, s, spec, frontR, Pin, o);
  otherwise
    error('unknown mode');
end
end

function s = sampleSpec(name)
switch name
  case 'GaP-800',    s = struct('name',name,'ti',800e-7, 'ganp',false,'N',60);
  case 'GaNP-800',   s = struct('name',name,'ti',800e-7, 'ganp',true, 'N',60);
  case 'GaNP-2000',  s = struct('name',name,'ti',2000e-7,'ganp',true, 'N',100);
  otherwise, error('unknown sample');
end
end

function dev = buildSukrittanonDev(s, o)
gap = gapParams();
vt = 0.025852;
optGaP = loadGaPOptical();
optI = optGaP;
EgI = gap.Eg; chiI = gap.chi;
if s.ganp
  pn = loadGaPNOptical();
  optI = pn;
  optI.E = pn.E + o.shift;            % shift the 2.1 %N edge to 1.8 %N
  optI.Emin = min(optI.E); optI.Emax = max(optI.E);
  EgI = 2.05;
  chiI = gap.chi + (gap.Eg - EgI);    % Delta Ev ~ 0, all offset in Ec
end
% thesis Table 2.1: emitter Le=190 nm, De=1.68 cm2/s; base Lb=4.1 um, Db=1.53
tauE = (190e-7)^2/1.68;     % ~2.1e-10 s
tauB = (4.1e-4)^2/1.53;     % ~1.1e-7 s
muPe = 1.68/vt;             % minority holes in n-GaP emitter
muNb = 1.53/vt;             % minority electrons in p-GaP base
if ~isnan(o.tauBase), tauB = o.tauBase; end
if ~isnan(o.tauEmit), tauE = o.tauEmit; end
Lb = 25e-4;                 % truncated base, cm (~6 Lb)

base = struct('chi',gap.chi,'eps',gap.eps,'NcNv',gap.NcNv,'Brad',0, ...
  'alphaMax',0,'edgeWidth',0.05,'useOpticalTable',true,'opticalTable',optGaP);
L1 = base; L1.t=15e-7;  L1.N=10; L1.Eg=gap.Eg; L1.Nd=8e18; L1.Na=0; L1.muN=gap.muN_doped; L1.muP=muPe; L1.tau=tauE;
L2 = base; L2.t=100e-7; L2.N=20; L2.Eg=gap.Eg; L2.Nd=4e18; L2.Na=0; L2.muN=gap.muN_doped; L2.muP=muPe; L2.tau=tauE;
L3 = base; L3.t=s.ti;   L3.N=s.N; L3.Eg=EgI; L3.chi=chiI; L3.Nd=0; L3.Na=o.NaI;
L3.muN=100*o.muScale; L3.muP=80*o.muScale; L3.tau=o.tauI; if ~s.ganp && ~isnan(o.tauGaPi), L3.tau=o.tauGaPi; end
L3.opticalTable=optI;
L3.NcNv = 1e19;
L4 = base; L4.t=Lb;     L4.N=80; L4.Eg=gap.Eg; L4.Nd=0; L4.Na=6.5e17; L4.muN=muNb; L4.muP=gap.muP_doped; L4.tau=tauB;
layers = [L1 L2 L3 L4];
if s.ganp, Sint = o.Sint; else, Sint = 0; end
junctions = struct('Sn_int',{0,Sint,Sint},'Sp_int',{0,Sint,Sint});
dev = buildDeviceML(layers, 3, junctions, struct('Sback',1e6));
dev.psiEq = solveEquilibrium(dev);
dev.alphaTable = optI;
end

function R = arcSukrittanon(lam)
% Si3N4 50 nm / SiO2 20 nm (thesis Sec. 4.5), on GaP.
nSiN = 1.95 + 0.02/(lam*1e-3)^2;
R = tmmReflectance(lam, [1.0, 1.46, nSiN, gapRefractiveIndex(lam)], [20, 50]);
end

function R = arcBareGaP(lam)
n = gapRefractiveIndex(lam);
R = ((n-1)/(n+1))^2;
end

function out = jscAndEqe(dev, s, spec, frontR, Pin)
c = constants();
N = dev.N;
prev = struct('psi',dev.psiEq,'phin',zeros(N,1),'phip',zeros(N,1));
% broadband Jsc
G = deviceGenAndTransmitML(dev, spec, frontR);
r = gummelSolveQFPAnalytic(dev, 0, G, struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip));
Jsc = mean(r.Jseg(N-20:N-1));
fprintf('%s  Jsc(V=0) = %.3f mA/cm2  accepted=%d\n', s.name, Jsc*1e3, qfpResultAccepted(r));
% EQE
lam = (300:10:700)';
eqe = zeros(size(lam));
for k = 1:numel(lam)
  irr = 1; % W/m2/nm, small enough; linear regime
  Gk = singleWavelengthGenML(dev, lam(k), irr*0.01, frontR); % 1e-2 W/m2/nm
  rk = gummelSolveQFPAnalytic(dev, 0, Gk, struct('psiInit',dev.psiEq,'phinInit',zeros(N,1),'phipInit',zeros(N,1)));
  if ~qfpResultAccepted(rk)
    rk = gummelSolveQFPHomotopy(dev, 0, Gk, struct('innerSolver',@gummelSolveQFPAnalytic));
  end
  if qfpResultAccepted(rk), Jk = mean(rk.Jseg(N-20:N-1)); else, Jk = NaN; end
  phi = photonFlux(lam(k), irr*0.01);
  eqe(k) = Jk/(c.Q*phi);
end
T = table(lam, eqe, 'VariableNames', {'lambda_nm','EQE'});
disp(T(1:2:end,:));
out = struct('Jsc',Jsc,'lam',lam,'eqe',eqe);
end

function out = jvSweep(dev, s, spec, frontR, Pin, o)
N = dev.N;
G = deviceGenAndTransmitML(dev, spec, frontR);
Vs = 0:0.02:1.55;
prev = struct('psi',dev.psiEq,'phin',zeros(N,1),'phip',zeros(N,1));
pts = []; nRej = 0; rejV = [];
hopts = @(pv) struct('innerSolver',@gummelSolveQFPAnalytic,'psiInit',pv.psi,'phinInit',pv.phin,'phipInit',pv.phip);
for k = 1:numel(Vs)
  V = Vs(k);
  r = gummelSolveQFPAnalytic(dev, V, G, struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip));
  if ~r.isFinite, break; end
  if ~qfpResultAccepted(r)
    rH = gummelSolveQFPHomotopy(dev, V, G, hopts(prev));
    if rH.isFinite && rH.accepted, r = rH; end
  end
  if qfpResultAccepted(r)
    pts(end+1,:) = [V, mean(r.Jseg(N-20:N-1))]; %#ok<AGROW>
    prev = r;
    if size(pts,1)>3 && pts(end,2) < -2e-3, break; end
  else
    nRej = nRej+1; rejV(end+1) = V; %#ok<AGROW>
  end
end
ptsInt = pts;
if o.Rs > 0, pts = [pts(:,1) - pts(:,2)*o.Rs, pts(:,2)]; end
m = pvMetrics(pts, Pin);
fprintf('%s tau=%.1e ns Sint=%g shift=%.2f muScale=%.2f NaI=%g AR=%d Rs=%g: Jsc=%.3f Voc=%.1f mV FF=%.1f%% eta=%.2f%%  (rej=%d %s)\n', ...
  s.name, o.tauI*1e9, o.Sint, o.shift, o.muScale, o.NaI, o.AR, o.Rs, m.Jsc*1e3, m.Voc*1e3, m.FF*100, m.PCE, nRej, mat2str(rejV(1:min(end,3))));
resultsDir = fullfile(fileparts(fileparts(mfilename('fullpath'))),'results','sukrittanon_ganp');
if ~exist(resultsDir,'dir'), mkdir(resultsDir); end
fn = sprintf('%s_tau%.0e_S%g_sh%.2f_mu%.2f_Na%g_AR%d_tB%.0e_tE%.0e_tG%.0e%s', s.name, o.tauI, o.Sint, o.shift, o.muScale, o.NaI, o.AR, o.tauBase, o.tauEmit, o.tauGaPi, o.tag);
writetable(array2table(ptsInt,'VariableNames',{'V','J_Acm2'}), fullfile(resultsDir, [fn '.csv']));
out = struct('m',m,'pts',pts,'ptsInt',ptsInt);
end

function out = darkSweep(dev, s, o)
% Dark J-V (intrinsic, no Rs): warm-started continuation from V=0.
N = dev.N; G = zeros(N,1);
Vs = 0.4:0.05:1.45;
prev = struct('psi',dev.psiEq,'phin',zeros(N,1),'phip',zeros(N,1));
pts = zeros(0,3);
for k = 1:numel(Vs)
  r = gummelSolveQFPAnalytic(dev, Vs(k), G, struct('psiInit',prev.psi,'phinInit',prev.phin,'phipInit',prev.phip));
  if ~r.isFinite, break; end
  Jv = mean(r.Jseg(N-20:N-1));
  pts(end+1,:) = [Vs(k), -Jv, double(qfpResultAccepted(r))]; %#ok<AGROW>
  prev = r;
end
fprintf('%s dark tau=%.1e: ', s.name, o.tauI); fprintf('[%.2f %.2e %d] ', pts'); fprintf(newline);
resultsDir = fullfile(fileparts(fileparts(mfilename('fullpath'))),'results','sukrittanon_ganp');
writetable(array2table(pts,'VariableNames',{'V','Jdark_Acm2','accepted'}), fullfile(resultsDir, sprintf('dark_%s_tau%.0e.csv', s.name, o.tauI)));
out = pts;
end
