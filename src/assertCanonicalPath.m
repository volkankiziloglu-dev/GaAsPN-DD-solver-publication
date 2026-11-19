function assertCanonicalPath()
% Dis inceleme madde 4.7: MATLAB path onceligi (startup.m'nin genpath ile
% TUM matlab/ agacini eklemesi, RearHJT/'deki AYNI-ISIMLI ama eps0-hatali
% eski solver dosyalariyla CAKISMA riski) alfabetik siralamaya
% (GaAsNP_RHJ_final < RearHJT) dayaniyordu -- calisiyor ama kirilgan bir
% garanti. Bu fonksiyon, ana sonuc uretecek her script'in basinda
% cagirilarak DOGRU (bu klasordeki) cozucu dosyalarinin gercekten
% kullanildigini ACIKCA dogrular -- yanlislikla eski RearHJT/ dosyasi
% yuklenirse SESSIZCE yanlis sonuc uretmek yerine hemen, anlasilir bir
% hatayla durur.
% DUZELTME (dis inceleme madde 1.1): startsWith(p,here) TEORIK olarak
% "here" ile ayni ONEKE sahip farkli bir klasorle (orn.
% GaAsNP_RHJ_final2/) yanlislikla eslesebilirdi. Klasor SINIRINI da
% iceren (filesep ile biten) beklenen kok kullaniliyor.
here = fileparts(mfilename('fullpath'));
expectedRoot = [here filesep];
fns = {'newtonPoisson','newtonPoissonML','gummelSolve','gummelSolveML', ...
       'sweepJV','sweepJVML','setupDevice','setupDeviceML','thomasSolve', ...
       'gummelSolveQFPAnalytic','gummelSolveQFP','assembleResidualQFP','interfaceRecombination', ...   % dd_solver v1.1.1: the QFP production path
       'runAcceptedJVSweep','indoorAcceptanceGate','pvMetrics','buildDeviceML', ...
       'sgBernoulli','qfpResultAccepted'};
% Publication package: legacy solver files (gummelSolve, sweepJV, setupDevice) are not part of the
% dependency closure and are therefore not shipped; only the functions present in src/ are checked.
fns = fns(cellfun(@(f) isfile(fullfile(here,[f '.m'])), fns));
for i = 1:numel(fns)
  p = which(fns{i});
  if isempty(p)
    error('assertCanonicalPath: %s hic bulunamadi (path''e eklenmemis olabilir)', fns{i});
  end
  if ~startsWith(p, expectedRoot)
    error(['assertCanonicalPath: %s BEKLENMEDIK yerden yukleniyor:\n  %s\n' ...
      'Beklenen klasor: %s\n' ...
      'Bu, eski (eps0-hatali) RearHJT/ kopyasinin path onceligi kazandigini gosterir -- ' ...
      'sonuclar SESSIZCE yanlis olabilir. addpath sirasini kontrol edin.'], ...
      fns{i}, p, here);
  end
end
end
