function n = gapRefractiveIndex(lambda_nm)
% GaP kirilma indisi n(lambda), GORECE KABA/YAKLASIK -- Aspnes & Studna
% (1983) turu olcumlerin bilinen genel egilimine (indirect-gap III-V,
% 550nm/2.26eV bant kenarina yaklastikca artan dispersiyon, uzun dalga
% boyunda ~3.2'ye duzlesme) DAYALI, ama TAM bir kaynak tablosuna
% birebir uydurulmus DEGIL -- RearHJT/gapParams.m'deki chi belirsizligiyle
% AYNI ruhta, "makul baslangic noktasi" olarak ele alinmali. Sadece ARC
% yansima hesaplamasi (Gorev #10) icin kullanilir; GaP hala alphaMax=0
% (optik olarak seffaf) varsayimini KORUR, bu sadece n (yansima icin) ekler.
lamTab = [400 450 500 550 600 650 700 750 800 850 900];
nTab   = [3.90 3.75 3.60 3.55 3.50 3.42 3.35 3.29 3.24 3.20 3.16];
n = interp1(lamTab, nTab, lambda_nm, 'linear', 'extrap');
end
