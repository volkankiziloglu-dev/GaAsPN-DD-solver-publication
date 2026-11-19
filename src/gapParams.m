function p = gapParams()
% GaP (galyum fosfit) malzeme parametreleri -- 4 katmanli hedef yapinin
% n+ pencere ve p+ arka emiter katmanlari icin (bkz. Gorev #9).
%
% ONEMLI: a-Si:H/c-Si kalibrasyonunda (Si/run_asi_csi_hjt_validate.m,
% MATHEMATICAL_MODEL.md Sec. I.2) ogrenilen ders burada da gecerli --
% literatur degerleri (ozellikle elektron ilgisi chi ve agir katkilamada
% mobilite) referans kaynaga gore ONEMLI OLCUDE degisebilir. Asagida her
% alan LITERATUR (dogrudan olcum/tablo degeri, kaynakli) veya TAHMINI/
% KALIBRASYON-GEREKLI (net bir kaynak yok, deneme-yanilma/ampirik ayar
% gerekebilir) olarak ACIKCA etiketlenmistir. Task #9'daki ilk calistirma
% bu yuzden bir DOGRULAMA degil, bir TUTARLILIK/KESIF calistirmasi olarak
% ele alinmalidir (literaturde bu TAM 4-katmanli yapi icin karsilastirma
% verisi yok).

% --- LITERATUR (Ioffe NSM Archive, ioffe.ru/SVA/NSM/Semicond/GaP/,
%     Levinshtein/Rumyantsev/Shur "Handbook Series on Semiconductor
%     Parameters" temelli derleme; 300 K) ---
p.Eg = 2.26;              % eV, temel (indirect, X-vadisi) bant araligi, 300K.
                           % NOT: direkt (Gamma-Gamma) gecis 2.78 eV'de --
                           % optik soğurma kenari icin (ince/seffaf pencere
                           % katmani oldugundan bu projede alphaMax=0
                           % varsayilarak soğurma ihmal ediliyor, bkz. asagi).
p.EgT_a = 2.34;            % Eg(T) = EgT_a - EgT_b*T^2/(T+EgT_c) [eV], 0-1200K
p.EgT_b = 6.0e-4;
p.EgT_c = 460;             % K
p.eps = 11.1;              % statik dielektrik sabiti (bazi kaynaklarda 9-11.1
                           % araliginda; yuksek-frekans eps_inf=9.11 ile
                           % karistirilmamali -- burada statik/DC deger kullanildi,
                           % Poisson denklemi icin dogru olan budur).
p.NcNv = 1.85e19;          % cm^-3, sqrt(Nc*Nv), Nc=1.8e19 Nv=1.9e19 (300K,
                           % Ioffe). Bu proje kongansiyonuyla (bkz. buildDevice.m
                           % NcNv_ETL/NcNv_Abs) tek deger olarak kullanilir.
p.muN_max = 250;           % cm^2/Vs, KATKISIZ/hafif katkili ust sinir (300K).
p.muP_max = 150;           % cm^2/Vs, ayni sekilde katkisiz ust sinir.

% --- LITERATUR (agir katkida mobilite, GUNCELLENDI -- Gorev "parametre
%     avi": Ma & Lundstrom tipi ab-initio iyonize-safsizlik sacilmasi
%     hesabi, arXiv:2301.02323, olculmus degerlerle capraz-dogrulanmis) ---
p.muN_doped = 100;         % cm^2/Vs, Nd~1e19 cm^-3 icin. Kaynak: 1e18'de
                           % hesap=157, olcum~100 cm2/Vs (hesap sistematik
                           % olarak olcumden yuksek); 1e19'da acik elektron
                           % degeri verilmedi ama holler icin 1e18->1e19
                           % gecisinde ~%24 dusus goruluyor (124->94) --
                           % elektron icin benzer oranli dusus + olcum-hesap
                           % farkini goz onune alarak 100 cm2/Vs SECILDI
                           % (eski deger 20 cm2/Vs idi -- 5 KAT DUSUK
                           % tahmin edilmisti).
p.muP_doped = 80;          % cm^2/Vs, Na~1e19 cm^-3 icin. Kaynak DOGRUDAN
                           % veriyor: hesap=94 cm2/Vs @ 1e19 (1e18'de
                           % hesap=124, olcum=90 -- yine hesap>olcum egilimi);
                           % bu farki goz onune alarak 80 cm2/Vs SECILDI
                           % (eski deger 10 cm2/Vs idi -- 8-9 KAT DUSUK
                           % tahmin edilmisti). Bu, RearHJT/run_4layer_*.m
                           % script'lerinde HENUZ uygulanmadi -- cagiran kod
                           % gap.muN_doped/gap.muP_doped'i otomatik cekiyor,
                           % bu yuzden bir sonraki calistirmada KENDILIGINDEN
                           % guncellenecek.

% --- TAHMINI / KALIBRASYON-GEREKLI (net, tek bir kaynak degeri YOK) ---
p.chi = 3.8;               % eV, elektron ilgisi -- BUYUK BELIRSIZLIK: farkli
                           % kaynaklar 3.8-4.3 eV araliginda degerler veriyor
                           % (klasik Sze-tipi tablo ~4.3 eV; bazi guncel
                           % derlemeler ~3.8 eV bildiriyor). Anderson kurali
                           % zaten yaklasik oldugundan (bkz. buildDevice.m /
                           % buildDeviceML.m yorumlari, MATHEMATICAL_MODEL.md
                           % Sec. A.2), a-Si:H/c-Si kalibrasyonunda oldugu gibi
                           % bu deger muhtemelen deneysel Voc/band-offset
                           % davranisina gore AYARLANACAK bir baslangic noktasi
                           % olarak ele alinmali, kesin fiziksel sabit degil.
                           %
                           % GOREV #11 BULGUSU (bkz. RearHJT/gummelSolveML.m
                           % arastirmasi): 4.3 (deltaEc=+0.3eV, n+GaP pencere/
                           % n-GaAsNP absorber arayuzunde elektronlar icin
                           % ENGELLEYICI bir sivri uc) ile bu modelin dejenere-
                           % OLMAYAN Boltzmann istatistigi (n=ni*exp(psiEff/VT),
                           % ustten sinirsiz) altinda ileri-yon onyaginda yerel
                           % elektron yigilmasi SINIRSIZ buyuyup Gummel dongusunu
                           % kararsizlastiriyordu (dogrulandi: n~5.8e23 cm^-3,
                           % NcNv=1.85e19'un ~30000 kati). 3.8 (deltaEc=-0.2eV,
                           % ENGELLEYICI DEGIL) ile ayni yapi TAM J(V) egrisini
                           % sorunsuz veriyor. Bu yuzden burada varsayilan olarak
                           % 3.8 SECILDI -- hem literaturde atifli, hem de bu
                           % mimarideki (n-n izotip on arayuz) elektron tasinimi
                           % icin fiziksel/elektriksel olarak TERCIH EDILEN yon.
                           % Gercekten engelleyici bir arayuz gerekirse (bazi
                           % gercek cihazlarda kasitli spike/notch tasarimlari
                           % olabilir), bu modelin dejenere-olmayan istatistigi
                           % o rejimde GUVENILMEZ olacaktir -- bkz. Gorev #11.
                           %
                           % NOT: bu "engelleyici" sorunu SADECE cogunluk-
                           % tasiyici (izotip, n-n veya p-p) arayuzlerde
                           % kritiktir. AZINLIK-tasiyici bir arayuzde (orn.
                           % gercek p-n eklemi olan absorber/arka-emiter
                           % sinirinda) ayni yonde bir engel BSF-benzeri
                           % FAYDALI bir etkidir ve ayni kararsizligi
                           % tetiklemez (dogrulandi, bkz. RearHJT/
                           % run_4layer_explore.m: pencere icin bu p.chi=3.8
                           % kullanilirken, arka emiter icin BILINCLI olarak
                           % chi=4.3 kullanilmasi tam J(V) egrisi veriyor).
                           % Yani p.chi TEK bir "GaP degeri" degil, konuma/
                           % role gore SECILMESI gereken bir baslangic noktasi
                           % olarak dusunulmeli.
p.tau = 1e-9;              % s, SRH omru -- TAHMIN, ince (20-150nm) pencere/emiter
                           % katmaninda tasiyicilarin buyuk kismi zaten toplanmadan
                           % once katmani gecer, bu deger hassas degil (dusuk
                           % duyarlilik beklenir, ama dogrulanmadi).

% --- OPTIK (bu projede pencere/emiter TRANSPARAN varsayiliyor) ---
p.alphaMax = 0;            % cm^-1, varsayilan: ince GaP pencere/emiter katmani
                           % OPTIK OLARAK SEFFAF kabul edilir (alphaMax=0),
                           % common/buildDevice.m'deki alphaETLmax=0 varsayilanina
                           % benzer sekilde -- gercek GaP dogrudan-olmayan bant
                           % gecisi nedeniyle zaten zayif soguruyor (bkz. Eg
                           % notundaki 2.78 eV direkt gecis), ama tam bir n(E),k(E)
                           % tablosu (RearHJT/data/ icin planlanan) olmadan bu
                           % basitlestirme kullaniliyor. Gorev #9'da net bir
                           % ihtiyac ortaya cikarsa RearHJT/data/ altina gercek
                           % optik tablo eklenmelidir.
end
