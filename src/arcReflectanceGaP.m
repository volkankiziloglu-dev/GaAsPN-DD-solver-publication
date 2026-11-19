function R = arcReflectanceGaP(lambda_nm)
% Cift katmanli ARC (Gorev #10): MgF2 (disarida, hava tarafi) / ZnS
% (icerde, GaP pencere tarafi) -- III-V (GaAs) gunes pilleri icin
% literaturde standart bir tasarim (kaynak: MgF2/ZnS cift-katman ARC,
% GaAs tek-eklem hucreler icin 450-1730nm araliginda R<%10 verdigi
% bildirilen d_MgF2=105nm / d_ZnS=89nm tasarimi -- bkz. arastirma notlari).
% GaP substratinin kirilma indisi de benzer mertebede (n~3.2-3.6) oldugundan
% bu tasarim GaP-onlu yapimiz icin de MAKUL bir baslangic noktasi, ama GaP'a
% OZEL olarak optimize edilmedi -- ileri bir adim olarak (Gorev #10'un
% "sonraki adimi") kalinliklar GaP'in kendi indisine gore yeniden
% optimize edilebilir.
n_air = 1.0;
n_MgF2 = 1.38; d_MgF2 = 105; % nm, disarida (hava tarafi)
n_ZnS = 2.35;  d_ZnS = 89;   % nm, icerde (GaP tarafi)
n_GaP = gapRefractiveIndex(lambda_nm);
R = tmmReflectance(lambda_nm, [n_air, n_MgF2, n_ZnS, n_GaP], [d_MgF2, d_ZnS]);
end
