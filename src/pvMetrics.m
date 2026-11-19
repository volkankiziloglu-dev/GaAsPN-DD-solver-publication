function m = pvMetrics(pts, Pin_Wm2)
% pts: Nx2 [V, J(A/cm2)], illuminated sweep, sign convention J(V=0)=Jsc>0 decreasing to 0 at Voc.
V = pts(:,1); J = pts(:,2);

Jsc = NaN;
for k = 1:numel(V)-1
  if V(k) <= 0 && V(k+1) >= 0
    t = (0-V(k))/(V(k+1)-V(k));
    Jsc = J(k) + t*(J(k+1)-J(k));
    break;
  end
end

Voc = NaN;
for k = 1:numel(V)-1
  if V(k) >= 0 && J(k) >= 0 && J(k+1) < 0
    t = (0-J(k))/(J(k+1)-J(k));
    Voc = V(k) + t*(V(k+1)-V(k));
    break;
  end
end

% DUZELTME (dis inceleme madde 2.1, KRITIK): ONCEDEN Voc/Jsc bulunamasa
% (NaN) bile Pmax TUM sweep araligindan araniyor ve PCE HER ZAMAN
% hesaplaniyordu -- bu, bu oturum boyunca birkac kez gorulen fiziksel
% olarak tutarsiz sonuclara yol acti (orn. RHJ ideal-arayuz: Voc=NaN,
% FF=NaN, ama eta=%4.16 raporlandi; BSF indoor: Voc=NaN ama eta=%43-48
% gibi guvenilmez sayilar). J-V egrisi sifir-akim gecisine ulasmamissa
% egri TAMAMLANMAMIStir ve Pmax gercek fotovoltaik guc maksimumu olarak
% GUVENILMEZ. Simdi Voc/Jsc gecersizse TUM cikislar (FF, PCE, Pmax) NaN
% ve m.valid=false -- cagiran kod (ozellikle otomatik karsilastirma/
% tarama scriptleri) artik bunu ACIKCA kontrol edebilir.
valid = ~isnan(Jsc) && ~isnan(Voc) && Voc > 0 && Jsc > 0;

if ~valid
  m = struct('Jsc',Jsc,'Voc',Voc,'FF',NaN,'PCE',NaN,'Pmax',NaN,'valid',false);
  return;
end

% Pmax SADECE 0<=V<=Voc, J>=0 bolgesinde aranir (Voc artik gecerli
% oldugu icin bu sinir guvenli).
Pmax = 0; kBest = -1;
for k = 1:numel(V)
  if V(k) < 0 || V(k) > Voc || J(k) < 0, continue; end
  P = V(k)*J(k);
  if P > Pmax, Pmax = P; kBest = k; end
end

% Dis inceleme madde 3.5: Pmax onceden SADECE mevcut (kaba) sweep grid
% noktalarindan seciliyordu -- gercek MPP grid noktalari arasinda
% kalirsa Pmax/FF hafifce dusuk/kararsiz cikabilirdi. Simdi kaba-grid
% maksimumunun ETRAFINDAKI 3 noktadan (varsa) yerel bir parabol
% gecirilip tepe noktasi kullaniliyor -- ucuz, ve grid adimindan
% bagimsiz bir MPP tahmini verir.
if kBest > 1 && kBest < numel(V)
  V1=V(kBest-1); V2=V(kBest); V3=V(kBest+1);
  P1=V1*J(kBest-1); P2=Pmax; P3=V3*J(kBest+1);
  % Standart parabolik-interpolasyon tepe formulu (Brent'in yontemindeki
  % gibi, sayisal optimizasyon ders kitaplarindan -- kendi turetmemden
  % daha az hataya acik).
  numer = (V2-V1)^2*(P2-P3) - (V2-V3)^2*(P2-P1);
  denom = (V2-V1)*(P2-P3) - (V2-V3)*(P2-P1);
  if abs(denom) > 0
    Vvertex = V2 - 0.5*numer/denom;
    if Vvertex > V1 && Vvertex < V3 % ekstrapolasyon degil, interpolasyon
      % Lagrange polinomu ile Vvertex'teki P degerini dogrudan degerlendir.
      Pvertex = P1*(Vvertex-V2)*(Vvertex-V3)/((V1-V2)*(V1-V3)) + ...
                P2*(Vvertex-V1)*(Vvertex-V3)/((V2-V1)*(V2-V3)) + ...
                P3*(Vvertex-V1)*(Vvertex-V2)/((V3-V1)*(V3-V2));
      if Pvertex > Pmax, Pmax = Pvertex; end
    end
  end
end

FF = Pmax/(Voc*Jsc);
PCE = Pmax * 1e4 / Pin_Wm2 * 100;

m = struct('Jsc',Jsc,'Voc',Voc,'FF',FF,'PCE',PCE,'Pmax',Pmax,'valid',true);
end
