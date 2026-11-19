function pts = buildCIEIndoorSpectrum(illum_mode, lux_target)
% Gercek CIE standart indoor aydinlatma spektrumu, deviceGenAndTransmit/
% spectrumPower'in bekledigi [lambda_nm, irradiance_Wm2nm] formatinda.
%
% Mantik get_CIE_irradiance.m (earlier GaAsNP project) ile AYNI -- veri dosyalari
% bu projenin data/ klasorune kopyalandi (loadRealAM15G.m ile ayni
% konvansiyon), boylece iki proje klasoru arasinda kirilgan bir bagimlilik
% olusmuyor.
%
% illum_mode: "2700K" | "4000K" | "6500K" (CIE LED-B1/B3/B5) | "4000F"
%             (CIE FL11 narrow-band fluorescent, ~4000 K)
% lux_target: hedef aydinlatma siddeti (lux). Tipik ofis: 300-750 lux.
%
% Donen pts(:,2) toplam W/m2/nm cinsinden GERCEK isinim yogunlugudur --
% buildIndoorLED.m'deki sentetik iki-Gaussian tahmininin YERINE gecer.

if nargin < 2 || isempty(lux_target), lux_target = 500; end

here = fileparts(mfilename('fullpath'));
V_data   = readmatrix(resolveDataFile('CIE_sle_photopic.csv'));
LED_data = readmatrix(resolveDataFile('CIE_illum_LEDs_1nm.csv'));
FL_data  = readmatrix(resolveDataFile('CIE_illum_FLs_1nm.csv'));

wl_nm = LED_data(:,1);
V_lambda = interp1(V_data(:,1), V_data(:,2), wl_nm, 'linear', 0);

switch illum_mode
  case "2700K"
    spec_raw = LED_data(:,2);  % LED-B1
  case "4000K"
    spec_raw = LED_data(:,4);  % LED-B3
  case "6500K"
    spec_raw = LED_data(:,6);  % LED-B5
  case "4000F"
    % Column 1 is wavelength; FL11 is therefore column 12. A legacy
    % implementation used column 20 (FL19) while labelling it F11-like.
    % Results made with that mapping are superseded and must not be mixed
    % with the production FL11 series.
    spec_raw = FL_data(:,12);  % CIE FL11 (formerly F11/TL84 class)
  otherwise
    error('buildCIEIndoorSpectrum: bilinmeyen illum_mode ''%s'' (2700K|4000K|6500K|4000F olmali)', illum_mode);
end

spec = spec_raw / trapz(wl_nm, spec_raw);
K = 683; % lm/W, standart fotopik luminous efficacy sabiti
lux_calc = K * trapz(wl_nm, spec .* V_lambda);
spec_lux = spec * (lux_target / lux_calc);

LER = lux_target / trapz(wl_nm, spec_lux); % = lux_calc (bkz. get_CIE_irradiance.m)
Pin_total = lux_target / LER;              % W/m2, hedef lux'u veren toplam radyant guc

irradiance = spec * Pin_total; % W/m2/nm
pts = [wl_nm, irradiance];
end
