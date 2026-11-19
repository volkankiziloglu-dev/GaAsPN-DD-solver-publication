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

% FIX (external review item 2.1, CRITICAL): PREVIOUSLY, even if Voc/Jsc were not found
% (NaN), Pmax was searched over the ENTIRE sweep range and PCE was ALWAYS
% computed -- this led to physically inconsistent results seen several times
% during this session (e.g. RHJ ideal interface: Voc=NaN,
% FF=NaN, but eta=%4.16 reported; BSF indoor: Voc=NaN but eta=%43-48
% and similar unreliable numbers). If the J-V curve has not reached the zero-current crossing,
% the curve is INCOMPLETE and Pmax is UNRELIABLE as the true photovoltaic power
% maximum. Now, if Voc/Jsc are invalid, ALL outputs (FF, PCE, Pmax) are NaN
% and m.valid=false -- calling code (especially automatic comparison/
% sweep scripts) can now check this EXPLICITLY.
valid = ~isnan(Jsc) && ~isnan(Voc) && Voc > 0 && Jsc > 0;

if ~valid
  m = struct('Jsc',Jsc,'Voc',Voc,'FF',NaN,'PCE',NaN,'Pmax',NaN,'valid',false);
  return;
end

% Pmax is searched ONLY in the 0<=V<=Voc, J>=0 region (this bound is safe
% since Voc is now valid).
Pmax = 0; kBest = -1;
for k = 1:numel(V)
  if V(k) < 0 || V(k) > Voc || J(k) < 0, continue; end
  P = V(k)*J(k);
  if P > Pmax, Pmax = P; kBest = k; end
end

% External review item 3.5: Pmax was previously chosen ONLY from the existing (coarse) sweep grid
% points -- if the true MPP lies between grid points,
% Pmax/FF could come out slightly low/unstable. Now a local parabola
% is fitted through the 3 points AROUND the coarse-grid maximum (if available)
% and its vertex is used -- cheap, and gives an MPP estimate independent
% of the grid step.
if kBest > 1 && kBest < numel(V)
  V1=V(kBest-1); V2=V(kBest); V3=V(kBest+1);
  P1=V1*J(kBest-1); P2=Pmax; P3=V3*J(kBest+1);
  % Standard parabolic-interpolation vertex formula (as in Brent's method,
  % from numerical optimization textbooks -- less error-prone than
  % my own derivation).
  numer = (V2-V1)^2*(P2-P3) - (V2-V3)^2*(P2-P1);
  denom = (V2-V1)*(P2-P3) - (V2-V3)*(P2-P1);
  if abs(denom) > 0
    Vvertex = V2 - 0.5*numer/denom;
    if Vvertex > V1 && Vvertex < V3 % interpolation, not extrapolation
      % Evaluate P at Vvertex directly with the Lagrange polynomial.
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
