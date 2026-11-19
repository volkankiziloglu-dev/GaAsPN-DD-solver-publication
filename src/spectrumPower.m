function P = spectrumPower(pts)
P = trapz(pts(:,1), pts(:,2));
end
