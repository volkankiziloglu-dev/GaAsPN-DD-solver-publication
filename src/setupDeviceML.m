function dev = setupDeviceML(layers, absorberIdx, junctions, opts)
% RearHJT/setupDeviceML.m'nin kopyasi -- gummelSolveML_v2 (eps0-duzeltilmis)
% cagiriyor.
c = constants();
dev = buildDeviceML(layers, absorberIdx, junctions, opts);

psiRough = solveEquilibrium(dev); % common/ -- SADECE baslangic tahmini
chiRel = dev.chi-dev.chiRef;
nRough = dev.Nc .* exp((psiRough+chiRel)/c.VT);
pRough = dev.Nv .* exp(-(psiRough+chiRel+dev.Eg)/c.VT);

dev.psiEq = psiRough;
r = gummelSolveML(dev, 0.0, zeros(dev.N,1), struct('psiInit',psiRough,'nInit',nRough,'pInit',pRough));
dev.psiEq = r.psi;
dev.nEq = r.n;
dev.pEq = r.p;
end
