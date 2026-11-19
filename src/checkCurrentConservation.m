function report = checkCurrentConservation(Jseg, Jterminal, JtolAbs, JtolRel)
%CHECKCURRENTCONSERVATION Verify steady-state total-current conservation.
if nargin < 3 || isempty(JtolAbs), JtolAbs = 1e-8; end
if nargin < 4 || isempty(JtolRel), JtolRel = 1e-3; end

Jscale = 1e-6;
deltaJabs = max(Jseg) - min(Jseg);
rJ = deltaJabs / max(abs(Jterminal), Jscale);

report.deltaJabs = deltaJabs;
report.rJ = rJ;
report.passAbs = deltaJabs < JtolAbs;
report.passRel = rJ < JtolRel;
report.pass = report.passAbs || report.passRel;
if report.pass
    report.reason = 'ok';
else
    report.reason = 'both_abs_and_rel_exceeded';
end
end
