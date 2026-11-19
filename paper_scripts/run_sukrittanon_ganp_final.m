% Final parameter set for the Sukrittanon GaNP replication (see results/sukrittanon_ganp/SUMMARY.md).
tB=1e-9; tau=5e-9; Rs=22;
fprintf('=== FINAL tau=%g Rs=%g tauBase=%g\n',tau,Rs,tB);
for ar=[false true]
 for s={'GaNP-2000','GaNP-800','GaP-800'}
  run_sukrittanon_ganp_validation('jv','sample',s{1},'tauBase',tB,'tauI',tau,'tauGaPi',tau,'Rs',Rs,'AR',ar,'tag','_final');
 end
end
lam=[400 450 500 530 560 580 600];
for s={'GaNP-2000','GaNP-800'}
  o=run_sukrittanon_ganp_validation('jsc','sample',s{1},'tauBase',tB,'tauI',tau,'AR',false);
  fprintf('EQEnoAR %s: ',s{1}); for l=lam, fprintf('%d:%.2f ',l,o.eqe(o.lam==l)); end; fprintf('\n');
end
