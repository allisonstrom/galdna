pro galdna_example

  cloudy = mrdfits('cloudy_strom2018.fits',1,/silent)

  galaxy = mrdfits('kbss_lm1_lines.fits',1,/silent)
  ;; galaxy.field = 'LEGAC'
  ;; galaxy.obj = ''
  ;; galaxy.fo2[0] = 0.0
  ;; galaxy.eo2[0] = 0.0
  ;; galaxy.fo2[1] = 0.0
  ;; galaxy.eo2[1] = 0.0
  ;; galaxy.fne3 = 0.0
  ;; galaxy.ene3 = 0.0
  ;; galaxy.fhb = 0.0
  ;; galaxy.ehb = 0.0
  ;; galaxy.fo3a = 0.0
  ;; galaxy.eo3a = 0.0
  ;; galaxy.fo3n[0] = 0.0
  ;; galaxy.eo3n[0] = 0.0
  ;; galaxy.fo3n[1] = 0.0
  ;; galaxy.eo3n[1] = 0.0
  ;; galaxy.fha = 0.0
  ;; galaxy.eha = 0.0
  ;; galaxy.fn2[0] = 0.0
  ;; galaxy.en2[0] = 0.0
  ;; galaxy.fn2[1] = 0.0
  ;; galaxy.en2[1] = 0.0
  ;; galaxy.fs2[0] = 0.0
  ;; galaxy.es2[0] = 0.0
  ;; galaxy.fs2[1] = 0.0
  ;; galaxy.es2[1] = 0.0

  filename = 'chain_'+galaxy.field+'_'+galaxy.obj+'.dat'
  galdna,galaxy,cloudy,filename=filename

  readcol,filename,scale,a,logl,zneb,no_ratio,logu,zstar,format='X,D,D,D,D,D,D,D',/silent
  plothist,zneb+8.69,bin=0.01,xtit='12+log(O/H)',ytit='N',xrange=8.69+[-1,0.3],/xs

end
