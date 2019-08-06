pro galdna_example

  ;; This is the default photoionization model grid that was used in
  ;; Strom et al. (2018)
  cloudy = mrdfits('cloudy_strom2018.fits',1,/silent)

  ;; This structure contains the measurements for the KBSS composite
  ;; from Steidel et al. (2016)
  galaxy = mrdfits('kbss_lm1_lines.fits',1,/silent)
  
  ;; If you would like to use measurements for a different galaxy, you
  ;; will need to create a structure with the following elements. If
  ;; you do not want to use a specific line in the MCMC
  ;; (because it isn't detected, etc.), make sure that the
  ;; error for that line is '0.0'.

  ;; galaxy = {obj:'', fHa:0.0, eHa:0.0, fHb:0.0, eHb:0.0, $
  ;;           fO3n:fltarr(2), eO3n:fltarr(2), fO3a:0.0, eO3a:0.0, fN2:fltarr(2), eN2:fltarr(2), $
  ;;           fS2:fltarr(2), eS2:fltarr(2), fO2:fltarr(2), eO2:fltarr(2), fNe3:0.0, eNe3:0.0}
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

  filename = 'chain_'+galaxy.obj+'.dat'
  galdna,galaxy,cloudy,filename=filename

  readcol,filename,scale,a,logl,zneb,no_ratio,logu,zstar,format='X,D,D,D,D,D,D,D',/silent
  plothist,zneb+8.69,bin=0.01,xtit='12+log(O/H)',ytit='N',xrange=8.69+[-1,0.3],/xs

end
