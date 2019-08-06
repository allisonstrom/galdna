function cloudy_match,gal,cloudy,err=err,logno=logno,noerr=noerr,logu=logu,uerr=uerr,zstar=zstar,zerr=zstarerr
  
  model = cloudy
  newgal = gal
  
  ebmv = 2.5*alog10((newgal.fha/newgal.fhb)/2.86)/(k_lambda_als(4861.36,/ccm)-k_lambda_als(6562.85,/ccm))
  if ebmv lt 0 then ebmv = 0
  
  r3a = newgal.fo3a/newgal.fhb*10^(ebmv/2.5*(k_lambda_als(5008.239,/ccm)-k_lambda_als(4364,/ccm)))
  if r3a le 0 then begin
     r3a_err = newgal.eo3a/newgal.fhb
     r3a = 0
  endif else r3a_err = abs(r3a*sqrt((newgal.eo3a/newgal.fo3a)^2.+(newgal.ehb/newgal.fhb)^2.))
  
  r3 = newgal.fo3n[1]/newgal.fhb*10^(ebmv/2.5*(k_lambda_als(5008.239,/ccm)-k_lambda_als(4861.36,/ccm)))
  if r3 le 0 then begin
     r3_err = newgal.eo3n[1]/newgal.fhb
     r3 = 0
  endif else r3_err = abs(r3*sqrt((newgal.eo3n[1]/newgal.fo3n[1])^2.+(newgal.ehb/newgal.fhb)^2.))
  
  n2 = newgal.fn2[1]/newgal.fhb*10^(ebmv/2.5*(k_lambda_als(6585.27,/ccm)-k_lambda_als(4861.36,/ccm)))
  if n2 le 0 then begin
     n2_err = newgal.en2[1]/newgal.fhb
     n2 = 0
  endif else n2_err = abs(n2*sqrt((newgal.en2[1]/newgal.fn2[1])^2.+(newgal.ehb/newgal.fhb)^2.))
  
  ne3 = newgal.fne3/newgal.fhb*10^(ebmv/2.5*(k_lambda_als(3869.81,/ccm)-k_lambda_als(4861.36,/ccm)))
  if ne3 le 0 then begin
     ne3_err = newgal.ene3/newgal.fhb
     ne3 = 0
  endif else ne3_err = abs(ne3*sqrt((newgal.ene3/newgal.fne3)^2.+(newgal.ehb/newgal.fhb)^2.))
  
  s2 = (newgal.fs2[0]+newgal.fs2[1])/newgal.fhb*10^(ebmv/2.5*(k_lambda_als(6732.68,/ccm)-k_lambda_als(4861.36,/ccm)))
  if s2 le 0 then begin
     s2_err = sqrt(newgal.es2[0]^2.+newgal.es2[1]^2.)/newgal.fhb
     s2 = 0
  endif else s2_err = abs(s2*sqrt((newgal.es2[0]^2.+newgal.es2[1]^2.)/(newgal.fs2[0]+newgal.fs2[1])^2.+(newgal.ehb/newgal.fhb)^2.))
  
  r2 = (newgal.fo2[0]+newgal.fo2[1])/newgal.fhb*10^(ebmv/2.5*(k_lambda_als(3726.03,/ccm)-k_lambda_als(4861.36,/ccm)))
  if r2 le 0 then begin
     r2_err = sqrt(newgal.eo2[0]^2.+newgal.eo2[0]^2.)/newgal.fhb
     r2 = 0
  endif else r2_err = abs(r2*sqrt((newgal.eo2[0]^2.+newgal.eo2[1]^2.)/(newgal.fo2[0]+newgal.fo2[1])^2.+(newgal.ehb/newgal.fhb)^2.))

  norange = indgen(146)/100.-1.8
  nscale = 10^(norange-(-0.86))
  
  ;; xi2 = -1000
  xi2arr = dblarr([n_elements(nscale),n_elements(model)])
  for n=0,n_elements(nscale)-1 do begin
     tmp = 0
     if r3a_err gt 0 then tmp += (r3a-model.fo3a)^2./r3a_err
     if r3_err gt 0 then tmp += (r3-model.fo3n[1])^2./r3_err
     if n2_err gt 0 then tmp += (n2-model.fn2[1]*nscale[n])^2./n2_err
     if ne3_err gt 0 then tmp += (ne3-model.fne3)^2./ne3_err
     if s2_err gt 0 then tmp += (s2-(model.fs2[0]+model.fs2[1]))^2./s2_err
     if r2_err gt 0 then tmp += (r2-model.fo2)^2./r2_err
     xi2arr[n,*] = tmp
  endfor
  wheretomulti,xi2arr,(sort(xi2arr))[0:4],noindex,metindex
  
  logu = model[metindex[0]].ion
  met = model[metindex[0]].met
  zstar = model[metindex[0]].teff ;/0.014
  logno = alog10(nscale[noindex[0]])-0.86
  return,met
end
  
