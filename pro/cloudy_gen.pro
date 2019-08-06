function cloudy_gen,cloudy,input

  interpmet = 10^input[0]
  interpion = input[2]
  interpzstar = 10^input[3]*0.014

  met = cloudy[where(cloudy.ion eq min(cloudy.ion) and cloudy.teff eq min(cloudy.teff))].met
  nmet = n_elements(met)  
  ion = cloudy[where(cloudy.met eq min(cloudy.met) and cloudy.teff eq min(cloudy.teff))].ion
  nion = n_elements(ion)
  zstar = cloudy[where(cloudy.met eq min(cloudy.met) and cloudy.ion eq min(cloudy.ion))].teff
  nzstar = n_elements(zstar)
  
  upindex = (where(met eq min(met[where(met-interpmet ge 0)])))[0]
  downindex = (where(met eq max(met[where(met-interpmet le 0)])))[0]
  if upindex eq downindex then interpmetindex = downindex else interpmetindex = (interpmet-met[downindex])/abs(met[upindex]-met[downindex])+min([upindex,downindex])
  ;; print,met[downindex],downindex
  ;; print,interpmet,interpmetindex
  ;; print,met[upindex],upindex

  upindex = (where(ion eq min(ion[where(ion-interpion ge 0)])))[0]
  downindex = (where(ion eq max(ion[where(ion-interpion le 0)])))[0]
  if upindex eq downindex then interpionindex = downindex else interpionindex = (ion[downindex]-interpion)/abs(ion[upindex]-ion[downindex])+max([upindex,downindex])
  ;; print,ion[downindex],downindex
  ;; print,interpion,interpionindex
  ;; print,ion[upindex],upindex

  upindex = (where(zstar eq min(zstar[where(zstar-interpzstar ge 0)])))[0]
  downindex = (where(zstar eq max(zstar[where(zstar-interpzstar le 0)])))[0]
  if upindex eq downindex then interpzstarindex = downindex else interpzstarindex = (interpzstar-zstar[downindex])/abs(zstar[upindex]-zstar[downindex])+min([upindex,downindex])
  ;; print,zstar[downindex],downindex
  ;; print,interpzstar,interpzstarindex
  ;; print,zstar[upindex],upindex

  newindices = [interpmetindex,interpionindex,interpzstarindex]

  if ~file_test('cloudy.sav') then begin
     ha = dblarr(nmet,nion,nzstar)
     hb = dblarr(nmet,nion,nzstar)
     o3 = dblarr(nmet,nion,nzstar)
     o4363 = dblarr(nmet,nion,nzstar)
     o1661 = dblarr(nmet,nion,nzstar)
     n2 = dblarr(nmet,nion,nzstar)
     s2a = dblarr(nmet,nion,nzstar)
     s2b = dblarr(nmet,nion,nzstar)
     ;; s3n1 = dblarr(nmet,nion,nzstar)
     ;; s3n2 = dblarr(nmet,nion,nzstar)
     ;; s3a = dblarr(nmet,nion,nzstar)
     o2tot = dblarr(nmet,nion,nzstar)
     ne3 = dblarr(nmet,nion,nzstar)
     ;; he2opt = dblarr(nmet,nion,nzstar)
     for i=0,nmet-1 do begin
        for j=0,nion-1 do begin
           ha[i,j,*] = cloudy[where(cloudy.met eq met[i] and cloudy.ion eq ion[j])].fha
           hb[i,j,*] = cloudy[where(cloudy.met eq met[i] and cloudy.ion eq ion[j])].fhb
           o3[i,j,*] = cloudy[where(cloudy.met eq met[i] and cloudy.ion eq ion[j])].fo3n[1]
           o4363[i,j,*] = cloudy[where(cloudy.met eq met[i] and cloudy.ion eq ion[j])].fo3a
           o1661[i,j,*] = cloudy[where(cloudy.met eq met[i] and cloudy.ion eq ion[j])].fo3uv
           n2[i,j,*] = cloudy[where(cloudy.met eq met[i] and cloudy.ion eq ion[j])].fn2[1]
           s2a[i,j,*] = cloudy[where(cloudy.met eq met[i] and cloudy.ion eq ion[j])].fs2[0]
           s2b[i,j,*] = cloudy[where(cloudy.met eq met[i] and cloudy.ion eq ion[j])].fs2[1]
           ;; s3n1[i,j,*] = cloudy[where(cloudy.met eq met[i] and cloudy.ion eq ion[j])].fs3n[0]
           ;; s3n2[i,j,*] = cloudy[where(cloudy.met eq met[i] and cloudy.ion eq ion[j])].fs3n[1]
           ;; s3a[i,j,*] = cloudy[where(cloudy.met eq met[i] and cloudy.ion eq ion[j])].fs3a
           o2tot[i,j,*] = cloudy[where(cloudy.met eq met[i] and cloudy.ion eq ion[j])].fo2
           ne3[i,j,*] = cloudy[where(cloudy.met eq met[i] and cloudy.ion eq ion[j])].fne3
           ;; he2opt[i,j,*] = cloudy[where(cloudy.met eq met[i] and cloudy.ion eq ion[j])].fhe2opt
        endfor
     endfor
     save,cloudy,filename='cloudy.sav'
     save,ha,filename='cloudy_ha.sav'
     save,hb,filename='cloudy_hb.sav'
     save,o3,filename='cloudy_o3.sav'
     save,o4363,filename='cloudy_o4363.sav'
     save,o1661,filename='cloudy_o1661.sav'
     save,n2,filename='cloudy_n2.sav'
     save,s2a,filename='cloudy_s2a.sav'
     save,s2b,filename='cloudy_s2b.sav'
     save,o2tot,filename='cloudy_o2tot.sav'
     save,ne3,filename='cloudy_ne3.sav'
  endif else begin
     restore,'cloudy.sav'
     restore,'cloudy_ha.sav'
     restore,'cloudy_hb.sav'
     restore,'cloudy_o3.sav'
     restore,'cloudy_o4363.sav'
     restore,'cloudy_o1661.sav'
     restore,'cloudy_n2.sav'
     restore,'cloudy_s2a.sav'
     restore,'cloudy_s2b.sav'
     restore,'cloudy_o2tot.sav'
     restore,'cloudy_ne3.sav'
  endelse
  
  interp_ha = interpolate(ha,newindices[0],newindices[1],newindices[2]);,/double)
  interp_hb = interpolate(hb,newindices[0],newindices[1],newindices[2]);,/double)
  interp_o3 = interpolate(o3,newindices[0],newindices[1],newindices[2]);,/double)
  interp_o4363 = interpolate(o4363,newindices[0],newindices[1],newindices[2]);,/double)
  interp_o1661 = interpolate(o1661,newindices[0],newindices[1],newindices[2]);,/double)
  interp_n2 = interpolate(n2,newindices[0],newindices[1],newindices[2]);,/double)
  interp_s2a = interpolate(s2a,newindices[0],newindices[1],newindices[2]);,/double)
  interp_s2b = interpolate(s2b,newindices[0],newindices[1],newindices[2]);,/double)
  interp_o2tot = interpolate(o2tot,newindices[0],newindices[1],newindices[2]);,/double)
  interp_ne3 = interpolate(ne3,newindices[0],newindices[1],newindices[2]);,/double)
  ;; interp_s3n1 = interpolate(s3n1,newindices[0],newindices[1],newindices[2]) ;,/double)
  ;; interp_s3n2 = interpolate(s3n2,newindices[0],newindices[1],newindices[2]) ;,/double)
  ;; interp_s3a = interpolate(s3a,newindices[0],newindices[1],newindices[2]) ;,/double)
  ;; interp_he2opt = interpolate(he2opt,newindices[0],newindices[1],newindices[2]) ;,/double)
  
  str={teff:0.0, ion:0.0, met:0.0, $
       fHa:0.0, fHb:0.0, fHg:0.0, fHd:0.0, $
       fO3n:fltarr(2), fO3a:0.0, fN2:fltarr(2), $
       fS2:fltarr(2), fO2:0.0, fNe3:0.0, $
       fO3uv:0.0, fHe2uv:0.0, fHe2opt:0.0, fSi3:0.0, fC3:0.0, fs3n:fltarr(2), fs3a:0.0 }

  cat=replicate(str,1)
  cat.teff=interpzstar
  cat.ion=interpion
  cat.met=interpmet
  cat.fha=interp_ha
  cat.fhb=interp_hb
  cat.fo3n[1]=interp_o3
  cat.fo3a=interp_o4363
  cat.fo3uv=interp_o1661
  cat.fn2[1]=interp_n2
  cat.fs2[0]=interp_s2a
  cat.fs2[1]=interp_s2b
  ;; cat.fs3n[0]=interp_s3n1
  ;; cat.fs3n[1]=interp_s3n2
  ;; cat.fs3a=interp_s3a
  cat.fo2=interp_o2tot
  cat.fne3=interp_ne3
  ;; cat.fhe2opt=interp_he2opt

  return,cat
  
end
