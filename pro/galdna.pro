function logl,newgal,model,logno,dered

  if ~keyword_set(dered) then dered=0
  
  ebmv = 2.5*alog10((newgal.fha/newgal.fhb)/model.fha)/(k_lambda(4861.36,/ccm)-k_lambda(6562.85,/ccm))
  ebmv_err = 2.5/(k_lambda(4861.36,/ccm)-k_lambda(6562.85,/ccm))*sqrt((newgal.eha/newgal.fha)^2.+(newgal.ehb/newgal.fhb)^2.)/alog(10.)
  if ebmv lt 0 or dered then ebmv = 0

  dust = 10^(ebmv/2.5*(k_lambda(4364,/ccm)-k_lambda(4861.36,/ccm)))
  dust_err = abs(dust*((k_lambda(4364,/ccm)-k_lambda(4861.36,/ccm))/2.5*alog(10.)*ebmv_err))
  r3a = newgal.fo3a/newgal.fhb*dust
  if r3a le 0 then begin
     r3a_err = newgal.eo3a/newgal.fhb
     r3a = 0
  endif else r3a_err = abs(r3a*sqrt((newgal.eo3a/newgal.fo3a)^2.+(newgal.ehb/newgal.fhb)^2.+(dust_err/dust)^2.))

  dust = 10^(ebmv/2.5*(k_lambda(5008.239,/ccm)-k_lambda(4861.36,/ccm)))
  dust_err = abs(dust*((k_lambda(5008.239,/ccm)-k_lambda(4861.36,/ccm))/2.5*alog(10.)*ebmv_err))
  r3 = newgal.fo3n[1]/newgal.fhb*dust
  if r3 le 0 then begin
     r3_err = newgal.eo3n[1]/newgal.fhb
     r3 = 0
  endif else r3_err = abs(r3*sqrt((newgal.eo3n[1]/newgal.fo3n[1])^2.+(newgal.ehb/newgal.fhb)^2.+(dust_err/dust)^2.))

  dust = 10^(ebmv/2.5*(k_lambda(6585.27,/ccm)-k_lambda(4861.36,/ccm)))
  dust_err = abs(dust*((k_lambda(6585.27,/ccm)-k_lambda(4861.36,/ccm))/2.5*alog(10.)*ebmv_err))
  n2 = newgal.fn2[1]/newgal.fhb*dust
  if n2 le 0 then begin
     n2_err = newgal.en2[1]/newgal.fhb
     n2 = 0
  endif else n2_err = abs(n2*sqrt((newgal.en2[1]/newgal.fn2[1])^2.+(newgal.ehb/newgal.fhb)^2.+(dust_err/dust)^2.))
  
  dust = 10^(ebmv/2.5*(k_lambda(3869.81,/ccm)-k_lambda(4861.36,/ccm)))
  dust_err = abs(dust*((k_lambda(3869.81,/ccm)-k_lambda(4861.36,/ccm))/2.5*alog(10.)*ebmv_err))
  ne3 = newgal.fne3/newgal.fhb*dust
  if ne3 le 0 then begin
     ne3_err = newgal.ene3/newgal.fhb
     ne3 = 0
  endif else ne3_err = abs(ne3*sqrt((newgal.ene3/newgal.fne3)^2.+(newgal.ehb/newgal.fhb)^2.+(dust_err/dust)^2.))
  
  dust = 10^(ebmv/2.5*(k_lambda(6732.68,/ccm)-k_lambda(4861.36,/ccm)))
  dust_err = abs(dust*((k_lambda(6732.68,/ccm)-k_lambda(4861.36,/ccm))/2.5*alog(10.)*ebmv_err))
  s2 = (newgal.fs2[0]+newgal.fs2[1])/newgal.fhb*dust
  if s2 le 0 then begin
     s2_err = sqrt(newgal.es2[0]^2.+newgal.es2[1]^2.)/newgal.fhb
     s2 = 0
  endif else s2_err = abs(s2*sqrt((newgal.es2[0]^2.+newgal.es2[1]^2.)/(newgal.fs2[0]+newgal.fs2[1])^2.+(newgal.ehb/newgal.fhb)^2.+(dust_err/dust)^2.))
  
  dust = 10^(ebmv/2.5*(k_lambda(3726.03,/ccm)-k_lambda(4861.36,/ccm)))
  dust_err = abs(dust*((k_lambda(3726.03,/ccm)-k_lambda(4861.36,/ccm))/2.5*alog(10.)*ebmv_err))
  r2 = (newgal.fo2[0]+newgal.fo2[1])/newgal.fhb*dust
  if r2 le 0 then begin
     r2_err = sqrt(newgal.eo2[0]^2.+newgal.eo2[0]^2.)/newgal.fhb
     r2 = 0
  endif else r2_err = abs(r2*sqrt((newgal.eo2[0]^2.+newgal.eo2[1]^2.)/(newgal.fo2[0]+newgal.fo2[1])^2.+(newgal.ehb/newgal.fhb)^2.+(dust_err/dust)^2.))

  logl = 0
  nscale = 10^(logno+0.86)
  if r3a_err gt 0 then logl -= (r3a-model.fo3a)^2./(2*r3_err^2.)
  if r3_err gt 0 then logl -= (r3-model.fo3n[1])^2./(2*r3_err^2.)
  if n2_err gt 0 then logl -= (n2-model.fn2[1]*nscale)^2./(2*n2_err^2.)
  if ne3_err gt 0 then logl -= (ne3-model.fne3)^2./(2*ne3_err^2.)
  if s2_err gt 0 then logl -= (s2-(model.fs2[0]+model.fs2[1]))^2./(2*s2_err^2.)
  if r2_err gt 0 then logl -= (r2-model.fo2)^2./(2*r2_err^2.)
  
  ;; return log likelihood
  return,logl
end

function cloudy_init,newgal,model,dered

  if ~keyword_set(dered) then dered=0

  common ofebounds, lowofe, highofe
  trimmodel = model[where(alog10(model.met)-alog10(model.teff/0.014) ge lowofe and alog10(model.met)-alog10(model.teff/0.014) le highofe)]

  norange = indgen(16)/10.-2.0
  loglarr = dblarr([n_elements(norange),n_elements(trimmodel)])
  for i=0,n_elements(norange)-1 do for j=0,n_elements(trimmodel)-1 do loglarr[i,j] = logl(newgal,trimmodel[j],norange[i],dered)
  wheretomulti,loglarr,(reverse(sort(loglarr)))[0:4],noindex,metindex
  
  logu = trimmodel[metindex[0]].ion
  met = trimmodel[metindex[0]].met
  zstar = trimmodel[metindex[0]].teff
  logno = norange[noindex[0]] & if logno lt -1.5 then logno = -1.5

  if logno gt (alog10(met)+8.69)-9.2 then begin
     newmet = logno+9.2-8.69+0.2
     zdiff = newmet-alog10(zstar/0.014)
     if zdiff lt lowofe or zdiff gt highofe then begin
        zdiff = alog10(met)-alog10(zstar/0.014)
        zstar = 10^(newmet-zdiff)*0.014
     endif
     met = 10^newmet
  endif

  return,[alog10(met),logno,logu,alog10(zstar/0.014)]
end
  
pro galdna,gal,cloudy,maxiter=maxiter,rhat_target=rhat_target,restart=restart,filename=filename,silent=silent,dered=dered

  if ~keyword_set(dered) then dered=0
  
  tic
  common ofebounds, lowofe, highofe
  lowofe = 0.0
  highofe = 0.73

  resolve_all,/continue_on_error,/quiet,skip='rsex'
  if not keyword_set(maxiter) then maxiter = 100000
  if not keyword_set(rhat_target) then rhat_target = 1.05 
  if keyword_set(restart) then restart = 1 else restart = 0
  if keyword_set(silent) then silent = 1 else silent = 0

  ;; save cloudy emission line tables
  met = cloudy[where(cloudy.ion eq min(cloudy.ion) and cloudy.teff eq min(cloudy.teff))].met
  nmet = n_elements(met)  
  ion = cloudy[where(cloudy.met eq min(cloudy.met) and cloudy.teff eq min(cloudy.teff))].ion
  nion = n_elements(ion)
  zstar = cloudy[where(cloudy.met eq min(cloudy.met) and cloudy.ion eq min(cloudy.ion))].teff
  nzstar = n_elements(zstar)
  
  ;; initialize chain
  nparam = 4
  newgal = gal
  format = '(%"%6d %f %10.3e %f '
  for n=0,nparam-1 do format += '%7.4f '
  format += '")'
  scale = 0.1

  if not restart then begin
     s = strsplit(systime(),/extract)
     if ~keyword_set(filename) then begin
        files = file_search('chain_'+s[2]+strlowcase(s[1])+strmid(s[-1],2,2)+'*.dat')
        if files[0] eq '' then filename = 'chain_'+s[2]+strlowcase(s[1])+strmid(s[-1],2,2)+'_1.dat'
        if files[0] ne '' then begin
           filenum = fix((strsplit((strsplit(files[-1],'_',/extract))[-1],'.',/extract))[0])
           filenum++
           filenum = strtrim(string(filenum),1)
           filename = 'chain_'+s[2]+strlowcase(s[1])+strmid(s[-1],2,2)+'_'+filenum+'.dat'
        endif
     endif
     openw,1,filename

     init = cloudy_init(newgal,cloudy)
     
     chain = init
     model = cloudy_gen(cloudy,init)
     logpchain = [logl(newgal,model,init[1],dered)]
     printf,1,0,scale,1,logpchain[0],chain[*,0],format=format
     n = long(1)
     clock = tic('for 100 iterations ('+string(100*float(100)/maxiter,format='(F0.1)')+'% complete)')
  endif else if restart then begin
     openr,1,filename
     lines = file_lines(filename)
     chain = dblarr(nparam,lines)
     logpchain = dblarr(lines)
     for i=0,lines-1 do begin
        s = ''
        readf,1,s
        restartnum = (strsplit(s,' ',/extract))[0]
        scale = (strsplit(s,' ',/extract))[1]
        logpchain[i] = (strsplit(s,' ',/extract))[3]
        chain[*,i] = (strsplit(s,' ',/extract))[4:-1]
     endfor
     close,1
     n = long((strsplit(s,' ',/extract))[0])+1
     openw,1,filename,/append
  endif
  
  ;; begin MCMC loop
  accept = long(0)
  while n lt maxiter do begin
     ;; normal proposal distribution with the current state as the mean
     newchain = chain[*,-1]+scale*randomn(seed,nparam)
     a = 1
     if newchain[0] lt alog10(min(met)) or newchain[0] gt alog10(max(met)) then a = 0
     if newchain[1] lt -2.0 then a = 0
     if newchain[2] lt min(ion) or newchain[2] gt max(ion) then a = 0
     if newchain[3] lt alog10(min(zstar)/0.014) or newchain[3] gt alog10(max(zstar)/0.014) then a = 0
     if newchain[0]-newchain[3] gt highofe then a = 0
     if newchain[0]-newchain[3] lt lowofe then a = 0
     ;; Metropolis-Hastings algorithm
     if a gt 0 then begin
        new_logp = logl(newgal,cloudy_gen(cloudy,newchain),newchain[1],dered)
        a = exp(new_logp-logpchain[-1])
     endif
     thresh = randomu(seed)
     if thresh lt a then begin
        chain = [[chain],[newchain]]
        logpchain = [logpchain,new_logp]
        accept++
     endif else begin
        chain = [[chain],[chain[*,-1]]]
        logpchain = [logpchain,logpchain[-1]]
     endelse
     
     ;; check recent acceptance rate, adjust width of the proposal distribution
     if (n+1) mod 100 eq 0 then begin
        if min(logpchain) eq max(logpchain) then begin
           ;; message,'Chains are stationary'
           print,'WARNING: Chains are stationary, skipping object.'
           n = maxiter
        endif           
        rate = float(accept)/100
        if rate gt 0.234 then scale *= 1.1 else scale /= 1.1
        accept = long(0)
        if ~silent then begin
           toc,clock
           clock = tic('for 100 iterations ('+string(100*float(n+100)/maxiter,format='(F0.1)')+'% complete, logL='+strtrim(string(logpchain[-1]),1)+', '+string(rate,'(F0.2)')+' acceptance)')
        endif
     endif

     ;; check modified Gelman-Rubin statistic, exit if all chains are converged
     write = 1
     stop = 0
     if (n+1) mod 10000 eq 0 or restart then begin
        current_rhat = fltarr(nparam)
        for m=0,nparam-1 do begin
           tmpchain = chain[m,*]
           current_rhat[m] = rhat(tmpchain)
        endfor
        index = where(current_rhat ge rhat_target,count)
        if ((count eq 0) or (count eq 1 and index[0] eq 1 and gal.en2[1] eq 0)) then begin
           if restart then write = 0
           stop = 1
           print,'Convergence criterion acheived for all chains, ending MCMC.'
        endif
        restart = 0
     endif     
     if write then printf,1,n,scale,a,logpchain[-1],chain[*,-1],format=format
     if stop then n = maxiter
     
     n++
  endwhile
  toc
  close,1
  spawn,'/bin/rm *sav'
  
end

  
