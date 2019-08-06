function rhat,chain,burnin=burnin

  if keyword_set(burnin) then burnin=burnin else burnin=0.1
   
  lchain = n_elements(chain)
  lburnin = long(burnin*lchain)
  nchain = max([ceil((lchain-lburnin)/32000.),10])
  lminichain = floor((lchain-lburnin)/nchain)
  chainarray = []
  means = []
  var = []
  for n=long(0),nchain-1 do begin
     index = long([n*lminichain+lburnin,(n+1)*lminichain-1+lburnin])
     chainarray = [[chainarray],[chain[index[0]:index[1]]]]
     means = [means,mean(chain[index[0]:index[1]],/double)]
     var = [var,variance(chain[index[0]:index[1]],/double)]
  endfor

  B = lminichain*variance(means,/double)
  W = mean(var,/double)
  R2 = (W*(lminichain-1)/lminichain+B/lminichain)/W
  return,sqrt(R2)

end
