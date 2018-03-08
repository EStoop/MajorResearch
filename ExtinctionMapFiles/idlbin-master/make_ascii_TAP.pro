pro make_ascii_TAP, bin_level

  planetname = 'WASP-15b'
  chname = '2'
  apradius = 1.75
  planetinfo = create_planetinfo()
  aorname= planetinfo[planetname, 'aorname_ch2'] 
  basedir = planetinfo[planetname, 'basedir']
  dirname = strcompress(basedir + planetname +'/') 
  
; read in the save files manually
  savefilename = strcompress(dirname+ '/corrflux1_noff.sav',/remove_all)
  restore, savefilename
  corrflux1_noff = corrflux1
  savefilename = strcompress(dirname+ '/corrflux2_noff.sav',/remove_all)
  restore, savefilename
  corrflux2_noff = corrflux2
  savefilename = strcompress(dirname+ '/corrflux3_noff.sav',/remove_all)
  restore, savefilename
  corrflux3_noff = corrflux3
  
  savefilename = strcompress(dirname+ '/corrflux1.sav',/remove_all)
  restore, savefilename
  savefilename = strcompress(dirname+ '/corrflux2.sav',/remove_all)
  restore, savefilename
  savefilename = strcompress(dirname+ '/corrflux3.sav',/remove_all)
  restore, savefilename
  
  savefilename = strcompress(dirname + 'pixphasecorr_ch'+chname+'_'+string(apradius)+'_xy.sav',/remove_all)
  restore, savefilename
  nnflux1 = flux
  savefilename = strcompress(dirname + 'pixphasecorr_ch'+chname+'_'+string(apradius)+'_np.sav',/remove_all)
  restore, savefilename
  nnflux2 = flux
  savefilename = strcompress(dirname + 'pixphasecorr_ch'+chname+'_'+string(apradius)+'_fw.sav',/remove_all)
  restore, savefilename
  nnflux3 = flux
  
;----------------------------------------------------------------------------------------------
;binning
  timearr =  (planethash[aorname(startaor),'timearr'] - (planethash[aorname(startaor),'timearr'])(0))
  bmjdarr = planethash[aorname(startaor),'bmjdarr'] 
  fluxarr = planethash[aorname(startaor),'flux'] 
  numberarr = findgen(n_elements(fluxarr))
  h = histogram(numberarr, OMIN=om, binsize = bin_level, reverse_indices = ri)
  
  bin_time = dblarr(n_elements(h))
  bin_flux = bin_time
  bin_bmjd = bin_time
  bin_corrflux1= bin_time
  bin_corrflux2= bin_time
  bin_corrflux3= bin_time
  bin_corrflux1_noff = bin_time
  bin_corrflux2_noff = bin_time
  bin_corrflux3_noff = bin_time
  bin_nn1 = bin_time
  bin_nn2 = bin_time
  bin_nn3 = bin_time
;  bin_corrfluxadd = bin_time
;  bin_corrfluxmult = bin_time

  c = 0
  for j = 0L, n_elements(h) - 1 do begin
     if (ri[j+1] gt ri[j] + 2)  then begin ;require 3 elements in the bin
        bin_time[c] = median(timearr[ri[ri[j]:ri[j+1]-1]])
        bin_bmjd[c] =  median(bmjdarr[ri[ri[j]:ri[j+1]-1]])
        meanclip_jk, fluxarr[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
        bin_flux[c] = meanx     ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
        meanclip_jk, corrflux1[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
        bin_corrflux1[c] = meanx ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
        meanclip_jk, corrflux2[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
        bin_corrflux2[c] = meanx ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
        meanclip_jk, corrflux3[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
        bin_corrflux3[c] = meanx ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
        meanclip_jk, corrflux1_noff[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
        bin_corrflux1_noff[c] = meanx ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
        meanclip_jk, corrflux2_noff[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
        bin_corrflux2_noff[c] = meanx ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
        meanclip_jk, corrflux3_noff[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
        bin_corrflux3_noff[c] = meanx ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
        meanclip_jk, nnflux1[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
        bin_nn1[c] = meanx      ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
        meanclip_jk, nnflux2[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
        bin_nn2[c] = meanx      ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
        meanclip_jk, nnflux3[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
        bin_nn3[c] = meanx      ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
        
        c = c + 1
     endif
  endfor
  bin_time = bin_time[0:c-1]
  bin_bmjd = bin_bmjd[0:c-1]
  bin_flux = bin_flux[0:c-1]
  bin_corrflux1 = bin_corrflux1[0:c-1]
  bin_corrflux2 = bin_corrflux2[0:c-1]
  bin_corrflux3 = bin_corrflux3[0:c-1]
  bin_corrflux1_noff = bin_corrflux1_noff[0:c-1]
  bin_corrflux2_noff = bin_corrflux2_noff[0:c-1]
  bin_corrflux3_noff = bin_corrflux3_noff[0:c-1]
  bin_nn1 = bin_nn1[0:c-1]
  bin_nn2 = bin_nn2[0:c-1]
  bin_nn3 = bin_nn3[0:c-1]

;----------------------------------------------------------------------------------------------
; print out ascii files for use in TAP
  outfilename =  [strcompress(dirname +'phot_TAP_pmap_xy.ascii',/remove_all), strcompress(dirname +'phot_TAP_pmap_np.ascii',/remove_all), strcompress(dirname +'phot_TAP_pmap_fw.ascii',/remove_all), strcompress(dirname +'phot_TAP_pmap_noff_xy.ascii',/remove_all), strcompress(dirname +'phot_TAP_pmap_noff_np.ascii',/remove_all), strcompress(dirname +'phot_TAP_pmap_noff_fw.ascii',/remove_all), strcompress(dirname +'phot_TAP_nn_xy.ascii',/remove_all), strcompress(dirname +'phot_TAP_nn_np.ascii',/remove_all), strcompress(dirname +'phot_TAP_nn_fw.ascii',/remove_all)]


  fluxname = [bin_corrflux1, bin_corrflux2, bin_corrflux3, bin_corrflux1_noff, bin_corrflux2_noff, $
              bin_corrflux3_noff, bin_nn1, bin_nn2, bin_nn3]
  endval = fix(0.3*n_elements(bin_nn1))
  start = 0
  print, 'n corrflux1', n_elements(bin_corrflux1)

  for fo = 0, n_elements(outfilename) - 1 do begin
     openw, outlun, outfilename(fo),/GET_LUN
     cb = 0
;     print, 'fo', fo, fo*n_elements(bin_corrflux1) , fo*n_elements(bin_corrflux1)  + endval
     normfactor = median(fluxname[fo*n_elements(bin_corrflux1):fo*n_elements(bin_corrflux1) + endval])
;     print, 'normfactor ', normfactor
        
     for te = 0,n_elements(bin_corrflux1) -1 do begin
;        print, 'te, ', te, fo*n_elements(bin_corrflux1) + te, fluxname[fo*n_elements(bin_corrflux1) + te]/normfactor
        printf, outlun, bin_bmjd(cb) , ' ', fluxname[fo*n_elements(bin_corrflux1) + te]/normfactor,  format = '(F0, A,D0)'
        cb = cb + 1
     endfor
     close, outlun
     free_lun, outlun
  endfor
  
  
end
