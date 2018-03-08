pro compare_stddev_snap
planetname = 'WASP-14b'
chname = '2'
apradius = 2.25
;run code to read in all the input planet parameters
  planetinfo = create_planetinfo()
  aorname= planetinfo[planetname, 'aorname_ch2'] 
  basedir = planetinfo[planetname, 'basedir']
  
  dirname = strcompress(basedir + planetname +'/');+'/hybrid_pmap_nn/')
  savefilename1 = strcompress(dirname + planetname +'_phot_ch'+chname+'_'+string(apradius)+'.sav',/remove_all)
  savefilename2 = strcompress(dirname + planetname +'_phot_ch'+chname+'_'+string(apradius)+'npmap.sav',/remove_all)
  savefilename3 = strcompress(dirname + planetname +'_phot_ch'+chname+'_'+string(apradius)+'dcorr.sav',/remove_all)
  savenames = [savefilename1, savefilename2, savefilename3]
  colorarr = ['black', 'blue', 'green']
  for s = 0, n_elements(savenames) - 1 do begin
     
     restore, savenames(s)
     print, 'restoring ', savenames(s)
     
     stddev_corrfluxarr = fltarr(n_elements(aorname))
     for a = 0,n_elements(aorname) - 1 do begin
        meanclip, planethash[aorname(a),'corrflux'] , mean_corrflux, stddev_corrflux
        stddev_corrflux = stddev(planethash[aorname(a),'corrflux'],/nan)
        stddev_corrfluxarr[a] = stddev_corrflux
     endfor
     
     plothist, stddev_corrfluxarr, xhist, yhist, bin = 1E-4, /noplot,/nan
     if s eq 0 then begin
        di = plot(xhist, yhist, thick = 2, xtitle = 'standard deviation of corrflux', ytitle = 'Number', title = planetname)
     endif else begin
        di = plot(xhist, yhist, thick = 2,color = colorarr(s), /overplot)
     endelse
     
  endfor
  
  
end
