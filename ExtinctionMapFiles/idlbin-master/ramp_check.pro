pro ramp_check

;for each AOR
;bin all fluxes into 4 bins (quarters)
;then bin all AORs together
;plot  data points.
;see if I can fit them with a slope.

  planetname = 'WASP-14b'
  planetinfo = create_planetinfo()
  aorname= planetinfo[planetname, 'aorname_ch2'] 
  basedir = planetinfo[planetname, 'basedir']
  period =  planetinfo[planetname, 'period']
  
  dirname = strcompress(basedir + planetname +'/')                                     
  ;;savefilename = strcompress(dirname + planetname +'_phot_ch2_2.25000_150226_bcdsdcorr.sav',/remove_all) ;
  savefilename = strcompress(dirname + planetname +'_phot_ch2_2.25000_160126.sav',/remove_all) ;
  restore, savefilename
  stareaor = 5
  startaor =  stareaor               
  stopaor =     n_elements(aorname) - 1
  

  for a = startaor, stopaor do begin
     corrflux = planethash[aorname(a),'corrflux']
     ;;get rid of NANs = uncorrectable positions
     n = where(finite(corrflux) gt 0, ngood)

     if ngood gt 0.2*(n_elements(corrflux)) then begin
        corrflux = corrflux(n)

        ;;dynamic arrays
        if a eq startaor then begin
           firstarr = corrflux[0:ngood/4 - 1]
           secondarr = corrflux[ngood/4:ngood/2 - 1]
           thirdarr = corrflux[ngood/2: 3*ngood /4 - 1]
           fourtharr = corrflux[3*ngood /4 - 1:ngood-1]
           nele = ngood/4

        endif else begin
           firstarr = [firstarr, corrflux[0:ngood/4 - 1]]
           secondarr = [secondarr, corrflux[ngood/4:ngood/2 - 1]]
           thirdarr = [thirdarr, corrflux[ngood/2: 3*ngood /4 - 1]]
           fourtharr = [fourtharr, corrflux[3*ngood /4 - 1:ngood-1]]
           nele = nele + ngood/4
           
        endelse
     endif

     
  endfor

  x = findgen(4) + 0.5
  meanclip, firstarr, m1, s1
  meanclip, secondarr, m2, s2
  meanclip, thirdarr, m3, s3
  meanclip, fourtharr, m4, s4
  y = [m1, m2, m3,m4]
  e = [s1, s2, s3, s4]
  e = e / sqrt(nele)
  p = errorplot(x, (y/mean(y) )  ,e/mean(y) , '1s', sym_filled = 1, xtitle = 'bins', $
                ytitle = 'Normalized Corrected Flux', margin = 0.2)

end


;     meanclip, corrflux[0:ngood/4 - 1], meanfirst, sigmafirst
;     firstarr[a] = meanfirst
;     meanclip, corrflux[ngood/4:ngood/2 - 1], meansecond, sigmasecond
;     secondarr[a] = meansecond
;     meanclip, corrflux[ngood/2: 3*ngood /4 - 1], meanthird, sigmathird
;     thirdarr[a] = meanthird
;     meanclip, corrflux[3*ngood /4 - 1:ngood-1], meanfourth, sigmafourth
;     fourtharr[a] = meanfourth
