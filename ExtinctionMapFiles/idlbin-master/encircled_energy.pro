pro encircled_energy, planetname
; this program does photometry to get the encircled energy as a function of aperture size for a staring mode data set
; the purpose of this code is to look for changes in the slope of that curve as a function of time


  print, systime()

  ;ap = [1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,6.5,7.,7.5,8.,8.5,9.,9.5,10.]
  ap =  findgen(45) / 5 + 1.   ;this will make a smoother curve
  back = [11.,15.]              ; [11., 15.5]

  planetinfo = create_planetinfo()
  chname = planetinfo[planetname, 'chname']
  ra_ref = planetinfo[planetname, 'ra']
  dec_ref = planetinfo[planetname, 'dec']
  aorname = planetinfo[planetname, 'aorname']
  basedir = planetinfo[planetname, 'basedir']

;---------------
dirname = strcompress(basedir + planetname +'/')
  for a = 0, n_elements(aorname) - 1 do begin

     print, 'working on aorname ', aorname(a)
     dir = dirname+ string(aorname(a) )      
     CD, dir                    ; change directories to the corrct AOR directory
     command  = strcompress( 'find ch'+chname+"/bcd -name '*_bcd.fits' > "+dirname+'bcdlist.txt')
     spawn, command
     command2 =  strcompress('find ch'+chname+"/bcd -name '*bunc.fits' > "+dirname + 'bunclist.txt')
     spawn, command2

     readcol,strcompress(dirname +'bcdlist.txt'),fitsname, format = 'A', /silent
     readcol,strcompress(dirname+'bunclist.txt'),buncname, format = 'A', /silent

     nfits = n_elements(fisname)
     nfits = 1700L  ; just run the first half for speed

     xarr = dblarr(64*nfits)
     yarr = dblarr(64*nfits)
     skyarr = dblarr(64*nfits)
     fluxarr = dblarr(64*nfits,n_elements(ap))
     backarr = dblarr(64*nfits)
     backsigarr = dblarr(64*nfits)
     sclkarr = dblarr(64*nfits)
     prfarr = dblarr(64*nfits)

     i = 0L
     print, 'n fits', nfits
      
     for f =0.D,  nfits - 1 do begin ;read each bcd file, find centroid, keep track
        fits_read, fitsname(f), data, header  
        gain = sxpar(header, 'GAIN')
        fluxconv = sxpar(header, 'FLUXCONV')
        exptime = sxpar(header, 'EXPTIME')
        ronoise = sxpar(header, 'RONOISE')
        aorlabel = sxpar(header, 'OBJECT')
        sclk = sxpar(header, 'SCLK_OBS')

        convfac = gain*exptime/fluxconv
        adxy, header, ra_ref, dec_ref, xinit, yinit
        for j = 1, 63 do begin  ;don't use the first frame
           indim = data[*,*,j]
           indim = indim*convfac
           irac_box_centroider, indim, xinit, yinit, 3, 6, 3, ronoise,xcen, ycen, box_f, box_sky, box_c, box_cb, box_np,/MMM
           xarr[i] = xcen
           yarr[i] = ycen
           skyarr[i] = box_sky
           
;now run aper at the centroids found on a whole string of apertures
           aper, indim,xcen, ycen, xf, xfs, xb, xbs, 1.0, ap, back, $
                 /FLUX, /EXACT, /SILENT, /NAN, READNOISE=ronoise
           ;print, 'fluxes', xf
           fluxarr(i,*) = xf
           backarr[i] = xb
           backsigarr[i] =xbs
           sclkarr[i] = sclk
           i = i + 1
        endfor   ;for each frame

     endfor ;for each fits image

  endfor ; for each AOR

;print, 'fluxarr[*,1]', fluxarr[1,*]
;print, 'fluxarr', fluxarr

;  for p = 0, 9 do begin
;     a = plot(ap, fluxarr[p,*], '1r', xtitle = 'Aperture size', ytitle = 'Encircled Energy',/overplot)
;  endfor

;ok but now bin by some number of images to reduce noise in the curves
  numberarr = findgen(n_elements(xarr))
  bin_level = 2*63L
  h = histogram(numberarr, OMIN=om, binsize = bin_level, reverse_indices = ri)
  print, 'omin', om, 'nh', n_elements(h)

  bin_flux = dblarr(n_elements(h),n_elements(ap))
  for j = 0L, n_elements(h) - 1 do begin
     if j eq 0 then print, 'binning together', n_elements(numberarr[ri[ri[j]:ri[j+1]-1]])
        ;print, 'binning', numberarr[ri[ri[j]:ri[j+1]-1]]
     
     for k = 0, n_elements(ap) - 1 do begin
        meanclip, fluxarr[ri[ri[j]:ri[j+1]-1], k], meanflux, sigmax
        if j eq 0 then print, 'meanflux', meanflux
        bin_flux(j,k) = meanflux   
     endfor

  endfor


;test output 
  print, 'bin_flux', bin_flux(0,*)


;normalize the binned curves
  for n = 0, n_elements(h) - 1 do begin
     bin_flux(n,*) = bin_flux(n,*) / bin_flux(n, n_elements(ap)-1)
  endfor

;test output 
  print, 'bin_flux', bin_flux(0,*)


;Plot the binned curve of growths
  min = 400
  for p = min, n_elements(h) - 1 do begin
     if p eq min then begin
        a = plot(ap, bin_flux(p,*), '1r', xtitle = 'Aperture size', ytitle = 'Encircled Energy') 
     endif
     if p gt 594 and p lt 670 then begin
        a = plot(ap, bin_flux(p,*), '1b',/overplot)
     endif else begin
        a = plot(ap, bin_flux(p,*), '1r', /overplot)
     endelse
     
  endfor


  save, /all, filename =strcompress(dirname + planetname +'_ee_ch'+chname+'.sav')
  print, systime()

end
