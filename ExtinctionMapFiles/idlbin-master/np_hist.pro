pro np_hist
  planetname = 'WASP-52b'
  chname = '2'
  apradius = 2.25
  planetinfo = create_planetinfo()
  if chname eq '2' then aorname= planetinfo[planetname, 'aorname_ch2'] else aorname = planetinfo[planetname, 'aorname_ch1'] 
  basedir = planetinfo[planetname, 'basedir']
  period =  planetinfo[planetname, 'period']
  
  dirname = strcompress(basedir + planetname +'/')
  savefilename = strcompress(dirname + planetname +'_phot_ch'+chname+'_'+string(apradius)+'.sav',/remove_all)
  print, 'restoring ', savefilename
  restore, savefilename
  print, 'aorname', aorname(0)
  
  colorarr = ['blue', 'red','black','green','grey','purple']

  for a = 0, n_elements(aorname) -1 do begin
     if a eq 0 then begin
        meannp = mean(planethash[aorname(a),'np'],/nan)
        meannpcent = mean(planethash[aorname(a),'npcentroids'],/nan)
        delta = meannp - meannpcent
     endif

     np = planethash[aorname(a),'np'] - delta
     plothist, np, xhist, yhist, /noplot, bin = 0.05
;     if a eq 0 then p = plot (xhist, yhist, color = colorarr(a), xtitle = 'Noise pixel', $
;                              ytitle = 'Number', title = planetname, xrange = [3.5, 5.5], thick = 3)
;     if a gt 0 then p = plot (xhist, yhist, color = colorarr(a),thick = 3, /overplot)
;     
     plothist, planethash[aorname(a),'npcentroids'], xhist, yhist, /noplot, bin = 0.05
;     p = plot (xhist, yhist, color = colorarr(a),/overplot, linestyle = 2, thick = 3)
  endfor
;----------------------------------------------------

;just work on one AOR
  a = 1
  dir = dirname+ string(aorname(a) ) 
  CD, dir                       ; change directories to the correct AOR directory
  command  = strcompress( 'find ch'+chname+"/bcd -name 'SPITZER*_bcd.fits' > "+dirname+'bcdlist.txt')
  spawn, command
  command2 =  strcompress('find ch'+chname+"/bcd -name '*bunc.fits' > "+dirname + 'bunclist.txt')
  spawn, command2
  
  
  readcol,strcompress(dirname +'bcdlist.txt'),fitsname, format = 'A', /silent
  readcol,strcompress(dirname+'bunclist.txt'),buncname, format = 'A', /silent

  xcen = planethash[aorname(a),'xcen']
  ycen = planethash[aorname(a),'ycen']
  ra_ref = planethash[aorname(a),'ra']
  dec_ref = planethash[aorname(a),'dec']

  npaper = fltarr(n_elements(fitsname)*64)
  np41012 = npaper
  npgc = npaper
  npaper3 = npaper
  print, 'n_elements(fitsname_', n_elements(fitsname), n_elements(xcen) / 63, n_elements(ycen) / 63
  for i = 0L, n_elements(fitsname) - 1 do begin
     fits_read, fitsname(i), im, header
     fits_read, buncname(i), unc, hunc

     frametime = sxpar(header, 'framtime')
     ch = sxpar(header, 'CHNLNUM')
     ronoise = sxpar(header, 'RONOISE')  ; this is not realistic
     ronoise = 12.1   ; appropriate for 2s subarray, I think
     gain = sxpar(header, 'GAIN')
     fluxconv = sxpar(header, 'FLUXCONV')
     exptime = sxpar(header, 'EXPTIME')
     naxis = sxpar(header, 'NAXIS')
     
           ;run noisepix with various permutations on aperture size
           ;or anything else
           ;to see what the effect on the overall levels are.
     ; need to do this right
                                ;run get_centroids first, then use
                                ;those centroids to run noisepix
                                ;maybe the last run from
                                ;phot_exoplant, uses a different
                                ;radius in get_centroids

;this isn't going to be exact because since np is so sensitive
;to the inner radius, then it will also be sensitive to square
;vs. circular apertures. - turn off exact in aper call

; doesn't box_centroider add one to the box size?


 ;run the centroiding and photometry
      get_centroids_for_calstar_jk,im, header, unc, ra_ref, dec_ref,  t, dt, hjd, xft, x3, y3, $
                                   x5, y5, x7, y7, xg, yg, xh, yh, f, b, x3s, y3s, x5s, y5s, $
                                   x7s, y7s, fs, bs, xp3, yp3, xp5, yp5, xp7, yp7, xp3s, yp3s, $
                                   xp5s, yp5s, xp7s, yp7s, fp, fps, np, flag, ns, sf, $
                                   xfwhm, yfwhm, /WARM
      xcenim = temporary(x3)
      ycenim = temporary(y3)
      npgc(i*64) = np          ;keep track of np from get_centroids
      
      npaper(i*64) = noisepix( im, xcenim, ycenim, ronoise, gain, exptime, fluxconv,naxis,3.5,3.5, 6.5)
      np41012(i*64) = noisepix( im, xcenim, ycenim, ronoise, gain, exptime, fluxconv,naxis, 3.5, 6.5,9.5)
      npaper3(i*64) = noisepix( im, xcenim, ycenim, ronoise, gain, exptime, fluxconv,naxis, 3.5, 9.5, 12.5)
           

  endfor
  binval = 0.01
  plothist, npaper, xhistaper, yhistaper,bin = binval,/noplot
  plothist, np41012, xhist41012, yhist41012, bin = binval,/noplot
  plothist, npaper3, xhistaper3, yhistaper3, bin = binval,/noplot
  p2 = plot(xhistaper, yhistaper, color = 'red', xtitle = 'noise pixel', ytitle = 'Number', thick = 5, name = 'aper 35,3.5, 6.5' )
  p3 = plot(xhist41012, yhist41012, color = 'blue', /overplot, thick = 3, name = 'aper 3.5,6.5, 9.5')
  p4 = plot(xhistaper3, yhistaper3, color = 'cyan', /overplot, thick = 2, name = 'aper 3.5,6.5, 12.5')

;compare to get_centroids version with the same aperture  (3,3,6)
  print, 'n check', n_elements(npgc), n_elements(npaper), n_elements(np41012)
  plothist, npgc, xhistgc, yhistgc, bin = binval, /noplot
  p5 = plot(xhistgc, yhistgc, color = 'green', /overplot, name = 'gc 3,6, 3')

  l = legend(TARGET=[p2,p3,p4,p5], POSITION=[4.6,1400], /DATA, /AUTO_TEXT_COLOR)
end
