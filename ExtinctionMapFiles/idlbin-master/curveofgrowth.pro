pro curveofgrowth, planetname, chname, corrflux = corrflux, selfcal=selfcal

  planetinfo = create_planetinfo()
  ra = planetinfo[planetname, 'ra']
  dec = planetinfo[planetname, 'dec']
  if chname eq '2'  then aorname= planetinfo[planetname, 'aorname_ch2'] else aorname = planetinfo[planetname, 'aorname_ch1'] 
  basedir = planetinfo[planetname, 'basedir']
  dirname = basedir+ planetname + '/'
  ; from get_centroids_for_calstar_jk
  aps1 = [ 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25]

  real_flux_ch1 = [8.68,647.38,38.20, 40.74,50.40]
  real_flux_ch2 = [5.66,421.19,24.74, 26.28,32.54]
  
  print, 'aorname', aorname

  for a = 1,  1 do begin;n_elements(aorname) - 1 do begin
     
     dirloc =dirname +string(aorname(a) ) 
     print, 'dirloc', dirloc
     cd, dirloc
     command  = strcompress( 'find ch'+chname+"/bcd -name 'sdcorr*_bcd.fits' > "+dirname+'bcdlist.txt')
     spawn, command
     command2 =  strcompress('find ch'+chname+"/bcd -name '*bunc.fits' > "+dirname + 'bunclist.txt')
     spawn, command2
        
     readcol,strcompress(dirname +'bcdlist.txt'),fitsname, format = 'A', /silent
     readcol,strcompress(dirname+'bunclist.txt'),buncname, format = 'A', /silent
        
;     fluxarr0 = fltarr((n_elements(fitsname)) - 1 * 63)

     for i =1, n_elements(fitsname) - 1 do begin
        fits_read,fitsname(i), im, h
        fits_read, buncname(i), unc, hunc
        chnlnum = sxpar(h, 'CHNLNUM')
        if planetname eq 'wasp14' then begin
           for j = 0, 63 do begin
              im[4:7, 13:16, j] = !Values.F_NAN ;mask region with nan set for bad regions                               
           endfor         
        endif
                                ;do the centroiding & photometry
        get_centroids_curveofgrowth,im, h, unc, ra, dec,  t, dt, hjd, xft, x3, y3, $
                                     x5, y5, x7, y7, xg, yg, xh, yh, f, b, x3s, y3s, x5s, y5s, $
                                     x7s, y7s, fs, bs, xp3, yp3, xp5, yp5, xp7, yp7, xp3s, yp3s, $
                                     xp5s, yp5s, xp7s, yp7s, fp, fps, np, flag, ns, sf, $
                                     xfwhm, yfwhm, /WARM

;        ch = 2

;run pmap correction
        if keyword_set(corrflux) then begin  ; run the pmap correction before making the signal to noise ratio
           f[*,0] = iracpc_pmap_corr( f[*,0] ,x3,y3,chnlnum,/threshold_occ, threshold_val = 20) 
           f[*,1] = iracpc_pmap_corr( f[*,1] ,x3,y3,chnlnum,/threshold_occ, threshold_val = 20) 
           f[*,2] = iracpc_pmap_corr( f[*,2] ,x3,y3,chnlnum,/threshold_occ, threshold_val = 20) 
           f[*,3] = iracpc_pmap_corr( f[*,3] ,x3,y3,chnlnum,/threshold_occ, threshold_val = 20) 
           f[*,4] = iracpc_pmap_corr( f[*,4] ,x3,y3,chnlnum,/threshold_occ, threshold_val = 20) 
           f[*,5] = iracpc_pmap_corr( f[*,5] ,x3,y3,chnlnum,/threshold_occ, threshold_val = 20) 
           f[*,6] = iracpc_pmap_corr( f[*,6] ,x3,y3,chnlnum,/threshold_occ, threshold_val = 20) 
           f[*,7] = iracpc_pmap_corr( f[*,7] ,x3,y3,chnlnum,/threshold_occ, threshold_val = 20) 
        endif 

        if keyword_set(selfcal) then begin
           ;XXXXX not finished coding, changing how to get these values outside of this code.
;get rid of position outliers  
           good = where(x3 lt mean(x3) + 3.0*stddev(x3) and x3 gt mean(x3) -3.0*stddev(x3) $
                        and y3 lt mean(y3) +3.0*stddev(y3) and y3 gt mean(y3) - 3.0*stddev(y3),$
                        ngood_pmap, complement=bad) 
           print, 'bad, good position',n_elements(bad), n_elements(good)
           x3[bad] = !Values.F_NAN
           y3[bad] = !Values.F_NAN
           f[bad] = !Values.F_NAN
           fs[bad] = !Values.F_NAN

  ;try getting rid of flux outliers.
  ;do some running mean with clipping
           start = 0
           for ni = 100, n_elements(f) -1,100 do begin
              meanclip,f[start:ni], m, s, subs = subs ;,/verbose
              
              if ni eq 100 then good_ni = subs+start else good_ni = [good_ni, subs+start]
              start = ni + 1
           endfor
           print, 'n good fluxclip', n_elements(good_ni), n_elements(f)
           x3 = x3[good_ni]
           y3 = y3[good_ni]
           f = f[good_ni]
           fs = fs[good_ni]

           amx = median(x3)         
           amy = median(y3)        
           adx = x3 - amx           
           ady = y3 - amy    
           
           nseg = 1
; Initial guesses	
           if keyword_set(sine) then begin ;use a fitting function with a sine curve or not
              known_period = period*24.    ; hours 
              pa0 = [median(y), 1.,1.,1.,1.,1.,0.001,known_period*(2*!Pi),0.9]
              func = 'fpa1_xfunc3'
           endif else begin
              pa0 = [median(y), 1.,1.,1.,1.,1.]
              func = 'fpa1_xfunc2'
              
           endelse
           
           parinfo = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0.D]}, n_elements(pa0))
                                ;limit the range of the sin cuve phase
           if keyword_set(sine) then begin ;use a fitting function with a sine curve
              parinfo[7].fixed = 1          
              parinfo[6].limited[0] = 1
              parinfo[6].limits[0] = 0.0
           endif
           
           
           afargs = {FLUX:y, DX:adx, DY:ady, T:x, ERR:yerr}
           pa = mpfit('fpa1_xfunc2', pa0, FUNCTARGS=afargs, PERROR=spa, BESTNORM=achi, DOF=adof, COVAR = COV, status = status, errmsg = errmsg, parinfo = parinfo) ;,/quiet)
           print, 'status', status
           print, errmsg
           achi = achi / adof
           print, 'reduced chi squared', achi
           if keyword_set(sine) then begin ;use a fitting function with a sine curve or not
              model_fit = pa[0] * (1+ pa[1] * adx + pa[2] * ady + pa[3] * adx * ady + pa[4] * adx * adx + pa[5] * ady * ady)+ pa[6]*sin(x/pa[7] + pa[8]) 
           endif else begin
              model_fit = pa[0] * (1+ pa[1] * adx + pa[2] * ady + pa[3] * adx * ady + pa[4] * adx * adx + pa[5] * ady * ady)
           endelse
           
        endif  ; keyword_set(selfcal)
        
        
        if i eq 1 then begin
           fluxarr0 = f[1:*,0] 
           fluxarr1 = f[1:*,1]
           fluxarr2 = f[1:*,2]
           fluxarr3 = f[1:*,3]
           fluxarr4 = f[1:*,4]
           fluxarr5 = f[1:*,5]
           fluxarr6 = f[1:*,6]
           fluxarr7 = f[1:*,7]
           uncarr0 = fs[1:*,0] 
           uncarr1 = fs[1:*,1]
           uncarr2 = fs[1:*,2]
           uncarr3 = fs[1:*,3]
           uncarr4 = fs[1:*,4]
           uncarr5 = fs[1:*,5]
           uncarr6 = fs[1:*,6]
           uncarr7 = fs[1:*,7]
        endif else begin
           fluxarr0 = [fluxarr0, f[1:*,0]]
           fluxarr1 = [fluxarr1, f[1:*,1]]
           fluxarr2 = [fluxarr2, f[1:*,2]]
           fluxarr3 = [fluxarr3, f[1:*,3]]
           fluxarr4 = [fluxarr4, f[1:*,4]]
           fluxarr5 = [fluxarr5, f[1:*,5]]
           fluxarr6 = [fluxarr6, f[1:*,6]]
           fluxarr7 = [fluxarr7, f[1:*,7]]
           uncarr0 = [uncarr0, fs[1:*,0]]
           uncarr1 = [uncarr1, fs[1:*,1]]
           uncarr2 = [uncarr2, fs[1:*,2]]
           uncarr3 = [uncarr3, fs[1:*,3]]
           uncarr4 = [uncarr4, fs[1:*,4]]
           uncarr5 = [uncarr5, fs[1:*,5]]
           uncarr6 = [uncarr6, fs[1:*,6]]
           uncarr7 = [uncarr7, fs[1:*,7]]
        endelse
        
     endfor                     ;end for each fits image
;     help, fluxarr0, fluxarr1, fluxarr2, fluxarr3
;     print, 'fluxarr0',fluxarr0[0:100]
     
     fluxarr = [[fluxarr0], [fluxarr1], [fluxarr2], [fluxarr3], [fluxarr4], [fluxarr5], [fluxarr6], [fluxarr7]]
     uncarr = [[uncarr0], [uncarr1], [uncarr2], [uncarr3], [uncarr4], [uncarr5], [uncarr6], [uncarr7]]
     help, fluxarr
     stdarr = fltarr(n_elements(aps1))
     meanarr = fltarr(n_elements(aps1))
     medarr = fltarr(n_elements(aps1))
     for ap = 0, n_elements(aps1) - 1 do begin
        meanclip, fluxarr[*,ap], meanflux, sigmacorr, clipsig = 3
        meanclip, uncarr[*,ap], meanunc, sigmaunc, clipsig = 3
        MMM, fluxarr[*,ap], skymde, sigma, skew
        stdarr[ap] = sigmacorr; meanunc
        meanarr[ap] = meanflux
        medarr[ap] = median(fluxarr[*,ap])
     endfor
           
     
  endfor                        ; for each aor
  
  print, 'stdarr', stdarr

  p = plot(aps1, meanarr/ stdarr, '1rs-', xtitle = 'aperture size', ytitle = 'SNR', title = planetname, xrange = [1.5, 3.5])
  p.save, dirname + 'curveofgrowth.png'

  save, aps1, meanarr, stdarr, filename = dirname + 'curveofgrowth_ch'+chname+'.sav' 
 ;      yerr = stdarr / (real_flux_ch1(r)*ratioarr) 
  ;       ploterror, aps1, ratioarr, yerr, yrange = [0.5, 1.2], psym = psymarr[r], xrange = [0,12], xtitle = 'aper radius', $
  ;                 ytitle = 'fraction of predicted flux', thick = 3, xthick = 3, ythick = 3, charthick = 3, title = 'PC5 - 20'
  ;      if r gt 0 then oploterror, aprad, ratioarr, yerr, psym = psymarr[r]
  
;---------------------------------------------------------------------------------------
;now redo the SNR and plotting for the corrected photometries (
;filename =strcompress(dirname + 'pixphasecorr_ch'+chname+'_'+aorname(a)+'.sav'
 
end
