pro test_bkgd, planetname, chname, apradius, nn

  apval = 2 ; aperture radius of 2.25
  pmapfile = '/Users/jkrick/irac_warm/pcrs_planets/pmap_phot/pmap_data_ch2_r2p25_s3_7_0p1s_x4_140904.sav'

  planetinfo = create_planetinfo()
  if chname eq '2' then aorname= planetinfo[planetname, 'aorname_ch2'] else aorname = planetinfo[planetname, 'aorname_ch1'] 
  basedir = planetinfo[planetname, 'basedir']
  period =  planetinfo[planetname, 'period']
  stareaor = planetinfo[planetname, 'stareaor']
  
  dirname = strcompress(basedir + planetname +'/')                                                            ;+'/hybrid_pmap_nn/')
  savefilename = strcompress(dirname + planetname +'_phot_ch'+chname+'_'+string(apradius)+'.sav',/remove_all) ;
  print, 'restoring ', savefilename
  restore, savefilename
  ra_ref = (planethash[aorname(0),'ra'])
  dec_ref = (planethash[aorname(0),'dec'])

  startaor =1
  stopaor = 1

;---------------
;re-run photometry at the stellar location but with a different background
  for a =startaor, stopaor do begin
     print, 'working on ',aorname(a)
     dir = dirname+ string(aorname(a) ) 
     CD, dir                    ; change directories to the correct AOR directory
     command  = strcompress( 'find ch'+chname+"/bcd -name 'SPITZER*bcd.fits' > "+dirname+'bcdlist.txt')
     spawn, command
     command2 =  strcompress('find ch'+chname+"/bcd -name '*bunc.fits' > "+dirname + 'bunclist.txt')
     spawn, command2
     command3 =  strcompress('find ch'+chname+"/raw -name '*dce.fits' > "+dirname + 'dcelist.txt')
     spawn, command3
     
     readcol,strcompress(dirname +'bcdlist.txt'),fitsname, format = 'A', /silent
     readcol,strcompress(dirname+'bunclist.txt'),buncname, format = 'A', /silent
     readcol,strcompress(dirname+'dcelist.txt'),rawname, format = 'A', /silent
     
     print,'n_elements(fitsname)', n_elements(fitsname)   
     
     xarr = fltarr(63*(n_elements(fitsname)))
     fluxarr3_7= xarr
     fluxarr5_10= xarr
     fluxarr7_15= xarr
     fluxerrarr = xarr
     corrfluxarr = xarr
     corrfluxerrarr = xarr
     backarr3_7 = xarr
     backarr5_10 = xarr
     backarr7_15 = xarr
     backerrarr = xarr
     npcentroidsarr = xarr
     xfwhmarr = xarr
     yfwhmarr = xarr
     
;---------------
;calculate the correction factor from the pmap dataset so I don't have to re-run hybrid
     corrfactor = (planethash[aorname(a),'corrflux']) / (planethash[aorname(a),'flux'])
     startfits = 0L


     for i =startfits,  n_elements(fitsname) - 1  do begin ;read each cbcd file, find centroid, keep track
        header = headfits(fitsname(i))                     ;
        frametime = sxpar(header, 'FRAMTIME')
        ch = sxpar(header, 'CHNLNUM')
        ronoise = sxpar(header, 'RONOISE') ; these are zeros
        gain = sxpar(header, 'GAIN')
        fluxconv = sxpar(header, 'FLUXCONV')
        exptime = sxpar(header, 'EXPTIME')
        aintbeg = sxpar(header, 'AINTBEG')
        atimeend = sxpar(header, 'ATIMEEND')
        naxis = sxpar(header, 'NAXIS')
        
        if ch eq '2' and frametime eq 2 then ronoise = 12.1
        
                                ;read in the files
        fits_read, fitsname(i), im, h
        fits_read, buncname(i), unc, hunc
        test = sxpar(h, 'RAWFILE')
        
                                ;apply mask file if necessary
        if planetname eq 'hd189733' then  im[13:16, 4:7, *] = !Values.F_NAN ;mask region with nan set for bad regions
        if planetname eq 'WASP-14b' then   im[4:7, 13:16, *] = !Values.F_NAN ;mask region with nan set for bad regions      
        if planetname eq 'HAT-P-22' and chname eq '2' then  im[19:25, 9:15, *] = !Values.F_NAN                              
        if planetname eq 'HAT-P-22' and chname eq '1' then  im[5:11, 11:17, *] = !Values.F_NAN                          
        if planetname eq 'HD93385' then im[17:21, 21:27, *] = !Values.F_NAN ;mask region with nan set for bad regions      
        
                                ;run the centroiding and photometry
        get_centroids_for_calstar_jk,im, h, unc, ra_ref, dec_ref,  t, dt, hjd, xft, x3, y3, $
                                     x5, y5, x7, y7, xg, yg, xh, yh, f, b, x3s, y3s, x5s, y5s, $
                                     x7s, y7s, fs, bs, xp3, yp3, xp5, yp5, xp7, yp7, xp3s, yp3s, $
                                     xp5s, yp5s, xp7s, yp7s, fp, fps, np, flag, ns, sf, $
                                     xfwhm, yfwhm, /WARM
        
        nanfound = where(FINITE(np) lt 1, nancount)
        if nancount gt 0 then print, 'NAN: ', fitsname(i), nanfound
        
        x_center = temporary(x3)
        y_center = temporary(y3)
                                ;background 3-7
        abcdflux3_7 = f[*,apval]      
        fs3_7 = fs[*,apval]                               
        back3_7 = b[*,0]
        backerr3_7 = bs[*,0]
        npcentroids = np       
                                 ;background 5-10
        abcdflux5_10= f[*,4]      
        fs5_10 = fs[*,4]                               
        back5_10 = b[*,1]
        backerr5_10 = bs[*,1]
                                ;background 3-7
        abcdflux7_15 = f[*,5]      
        fs7_15 = fs[*,5]                               
        back7_15 = b[*,2]
        backerr7_15 = bs[*,2]
       
        ;remove the first frame, store them in a big array
        fluxarr3_7[i*63] = abcdflux3_7[1:63]
        fluxarr5_10[i*63] = abcdflux5_10[1:63]
        fluxarr7_15[i*63] = abcdflux7_15[1:63]
        backarr3_7[i*63] = back3_7[1:63]
        backarr5_10[i*63] = back5_10[1:63]
        backarr7_15[i*63] = back7_15[1:63]

                                ;find the pmap correction
                                ;this won't work on the different background levels because the pmap has been measured with [3,7]
;        corrflux3_7 = pmap_correct(x_center,y_center,abcdflux3_7,$
;                                ch,xfwhm,yfwhm,NP = npcentroids,$
;                                FUNC=fs,CORR_UNC=corrfluxerr,$
;                                DATAFILE=pmapfile,NNEAREST=nn)

       endfor                   ; for each image
 
  endfor                        ; for each AOR
  
    savename = '/Users/jkrick/irac_warm/HD209458/testbkgd.sav'
     save, /all, filename = savename

end
