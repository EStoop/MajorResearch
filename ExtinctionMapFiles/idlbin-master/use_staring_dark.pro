pro use_staring_dark

  sdarkdir = '/Users/jkrick/irac_warm/darks/staring/' ; staringdark_s0p4s_ch1.fits

;setup: for each planet
  planetname ='HAT-P-8';, 'pmap_star_ch2'; 'HD209458';, 'HD158460'
  chname = '2'
 
;00000000000000000000000000000000000000000000000000000
;read in the data in fits files
  planetinfo = create_planetinfo()
  if chname eq '2'  then aorname= planetinfo[planetname, 'aorname_ch2'] else aorname = planetinfo[planetname, 'aorname_ch1'] 
  basedir = planetinfo[planetname, 'basedir']
  dirname = strcompress(basedir + planetname +'/')
  ra_ref = planetinfo[planetname, 'ra']   & dec_ref = planetinfo[planetname, 'dec']
  
  print, 'ra, dec', ra_ref, dec_ref
  planethash = hash()
  
  startaor =0
  stopaor = n_elements(aorname) - 1
  for a =startaor, stopaor do begin
     print, 'working on ',aorname(a)
     dir = dirname+ string(aorname(a) ) 
     CD, dir                    ; change directories to the correct AOR directory
     command  = strcompress( 'find ch'+chname+"/bcd -name 'SPITZER*bcd.fits' > "+dirname+'bcdlist.txt')
     print, 'command', command
     spawn, command
     command2 =  strcompress('find ch'+chname+"/bcd -name '*bunc.fits' > "+dirname + 'bunclist.txt')
     spawn, command2
     command3 =  strcompress('find ch'+chname+"/raw -name '*dce.fits' > "+dirname + 'dcelist.txt')
     spawn, command3
     
     readcol,strcompress(dirname +'bcdlist.txt'),fitsname, format = 'A', /silent
     readcol,strcompress(dirname+'bunclist.txt'),buncname, format = 'A', /silent
     readcol,strcompress(dirname+'dcelist.txt'),rawname, format = 'A', /silent
     
     print,'n_elements(fitsname)', n_elements(fitsname)
;     aparr = dblarr(n_elements(fitsname))  ;keep the aperture sizes used
     
     ;get skyflat ; only need to do this once
     skyflatname = strcompress( dirname + aorname(a) + '/ch' + chname + '/cal/*superskyflat*.fits',/remove_all)
     if a eq startaor then begin
        fits_read, skyflatname, flatdata, flatheader
                                ;need to make this [32,32,64]
        flat64 = fltarr(32,32,64)
        flatsingle = flatdata[*,*,0]
        for f = 0, 63 do flat64[*,*,f] = flatsingle
     endif

     startfits = 0.D
     stopfits =  n_elements(fitsname) - 1
     for i =startfits, stopfits  do begin ;read each cbcd file, find centroid, keep track
        header = headfits(fitsname(i)) ;
        sclk_obs= sxpar(header, 'SCLK_OBS')
        frametime = sxpar(header, 'FRAMTIME')
        bmjd_obs = sxpar(header, 'BMJD_OBS')
                                ;utcs_obs = sxpar(header, 'UTCS_OBS')
        ch = sxpar(header, 'CHNLNUM')
        ronoise = sxpar(header, 'RONOISE') ; these are zeros
        gain = sxpar(header, 'GAIN')
        fluxconv = sxpar(header, 'FLUXCONV')
        exptime = sxpar(header, 'EXPTIME')
        aintbeg = sxpar(header, 'AINTBEG')
        atimeend = sxpar(header, 'ATIMEEND')
        naxis = sxpar(header, 'NAXIS')
        framedly = sxpar(header, 'FRAMEDLY')
        if ra_ref gt 359. then begin
           print, 'inside ra_ref gt 359'
           ra_ref = sxpar(header, 'RA_RQST')
           dec_ref = sxpar(header, 'DEC_REF')
        endif
        
        if ch eq '2' and frametime eq 2 then ronoise = 12.1
        if i eq startfits then sclk_0 = sclk_obs
        
      ;read in the files
      fits_read, fitsname(i), data, h
      fits_read, buncname(i), unc, hunc

      ;apply mask file if necessary
      if planetname eq 'hd189733' then  data[13:16, 4:7, *] = !Values.F_NAN ;mask region with nan set for bad regions
      if planetname eq 'WASP-14b' then  data[4:7, 13:16, *] = !Values.F_NAN ;mask region with nan set for bad regions      
      if planetname eq 'HAT-P-22' and chname eq '2' then  data[19:25, 9:15, *] = !Values.F_NAN                              
      if planetname eq 'HAT-P-22' and chname eq '1' then  data[5:11, 11:17, *] = !Values.F_NAN            
      if planetname eq 'HD93385' then data[17:21, 21:27, *] = !Values.F_NAN ;mask region with nan set for bad regions                         

        
;00000000000000000000000000000000000000000000000000000
;remove dithered correction and put in staring mode correction

      data = data / fluxconv    ; now in DN/s
      data = data* exptime      ; now in DN
                                ; flip the image
      data = reverse(data, 2)
      
                                ;undo the flat removal
      data = data * flat64
      
                                ;remove the dark that was already used in the image
      darkname = sxpar(header, 'SKDKRKEY')
      darkepid = sxpar(header, 'SDRKEPID')
      framedelay = sxpar(header, 'FRAMEDLY')
      aorkey = sxpar(header, 'AORKEY')
      fits_read, strcompress('ch' + chname+'/cal/SPITZER_I'+chname+'_'+string(darkname)+ '_0000_*_C'+string(darkepid)+'_sdark.fits',/remove_all), darkdata, darkheader
      data = data + darkdata

      ;------------
                                ;put in the staring mode dark
      if i eq startfits then begin
         if exptime gt 0.35 and exptime lt 0.37 then staredarkname = strcompress(sdarkdir + 'staringdark_s0p4s_ch2.fits')
         if exptime gt 1.8 and exptime lt 2.0 then staredarkname = strcompress(sdarkdir + 'staringdark_s2s_ch2.fits')
         if exptime gt 0.07 and exptime lt 0.1 then staredarkname = strcompress(sdarkdir + 'staringdark_s0p1s_ch2.fits')
         if exptime lt 0.05 then staredarkname = strcompress(sdarkdir + 'staringdark_s0p02s_ch2.fits')
         fits_read, staredarkname, staredarkdata, staredarkheader
      endif

      data = data - staredarkdata
      
                                ;remove the flat back again
      data = data / flat64
      
      ;re-flip
      data = reverse(data, 2)

      ;convert back to MJy/sr
      data = data/ exptime      ; now in DN
      data = data * fluxconv    ; now in DN/s

      ;write out the fits files
      outname = strcompress(strmid(fitsname(i), 0, 43) + '_sdark.fits',/remove_all)
      fits_write, outname, data, h

   endfor                       ; for each fits file

  endfor                        ; for each AOR 


;00000000000000000000000000000000000000000000000000000
;run photometry

;00000000000000000000000000000000000000000000000000000
;save photometry in a save file

;00000000000000000000000000000000000000000000000000000
;run plot_nonlin.pro


end
