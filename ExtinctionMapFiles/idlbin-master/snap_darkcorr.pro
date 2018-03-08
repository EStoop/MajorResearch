pro snap_darkcorr, planetname, frametime, chname
;this code backs the dark correction out of the snapshot observations,
;and then puts back the same dark for all snaps.

;assumes subarray - aka naxis = 3
  

;  basedir = '/Users/jkrick/irac_warm/'+planetname + '/'
  basedir = '/Users/jkrick/irac_warm/pcrs_planets/'+planetname + '/'

 
; the superdark
  if frametime eq 2 then frametimename = 's2s'
  if frametime eq 0.1 then frametimename = 's0p1s'
  if frametime eq 0.4 then frametimename = 's0p4s'
  fits_read,strcompress( '/Users/jkrick/external/irac_warm/darks/superdarks/superdark_ch'+chname+'_'+frametimename +'_S19.fits',/remove_all), superdark, superheader

  if planetname eq 'HD158460' then aorname = [ 'r45184256','r45184512','r45184768','r45185024','r45185280','r45185536','r45185792','r45186048','r45186304','r45186560','r45186816','r45187072','r45187328','r45187584','r45187840','r45188096','r45188352','r45188608'] ;0,1s snaps ch2

  if planetname eq 'WASP-14b' then aorname = ['r45428992', 'r45428736', 'r45428480', 'r45428224', 'r45426688','r45840128','r45839616','r45839104','r45838592','r45847040','r45846784','r45846528','r45846272','r45846016','r45843200','r45842944','r45842688','r45842432','r45842176','r45841920','r45841664','r45841408','r45845504','r45840896','r45845760','r45845504','r45845248','r45844992','r45844736','r45844480','r45844224','r45843968','r45843712','r45843456','r45840640','r45840384','r45839872','r45839360','r45838848','r45838336','r48688384','r48688128','r48687872','r48687616','r48683776','r48683264','r48682752','r48682240','r48681472','r48681216','r48680704']  ; ch2 stares and snaps  s2s
  

  if planetname eq 'HD209458' then aorname =  ['r38704128','r38704384','r38701312','r38703360','r38703616'] ;'r38703872',
;['r45188864','r45189120','r45189376','r45189632','r45189888','r45190144','r45190400','r45190656','r45190912','r45191168','r45191424','r45191680','r45191936','r45192192','r45192704','r45195264','r45192960','r45193216','r45193472','r45193984','r45193728','r45195520','r45194240','r45194496','r45194752','r45195008','r45196288','r45195776','r45197312','r45196032','r45196544','r45196800','r45197056','r45197568','r45197824','r45198080','r45192448'] ; just snaps sp4s
  prefluxarr = fltarr(n_elements(aorname) * 3200* 64)  ;just guesses at size for now
  superfluxarr = prefluxarr
  c = 0L

  for a = 0,   n_elements(aorname) - 1 do begin
     print, 'working on ',aorname(a)
     dir = basedir + string(aorname(a) ) 
     print, 'dir', dir
     CD, dir                    ; change directories to the correct AOR directory
     command  = strcompress( 'find ch'+chname+"/bcd -name 'SPITZER*_bcd.fits' > "+dir+'bcdlist.txt')
     spawn, command
     command2 =  strcompress('find ch'+chname+"/bcd -name 'SPITZER*_bunc.fits' > "+dir+'bunclist.txt')
     spawn, command2
   
     readcol,strcompress(dir +'bcdlist.txt'),fitsname, format = 'A', /silent
     readcol,strcompress(dir +'bunclist.txt'),buncname, format = 'A', /silent
     print,'n_elements(fitsname)', n_elements(fitsname)
     
       ;need to read in the flat and make this [32,32,64]
     skyflatname = strcompress( 'ch' + chname + '/cal/*superskyflat*.fits',/remove_all)
     fits_read, skyflatname, flatdata, flatheader
     flat64 = fltarr(32,32,64)
     flatsingle = flatdata[*,*,0]
     for f = 0, 63 do flat64[*,*,f] = flatsingle
     skyflatout = strcompress( 'ch' + chname + '/cal/sub_skyflat.fits',/remove_all)
     fits_write, skyflatout, flat64  ; give it junk header for now

     for i = 0,  n_elements(fitsname) - 1 do begin ; for each image

                                ;read in the data
        fits_read, fitsname(i), data, header
        fits_read, buncname(i), unc, uncheader

        ra_ref = sxpar(header, 'RA_RQST')
        dec_ref = sxpar(header, 'DEC_RQST')

        ; do some baseline aperture photometry
;        get_centroids_for_calstar_jk,data, header, unc, ra_ref, dec_ref,  t, dt, hjd, xft, x3, y3, $
;                                     x5, y5, x7, y7, xg, yg, xh, yh, f, b, x3s, y3s, x5s, y5s, $
;                                     x7s, y7s, fs, bs, xp3, yp3, xp5, yp5, xp7, yp7, xp3s, yp3s, $
;                                     xp5s, yp5s, xp7s, yp7s, fp, fps, np, flag, ns, sf, $
;                                     xfwhm, yfwhm, /WARM
;        prefluxarr(c) = f[*,1]      
        
                                ;back out the flux conversion
        fluxconv = sxpar(header, 'FLUXCONV')
        exptime = sxpar(header, 'EXPTIME')
        data = data / fluxconv  ; now in DN/s
        data = data* exptime    ; now in DN
        
                                ; flip the image
        data = reverse(data, 2)
        
                                ;remove the flat
        data = data * flat64
        
                                ;remove the dark that was already used in the image
        command  = strcompress( 'find ch'+chname+"/cal -name '*_sdark.fits' > "+dir+'dark.txt')
        spawn, command
        readcol,strcompress(dir +'dark.txt'),darkname, format = 'A', /silent
        fits_read, darkname, dark, darkheader
        
        data = data + dark
       
                                ;--------------------------
        
                                ;put the same dark back into every frame
        data = data - superdark
        
                                ;divde back out the flat
        data = data / flat64
        
                                ;flip the image
        data = reverse(data, 2)
        
                                ;put the fluxconv back in
        data = data / exptime   ; DN/s
        data = data * fluxconv  ;MJy/sr
        
                                ;write out the new fits file
        outfilename =strmid(fitsname(i),0, 8) + 'sdcorr' + strmid(fitsname(i), 8)
        fits_write, outfilename, data, header
        

        ;do some aperture photometry after changing the dark
;        get_centroids_for_calstar_jk,data, header, unc, ra_ref, dec_ref,  t, dt, hjd, xft, x3, y3, $
;                                     x5, y5, x7, y7, xg, yg, xh, yh, f, b, x3s, y3s, x5s, y5s, $
;                                     x7s, y7s, fs, bs, xp3, yp3, xp5, yp5, xp7, yp7, xp3s, yp3s, $
;                                     xp5s, yp5s, xp7s, yp7s, fp, fps, np, flag, ns, sf, $
;                                     xfwhm, yfwhm, /WARM
;        superfluxarr(c) = f[*,1]      
;        c = c + 63

     endfor                     ; for each image
  endfor                        ; for each AOR

;  prefluxarr = prefluxarr[0:c-63]
;  superfluxarr = superfluxarr[0:c-63]
;  plothist, prefluxarr, prexhist, preyhist, /noplot, bin = 0.0001
;  h = plot(prexhist, preyhist, thick = 2, xtitle = 'Aperture flux', ytitle = 'Number', title = 'WASP-14b', $
;           color = 'Navy', name = 'campaign dark', xrange = [0.057, 0.0595])
  
;  plothist, superfluxarr, superxhist, superyhist, /noplot, bin = 0.0001
;  h2 =  plot(superxhist, superyhist, thick = 2, /overplot, color = 'red', name = 'superdark')

;  l = legend(target = [h, h2], position = [0.0595, 1000], /data, /auto_text_color)

end
