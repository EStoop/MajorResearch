pro read_staring_dark
  
;what is the mean level as a function of frame number and amplifier and frame delay

;  aorname = 'r40151040'        ; this one has a dark change in the middle of it.
  ;40151552 has a star in it, 
  aorname = ['/Users/jkrick/irac_warm/hd189733/r40152064/' ,'/Users/jkrick/irac_warm/hd189733/r40151040/'  , '/Users/jkrick/irac_warm/hd189733/r40150016/', '/Users/jkrick/irac_warm/hd189733/r40150528/', '/Users/jkrick/irac_warm/hd189733/r40151552/', '/Users/jkrick/irac_warm/HD209458/r44201728/', '/Users/jkrick/irac_warm/HD209458/r44201216/', '/Users/jkrick/irac_warm/HD209458/r44201472/', '/Users/jkrick/irac_warm/pcrs_planets/HD158460/r42051584/', '/Users/jkrick/irac_warm/pcrs_planets/HD158460/r42506496/']
  
  chname =[  '2',  '2', '2', '2', '2', '1',  '1', '1', '1', '1']
  
  a = findgen(32)
  b = a mod 4
  amp0 = where(b eq 0)
  amp1 = where(b eq 1)
  amp2 = where(b eq 2)
  amp3 = where(b eq 3)
  
                                ;for each aor I need to carry around
 ; fitsname, framenumber(0: n_elements(fitsname)), imagenum[0:63,0:63,...)], mean & stddev in each amplifier[0:63], framedelay[0:n_elements(fitsname)]
  biashash = hash()
 
  for a = 0, n_elements(aorname) - 1 do begin
     print, 'working on aor', aorname(a)
     skyflatname = strcompress( aorname(a) + 'ch' + chname(a) + '/cal/*superskyflat*.fits',/remove_all)
     print, 'skyflatname', skyflatname
     fits_read, skyflatname, flatdata, flatheader
                                ;need to make this [32,32,64]
     flat64 = fltarr(32,32,64)
     flatsingle = flatdata[*,*,0]
     for f = 0, 63 do flat64[*,*,f] = flatsingle
     
    
     CD, aorname(a)             ; change directories to the correct AOR directory
     command  = strcompress( 'find ch'+chname(a)+"/bcd -name 'SPITZER*_bcd.fits' > "+aorname(a)+'bcdlist.txt')
     spawn, command
     readcol,strcompress(aorname(a) +'bcdlist.txt'),fitsname, format = 'A', /silent
     
     ;arrays over all subframes in the AOR
     meanarr_0 = fltarr(n_elements(fitsname)*64)
     meanarr_1 = meanarr_0
     meanarr_2 = meanarr_0
     meanarr_3 = meanarr_0
     stddevarr_0 = meanarr_0
     stddevarr_1 = meanarr_1
     stddevarr_2 = meanarr_2
     stddevarr_3 = meanarr_3
     framedelayarr = meanarr_0
     subimagenum = meanarr_0
     medarr_0 = meanarr_0
     medarr_1  = meanarr_0
     medarr_2 = meanarr_0
     medarr_3 = meanarr_0
     ;want this to run from 0 to 63 for n_elements(fitsname)
     for si = 0, n_elements(fitsname)*63, 64 do subimagenum(si) = findgen(64)
    
     c = 0L
     
     for i = 0,  n_elements(fitsname) - 1 do begin
        fits_read, fitsname(i), data, header
        
                                ;back out the flux conversion
        fluxconv = sxpar(header, 'FLUXCONV')
        exptime = sxpar(header, 'EXPTIME')
        data = data / fluxconv  ; now in DN/s
        data = data* exptime  ; now in DN
                                ; flip the image
        data = reverse(data, 2)
        
                                ;remove the flat
        data = data * flat64
        
        darkname = sxpar(header, 'SKDKRKEY')
        darkepid = sxpar(header, 'SDRKEPID')
        framedelay = sxpar(header, 'FRAMEDLY')
        aorkey = sxpar(header, 'AORKEY')

                                ;remove the dark that was already used in the image
        fits_read, strcompress(aorname(a) + 'ch' + chname(a)+'/cal/SPITZER_I'+chname(a)+'_'+string(darkname)+ '_0000_*_C'+string(darkepid)+'_sdark.fits',/remove_all), darkdata, darkheader
        data = data + darkdata
        
        for j = 0, 63 do begin
           meanclip, data[amp0,*,j], m0, sm0
           meanarr_0(c) = m0
           stddevarr_0(c) = sm0
           medarr_0(c) = median(data[amp0,*,j])
           meanclip, data[amp1,*,j], m1, sm1
           meanarr_1(c) = m1
           stddevarr_1(c) = sm1
           medarr_1(c) = median(data[amp1,*,j])
           meanclip, data[amp2,*,j], m2, sm2
           meanarr_2(c) = m2
           stddevarr_2(c) = sm2
           medarr_2(c) = median(data[amp2,*,j])
           meanclip, data[amp3,*,j], m3, sm3
           meanarr_3(c) = m3
           stddevarr_3(c) = sm3
           medarr_3(c) = median(data[amp3,*,j])
           
           framedelayarr(c) = framedelay;keep track of  delaytime for later use
           c = c + 1
        endfor
                                

        ;plot some distributions that are going into these means
        if i eq 200 then begin
           ;just look at a few subframes
           plothist, data[amp0,*,10], xhist, yhist, /noplot,bin = 0.5
           print, 'data values image 200 frame 10', data[amp0, 10, 10]
           testplot = plot(xhist, yhist, thick = 2, title = strcompress('200th frame' + aorkey), xtitle = 'bias level per subframe', ytitle = 'Number', xrange = [0,20], yrange = [0,40])
           plothist, data[amp0,*,30], xhist, yhist, /noplot,bin = 0.5
           print, 'data values image 200 frame 30', data[amp0, 10, 30]
;           testplot = plot(xhist, yhist, thick = 2,/overplot, color = 'red')
           plothist, data[amp0,*,50], xhist, yhist, /noplot,bin = 0.5
           testplot = plot(xhist, yhist, thick = 2,/overplot, color = 'blue')
           print, 'data values image 200 frame 50', data[amp0, 10, 50]
           
        endif

        ;plot some distributions that are going into these means
        if i eq 500 then begin
           ;just look at a few subframes
           plothist, data[amp0,*,10], xhist, yhist, /noplot,bin = 0.5
           print, 'data values image 500 frame 10', data[amp0, 10, 10]
           testplot = plot(xhist, yhist, thick = 2, title =strcompress('500th frame' + aorkey), xtitle = 'bias level per subframe', ytitle = 'Number', xrange = [0,20], yrange = [0,40])
           plothist, data[amp0,*,30], xhist, yhist, /noplot,bin = 0.5
           print, 'data values image 500 frame 30', data[amp0, 10, 30]
;           testplot = plot(xhist, yhist, thick = 2,/overplot, color = 'red')
           plothist, data[amp0,*,50], xhist, yhist, /noplot,bin = 0.5
           testplot = plot(xhist, yhist, thick = 2,/overplot, color = 'blue')
           print, 'data values image 500 frame 50', data[amp0, 10, 50]

        endif

     endfor ; for each fitsname
;  print, 'n', n_elements(fitsname), n_elements(fitsname)*64, c - 1
     framnumber = findgen(n_elements(meanarr_0))
     framtime = sxpar(header, 'FRAMTIME')
     aorlabel = sxpar(header, 'AORLABEL')

     ;loooking for correlations
     AVRSTBEG = sxpar(header, 'AVRSTBEG')
     AVDETBEG = sxpar(header, 'AVDETBEG')
     AVGG1BEG  = sxpar(header, 'AVGG1BEG')
     AVDDUBEG = sxpar(header, 'AVDDUBEG')
     AHTRIBEG= sxpar(header, 'AHTRIBEG')
     AFPAT2B= sxpar(header, 'AFPAT2B')
     AFPAT2E= sxpar(header, 'AFPAT2E')
     ACTENDT= sxpar(header, 'ACTENDT')
     AFPECTE= sxpar(header, 'AFPECTE')
     AFPEATE= sxpar(header, 'AFPEATE')
     ASHTEMPE= sxpar(header, 'ASHTEMPE')
     ATCTEMPE= sxpar(header, 'ATCTEMPE')
     APDTEMPE= sxpar(header, 'APDTEMPE')
     ACATMP1E= sxpar(header, 'ACATMP1E')            
     ACATMP2E= sxpar(header, 'ACATMP2E')
     ACATMP3E= sxpar(header, 'ACATMP3E')
     ACATMP4E= sxpar(header, 'ACATMP4E')
     ACATMP5E= sxpar(header, 'ACATMP5E')
     ACATMP6E= sxpar(header, 'ACATMP6E')
     ACATMP7E= sxpar(header, 'ACATMP6E')
     ACATMP8E= sxpar(header, 'ACATMP8E')
     
                                ;remove outliers from the y range
     meansort = meanarr_0(sort(meanarr_0))
     ns= .001*n_elements(meansort)
     clipsort = meansort[ns:n_elements(meansort) - ns]
     myrange = [min(clipsort), max(clipsort)]
     
     if chname(a) eq '1' then myrange =  [115,180]
     if chname(a) eq '2' then myrange =  [9,22]
     ;plot whole time range for each amplifier
;     p = plot(framnumber, meanarr_0, '1s', xtitle = 'Frame Number', ytitle = 'Mean level in Frame (DN)', $
;              title = aorlabel+ ' ch' + chname(a)+'  '+string(framtime), $
;              sym_size = 0.1,   sym_filled = 1, yrange =myrange)
;     p1= plot(framnumber, meanarr_1, '1s', sym_size = 0.1,   sym_filled = 1, color = 'gray',/overplot)
;     p2= plot(framnumber, meanarr_2, '1s', sym_size = 0.1,   sym_filled = 1, color = 'blue',/overplot)
;     p3= plot(framnumber, meanarr_3 ,'1s', sym_size = 0.1,   sym_filled = 1, color = 'cyan',/overplot)
;     t1 = text(2E4, 170, aorlabel)

;     plotsavename = strcompress('/Users/jkrick/irac_warm/darks/staring/r'+string(aorkey)+ '_ch'+chname(a) +'_'+string(framtime)+'.png',/remove_all)
;     p.save, plotsavename
  
     keys =['aorname','aorlabel','aorkey', 'framtime', 'chname', 'fitsname', 'framenumber', 'subimagenum', 'mean_0', 'stddev_0','mean_1', 'stddev_1','mean_2', 'stddev_2','mean_3', 'stddev_3','framedelay','med_0', 'med_1', 'med_2', 'med_3','AVRSTBEG', 'AVDETBEG','AVGG1BEG','AVDDUBEG','AHTRIBEG','AFPAT2B','AFPAT2E','ACTENDT','AFPECTE','AFPEATE','ASHTEMPE','ATCTEMPE','APDTEMPE','ACATMP1E','ACATMP2E','ACATMP3E','ACATMP4E','ACATMP5E','ACATMP6E','ACATMP7E','ACATMP8E']
     values=list(aorname(a),aorlabel, aorkey, framtime,chname(a), fitsname, findgen(n_elements(fitsname)), subimagenum, meanarr_0, stddevarr_0, meanarr_1, stddevarr_1,meanarr_2, stddevarr_2,meanarr_3, stddevarr_3,framedelayarr, medarr_0, medarr_1, medarr_2, medarr_3,AVRSTBEG, AVDETBEG,AVGG1BEG,AVDDUBEG,AHTRIBEG,AFPAT2B,AFPAT2E,ACTENDT,AFPECTE,AFPEATE,ASHTEMPE,ATCTEMPE,APDTEMPE,ACATMP1E,ACATMP2E,ACATMP3E,ACATMP4E,ACATMP5E,ACATMP6E,ACATMP7E,ACATMP8E)
;     biashash[aorname(a)] = HASH(keys, values)

     print, 'end for each aorname'
  endfor                        ; for each aorname
  
;  save, biashash, filename='/Users/jkrick/irac_warm/darks/staring/staring_darks.sav'
end
