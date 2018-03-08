pro plot_staring_darks, ch,ignoreearly = ignoreearly

restore,  '/Users/jkrick/irac_warm/darks/staring/staring_darks.sav'
 aorname = ['/Users/jkrick/irac_warm/hd189733/r40152064/' ,'/Users/jkrick/irac_warm/hd189733/r40151040/'  , '/Users/jkrick/irac_warm/hd189733/r40150016/', '/Users/jkrick/irac_warm/hd189733/r40150528/', '/Users/jkrick/irac_warm/hd189733/r40151552/', '/Users/jkrick/irac_warm/HD209458/r44201728/', '/Users/jkrick/irac_warm/HD209458/r44201216/', '/Users/jkrick/irac_warm/HD209458/r44201472/', '/Users/jkrick/irac_warm/pcrs_planets/HD158460/r42051584/', '/Users/jkrick/irac_warm/pcrs_planets/HD158460/r42506496/']
 
colorarr= ['blue', 'red', 'purple', 'orange', 'gray','blue', 'red', 'purple', 'orange', 'gray']
;for reference 
;keys =['aorname','aorlabel', 'framtime', 'chname', 'fitsname', 'framenumber', 'subimagenum', 'mean_0', 'stddev_0','mean_1', 'stddev_1','mean_2', 'stddev_2','mean_3', 'stddev_3','framedelay']
if ch eq 1 then begin
   startaor = 5
   stopaor = n_elements(aorname) - 1
   meanyrange = [116, 130]
   yr = [120,134]
   xr = [0,1300]
endif
if ch eq 2 then begin
   startaor = 0
   stopaor = 4      
   meanyrange = [9.5,11.5]
   yr = [9.5, 12]
   xr = [0,650]
endif

na = n_elements(aorname)
AVRSTBEGarr = fltarr(na)
AVDETBEGarr  = avrstbegarr
AVGG1BEGarr  = avrstbegarr
AVDDUBEGarr  = avrstbegarr
AHTRIBEGarr  = avrstbegarr
AFPAT2Barr = avrstbegarr
AFPAT2Earr  = avrstbegarr
ACTENDTarr  = avrstbegarr
AFPECTEarr  = avrstbegarr
AFPEATEarr  = avrstbegarr
ASHTEMPEarr  = avrstbegarr
ATCTEMPEarr  = avrstbegarr
APDTEMPEarr  = avrstbegarr
ACATMP1Earr  = avrstbegarr
ACATMP2Earr  = avrstbegarr
ACATMP3Earr  = avrstbegarr
ACATMP4Earr  = avrstbegarr
ACATMP5Earr  = avrstbegarr
ACATMP6Earr = avrstbegarr
ACATMP7Earr  = avrstbegarr
ACATMP8Earr  = avrstbegarr

levelarr = fltarr(na)

 for a = startaor,  stopaor do begin
;----------------------------------------------------------------------------
;----------------------------------------------------------------------------
;------------------------
; would like to plot mean with error bars as a function of time color coded by frame delay for one amplifier

    framnumber = findgen(n_elements(biashash[aorname(a), 'mean_0']))
;Take out the first three images, mean together the rest of the 64 subimages 
    rmnum = 2.

    if keyword_set(ignoreearly) then begin
       f1 = where(biashash[aorname(a),'framedelay'] ge 1.5 and biashash[aorname(a),'framedelay'] lt 1.81 and biashash[aorname(a), 'subimagenum'] gt rmnum, f1count)
       f2 = where(biashash[aorname(a),'framedelay'] ge 1.81 and biashash[aorname(a),'framedelay'] lt 1.9 and biashash[aorname(a), 'subimagenum'] gt rmnum, f2count)
    endif else begin
       f1 = where(biashash[aorname(a),'framedelay'] ge 1.5 and biashash[aorname(a),'framedelay'] lt 1.81 , f1count)
       f2 = where(biashash[aorname(a),'framedelay'] ge 1.81 and biashash[aorname(a),'framedelay'] lt 1.9 , f2count)
    endelse

;    f3 = where(biashash[aorname(a),'framedelay'] ge 1.7 and biashash[aorname(a),'framedelay'] lt 1.8, f3count)
;    f4 = where(biashash[aorname(a),'framedelay'] ge 1.82 and biashash[aorname(a),'framedelay'] lt 1.9, f4count)

    print, 'counts, ', f1count/(64.-rmnum - 1. ), f2count/(64.-rmnum - 1. );, f3count/64, f4count/64;, f5count/64, f6count/64

    mean0_f1 = (biashash[aorname(a), 'med_1'])(f1)
    mean0_f2 = (biashash[aorname(a), 'med_1'])(f2)
;    mean0_f3 = (biashash[aorname(a), 'mean_0'])(f3)
;    mean0_f4 = (biashash[aorname(a), 'mean_0'])(f4)
    stddev0_f1 = (biashash[aorname(a), 'stddev_1'])(f1)
    stddev0_f2 = (biashash[aorname(a), 'stddev_1'])(f2)
;    stddev0_f3 = (biashash[aorname(a), 'stddev_0'])(f3)
;    stddev0_f4 = (biashash[aorname(a), 'stddev_0'])(f4)
    subimage_f1 = (biashash[aorname(a), 'subimagenum'])(f1)
    subimage_f2 = (biashash[aorname(a), 'subimagenum'])(f2)
;    subimage_f3 = (biashash[aorname(a), 'subimagenum'])(f3)
;    subimage_f4 = (biashash[aorname(a), 'subimagenum'])(f4)
    framnumber_f1 = framnumber(f1)
    framnumber_f2 = framnumber(f2)
;    framnumber_f3 = framnumber(f3)
;    framnumber_f4 = framnumber(f4)
    

    if keyword_set(ignoreearly) then begin
       count = 0
       bin_framnumber_f1 = fltarr(f1count/(64 - rmnum - 1),/nozero)
       for j = 0, n_elements(framnumber_f1) - 1, (64 - rmnum - 1)  do begin
          bin_framnumber_f1[count]= framnumber_f1(j) /64 
          count = count + 1
       endfor
       
       count = 0
       bin_framnumber_f2 = fltarr(f2count/(64 - rmnum - 1),/nozero)
       for j = 0, n_elements(framnumber_f2) - 1, (64 - rmnum - 1)  do begin
          bin_framnumber_f2[count]= framnumber_f2(j) /64 
          count = count + 1
       endfor
       
       
       
       h = histogram(framnumber_f1, OMIN=om, binsize = 64 , reverse_indices = ri)
       
       bin_mean0_f1 = dblarr(n_elements(h))
       bin_stddev0_f1 = bin_mean0_f1
       c = 0
       for j = 0L, n_elements(h) - 1 do begin
          
          if (ri[j+1] gt ri[j] + 2)  then begin ;require 3 elements in the bin
             n = where(finite( total(mean0_f1[ri[ri[j]:ri[j+1]-1]])) lt 1, ncount)
             if ncount gt 0 then print, 'nan!'
             meanclip, mean0_f1[ri[ri[j]:ri[j+1]-1]], mean, sigma
             bin_mean0_f1[c] = mean
             bin_stddev0_f1[c] = sigma
             c = c + 1
          endif
       endfor
       
       bin_mean0_f1 = bin_mean0_f1[0:c-1]
       bin_stddev0_f1 = bin_stddev0_f1[0:c-1]
     print, 'n_ele', n_elements(bin_framnumber_f1), n_elements(bin_mean0_f1), n_elements(bin_stddev0_f1)
      pe = errorplot(bin_framnumber_f1, bin_mean0_f1, bin_stddev0_f1, '1s', sym_size = 0.5, $
                     sym_filled = 1, xtitle = 'Image number', ytitle = 'Mean Level (DN)', xrange = xr, yrange =yr,$ ;
                     title =   biashash[aorname(a), 'aorlabel']+ ' ch' + biashash[aorname(a), 'chname']+'  '+string(biashash[aorname(a),'framtime']))
;-------
       h = histogram(framnumber_f2, OMIN=om, binsize = 64 , reverse_indices = ri)
;     print, 'omin', om, 'nh', n_elements(h)
       
       bin_mean0_f2 = dblarr(n_elements(h))
       bin_stddev0_f2 = bin_mean0_f2
       c = 0
       for j = 0L, n_elements(h) - 1 do begin
          
          if (ri[j+1] gt ri[j] + 2)  then begin ;require 3 elements in the bin
             
             meanclip, mean0_f2[ri[ri[j]:ri[j+1]-1]], mean, sigma
             bin_mean0_f2[c] = mean
;             sigma = stddev(mean0_f2[ri[ri[j]:ri[j+1]-1]])
             bin_stddev0_f2[c] = sigma
             c = c + 1
          endif
       endfor
       
       bin_mean0_f2 = bin_mean0_f2[0:c-1]
       bin_stddev0_f2 = bin_stddev0_f2[0:c-1]
;  ;   print, 'n_ele', n_elements(bin_framnumber_f2), n_elements(bin_mean0_f2), n_elements(bin_stddev0_f2)
      pe = errorplot(bin_framnumber_f2, bin_mean0_f2, bin_stddev0_f2, '1s', sym_size = 0.5, $
                    sym_filled = 1,/overplot, color = 'red', errorbar_color = 'blue')
       
      t = text(450, 11.5, 'framedly 1.82 - 1.9', color = 'red',/overplot,/data)
      t = text(450, 11.6, 'framedly 1.5 - 1.82', color = 'black',/overplot,/data)
       
       
     saven = strcompress('/Users/jkrick/irac_warm/darks/staring/framedelay' + string(biashash[aorname(a),'aorkey'])+ ' ch' + biashash[aorname(a),'chname']+'  '+string(biashash[aorname(a),'framtime']) + '_amp0.png',/remove_all)
      pe.save, saven
       
    endif  ; keyword set ignoreearly

;------------------------
;want to look at the per amplifier output for just one framedelay
    meansub = fltarr(64)
    sigsub = meansub
                                ;mean within each AOR on subimagenumber.
                                ;note the first three subimages are already removed, so those will get zero means
    print, 'numbers', n_elements(subimage_f1), n_elements(mean0_f1)
    for si =0, 63 do begin
       good = where(subimage_f1 eq si, goodcount)
;       print, si, goodcount
       if goodcount gt 0 then begin
          meanclip, mean0_f1(good), mc, sc
;          if a eq 1 and si gt 57 and si lt 62 then begin
;             plothist, mean0_f1(good), xhist, yhist, /noplot, bin = 0.1
;             p58 = plot(xhist, yhist, xtitle = ' mean level, single amplifier, frame delay 1.5 - 1.81', ytitle = 'Number'$
;                        ,/overplot, thick = 2, color = colorarr(si - 57), title ='distribution of mean levels in the 58 - 61 frames')
;          endif
       endif else begin
          mc = 0
          sc = 1000.
       endelse

       meansub(si) = mc
       sigsub(si) = sc
    endfor
    
    titlename = string(biashash[aorname(a),'aorkey'])+ ' ch' + biashash[aorname(a),'chname']+'  '+string(biashash[aorname(a),'framtime'])
    print, 'titlename', titlename
;    ms = errorplot(indgen(64), meansub,sigsub, '1s', xtitle = 'Sub image Number', ytitle = 'Mean level in Frame (DN)', $
;             title =titlename, sym_size = 0.3,   sym_filled = 1, xrange = [0,65], yrange = meanyrange,/overplot); yrange = [10, 11.5])

    ; fit this with a  function
    start = [118.,2.0, 4.0];, 0.5, 0.5, 0.5]
    if a gt startaor then start = result  ; once the first AOR is done, start with the results from that run
    x = indgen(64)
    x = x[3:*]
    y = meansub[3:*]
    err =sigsub[3:*]
    pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},n_elements(start))
    ;fix the parameters to that found in the first AOR, except the levels
    if a gt startaor then begin
       pi(1).fixed = 1
       pi(2).fixed = 1
;       pi(3).fixed = 1
    endif

    result= MPFITFUN('x_1plusx',x,y, err, start,/nan, parinfo = pi,/quiet) 
;    ms2 = plot(x,result(0) +result(1)*x/(result(2) + x), thick = 3,color = 'blue',/overplot)
    resid = y - (result(0) +result(1)*x/(result(2) + x))
;    ms3 = plot(x, resid, xtitle = 'Sub Image Number', ytitle = 'Residuals',/overplot, xrange = [0, 65])
    meansub = fltarr(64)
    sigsub = meansub
    levelarr(a) = result(0)
                                ;mean within each AOR on subimagenumber.
                                ;note the first three subimages are already removed, so those will get zero means
    for si =0, 63 do begin
       good = where(subimage_f2 eq si, goodcount)
       if goodcount gt 0 then begin
          meanclip, mean0_f2(good), mc, sc
       endif else begin
          mc = 0
          sc = 0
       endelse

       meansub(si) = mc
       sigsub(si) = sc
    endfor
    
;   ms = plot(indgen(64), meansub, '1s', sym_size = 0.3, sym_filled = 1, color = 'red',/overplot)


;    t = text(40, 120, 'framedly 1.81- 1.9', color = 'red',/overplot,/data)
;    t = text(40, 121, 'framedly 1.5 - 1.81', color = 'black',/overplot,/data)
;    t = text(40, 122, 'framedly 1.7 - 1.8', color = 'slate grey',/overplot,/data)
;    t = text(40, 123, 'framedly 1.8 - 1.9', color = 'red',/overplot,/data)
    
    saven = strcompress('/Users/jkrick/irac_warm/darks/staring/framedelay' + string(biashash[aorname(a),'aorkey'])+ ' ch' + biashash[aorname(a),'chname']+'  '+string(biashash[aorname(a),'framtime']) + '.png',/remove_all)
;    ms.save, saven


;--------------------------------------------------------
     ;plot whole time range for each amplifier
                                 ;remove outliers from the y range
     meansort = meansub(sort(meansub))
     ns= .0005*n_elements(meansort)
     clipsort = meansort[ns:n_elements(meansort) - ns]
     myrange = [min(clipsort), max(clipsort)]
     framnumber = findgen(n_elements(biashash[aorname(a), 'mean_0']))


;    p = plot(framnumber, biashash[aorname(a), 'mean_0'], '1s', xtitle = 'Frame Number', ytitle = 'Mean level in Frame (DN)', $
;              title = biashash[aorname(a), 'aorlabel']+ ' ch' + biashash[aorname(a), 'chname']+'  '+string(biashash[aorname(a),'framtime']), $
;              sym_size = 0.1,   sym_filled = 1, yrange = [8,21]);, yrange = [115, 180])
;     p1= plot(framnumber, biashash[aorname(a), 'mean_1'], '1s', sym_size = 0.1,   sym_filled = 1, color = 'gray',/overplot)
;     p2= plot(framnumber, biashash[aorname(a), 'mean_2'], '1s', sym_size = 0.1,   sym_filled = 1, color = 'dark slate grey',/overplot)
;     p3= plot(framnumber, biashash[aorname(a), 'mean_3'] ,'1s', sym_size = 0.1,   sym_filled = 1, color = 'light cyan',/overplot)

;     plotsavename = strcompress('/Users/jkrick/irac_warm/darks/staring/r'+string(biashash[aorname(a),'aorkey'])+ '_ch'+biashash[aorname(a),'chname'] +'_'+string(biashash[aorname(a),'framtime'])+'.png',/remove_all)
;     p.save, plotsavename
 
;================================================
endfor


;make arrays with the bias and temp values from the headers
;plot against levelarr = bias level in the images.

end
