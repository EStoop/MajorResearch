pro track_calstar_photometry, savefilename,ch, make_plot = make_plot, binning = binning
;;track_calstar_photometry,'/Users/jkrick/irac_warm/calstars/track_calstar_ch2.sav',2, /make_plot,/binning
;;  track_calstar_photometry,'/Users/jkrick/irac_warm/calstars/track_calstar_ch1.sav', 1, /make_plot,/binning
 st1 = systime(1)

                                ;am I using the right corrections for
                                ;the right processing versions?  no,
                                ;too many processing versions.

COMMON track_block, photcordata, photcorhead, xarr,  yarr,  starnamearr,  timearr,  fluxarr,  fluxerrarr,  backarr,  corrfluxarr ,  raarr ,  decarr, bmjdarr, procverarr, campaignarr

;look for a save file with photometry
;for array dependant photometric correction warm
if ch eq 1 then fits_read, '/Users/jkrick/irac_warm/calstars/arrayloccorr/ch1_photcorr_ap_5.fits', photcordata, photcorhead
if ch eq 2 then fits_read, '/Users/jkrick/irac_warm/calstars/arrayloccorr/ch2_photcorr_ap_5.fits', photcordata, photcorhead
 
   
result = file_test(savefilename)
orig_num = 33500                ;33500  38500

if result eq 1 then begin  ;file exists
;if that file exists then figure out which campaign was last done and
;incrememnt by 1.  + turn into string
;                        to find the right directory
   restore, savefilename
   print, 'ti from file', ti
    
   ;;ti = ti - 3300                ; 100;  2300          ; 100 set back since the first thing the while loop does is to increment
endif else begin
;if that file doesn't exist, then start from old campaign
;and increment until there are no more directories
   ti = orig_num
endelse
; start from old campaign
;and increment until there are no more directories
;   ti = orig_num

;error check on existance of directory
good_dir = 1
;once I successfully make the save file with all S19.2 calstars up to
;June 2016, can just start there and move forward.
;ti = 42600
print, 'ti to start', ti
while good_dir eq 1 do begin    ; while there is data to work with, and not an empty directory
;;   if (ti ) eq 41000 then ti = 41100
   if (ti ) eq 42700 then ti = 43100
;;   endif else begin
;;      ti =ti + 100
;;   endelse
   
   print, 'ti after changing ti', ti
   dirname = strcompress('/Volumes/iracdata/flight/IWIC/IRAC0' +string(ti),/remove_all)
   print, 'dirname ', dirname
   good_dir = chk_dir(dirname)
   if good_dir eq 1 then begin  ;got a live one
;      print, 'starting photometry', systime(1) - st1
      junk = do_calstar_photometry(ch, dirname,'iracdata')
;      print, 'finished photometry', systime(1) - st1

      if ti eq orig_num then begin ;;+ 100
      ;need to dynamically assaign these arrays since I don't
      ;know how many there are
         
         bigxcen = xarr
         bigycen = yarr
         bigstarnamearr = starnamearr
         bigtimearr = bmjdarr+ 2400000.5D
         bigfluxarr = fluxarr
         bigfluxerrarr = fluxerrarr
         bigbackarr = backarr
         bigcorrfluxarr = corrfluxarr
         bigraarr = raarr
         bigdecarr = decarr
      endif else begin
         bigxcen = [bigxcen, xarr]
         bigycen = [bigycen,yarr]
         bigstarnamearr = [bigstarnamearr,starnamearr]
         bigtimearr = [bigtimearr,bmjdarr + 2400000.5D]
         bigfluxarr = [bigfluxarr, fluxarr]
         bigfluxerrarr = [bigfluxerrarr,fluxerrarr]
         bigbackarr = [bigbackarr,backarr]
         bigcorrfluxarr = [bigcorrfluxarr, corrfluxarr]
         bigraarr = [bigraarr, raarr]
         bigdecarr = [bigdecarr, decarr]
      endelse
   ti = ti + 100
   print, 'ti at end of while', ti
   junk12 = where(bigtimearr lt 2400001, wcount)
  ;; print, 'bigtimearr', n_elements(bigtimearr), n_elements(bmjdarr), wcount
   endif
endwhile

;convert bmjd to month, day, year
;bigtimearr = bigtimearr + 2400000.5 ; back to JD
;caldat, bigtimearr, montharr, dayarr, yeararr

print, 'after while loop ', ti
print, 'saving', savefilename
save,  ti, campaignarr, bigxcen,  bigycen,  bigstarnamearr,  bigtimearr,  bigfluxarr,  bigfluxerrarr,  bigcorrfluxarr, bigbackarr,  bigraarr,  bigdecarr , filename =savefilename

;;think about sorting
;;s1 = sort(bigtimearr)
;;bigtimearr = bigtimearr(s1)
;;bigcorrfluxarr = bigcorrfluxarr(s1)
;;bigfluxerrarr = bigfluxerrarr(s1)
;;bigfluxarr = bigfluxarr(s1)
;;print, min(bigtimearr), 'min(bigtimearr)', bigtimearr[0:10]
w = where(bigtimearr lt 2400001, wcount)
print, 'wcount ', wcount, n_elements(bigtimearr)
print, bigstarnamearr(w)

;for bsn =0, 500 do print, bsn, bigstarnamearr(bsn)

if keyword_set(make_plot) then begin
    colorarr = ['blue', 'red','black','green','grey','purple', 'deep_pink','fuchsia', 'medium_purple','medium_orchid', 'orchid', 'violet', 'plum', 'thistle', 'pink', 'light_pink', 'rosy_brown','pale_violet_red',  'chocolate', 'saddle_brown', 'maroon', 'hot_pink', 'dark_orange', 'peach_puff', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'dodger_blue', 'dark_blue', 'indigo','dark_slate_blue', 'blue_violet', 'purple','dim_grey', 'slate_grey', 'dark_slate_grey', 'khaki', 'tomato', 'lavender','gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green','blue', 'red','deep_pink', 'magenta', 'medium_purple','light_sea_green', 'teal']

    
colornames = ['blue', 'red','black','green','grey','purple','deep pink', 'thistle', 'indigo', 'blue_violet', 'light_sky_blue','dodger_blue', 'Navy', 'turquoise', 'pale_green', 'lime', 'orange', 'yellow', 'dark_khaki', 'sienna','lemon_chiffon']
     
   ;need to sort by star
     names = [  'KF08T3_', 'KF06T2_', 'KF06T1_', 'KF09T1_',  'NPM1p67','NPM1p60','BP20_417', 'HD4182_', 'HD37725', 'HD55677','HD77823', 'HD89562', 'HD99253', 'HD11374', 'HD13176', 'HD13742', 'HD15689', 'HD18483', 'HD19506', 'HD21852', 'HD28437'];, '1812095','NPM1p68'] ;, 'HD16545']

;;     loadct, 42, ncolors = n_elements(names), RGB_TABLE = colornames
     normvals = fltarr(n_elements(names))
     for n = 0,  n_elements(names) - 1 do begin
        
        an = where(bigstarnamearr eq names(n), count)
        if count gt 10 then begin
           print, 'n, name, color', n, names(n), colornames(n), n_elements(bigtimearr(an))
           
           ;;if n eq 0 then print, 'bigtimearr', bigtimearr(an)
                                              
           ;check normalization
           junkcorr = bigcorrfluxarr(an)
           junkflux = bigfluxarr(an)
           junkname = bigstarnamearr(an)
           print, 'stddev flux', colornames(n), stddev(junkflux)
           print, 'stddev corrflux', colornames(n), stddev(junkcorr)
;           print, 'bigstarname', junkname[0:50]
;           print, 'bigcorrfluxarr', junkcorr[0:50] 
;           print, 'bigfluxarr', junkflux[0:50]
           print, median(bigcorrfluxarr(an)), median(bigfluxarr(an))
           normcorr = median(bigcorrfluxarr(an))  ; median(bigcorrfluxarr(an[0:10]))
           normflux = median(bigfluxarr(an)); median(bigfluxarr(an[0:10]))
           p = errorplot(bigtimearr(an) , bigcorrfluxarr(an)/normcorr, $ ;- bigtimearr(0))/60./60./ 24.
                         bigfluxerrarr(an)/normcorr, '1s', sym_size = 0.5, ERRORBAR_COLOR = colornames(n),$
                         sym_filled = 1, ytitle = 'Corrected Flux',color = colornames(n), overplot = p,$
                         yrange = [0.85,1.15], XTICKFORMAT='(C(CDI,1x,CMoA,1x,CYI))', xtickunits = 'months', $
                         title = 'Ch' + string(ch), xrange = [min(bigtimearr),max(bigtimearr)],/buffer)

           pz = errorplot(bigtimearr(an) , bigfluxarr(an)/normflux, $
                         bigfluxerrarr(an)/normflux, '1s', sym_size = 0.5, ERRORBAR_COLOR = colornames(n),$
                         sym_filled = 1, ytitle = 'Flux',color = colornames(n), overplot = pz,$
                         yrange = [0.85,1.15], XTICKFORMAT='(C(CDI,1x,CMoA,1x,CYI))', xtickunits = 'months', $
                          title = 'Ch' + string(ch), xrange = [min(bigtimearr),max(bigtimearr)],/buffer)

          ;-----------------

           if keyword_set(binning) then begin
              
              s = sort(bigtimearr(an))
              sort_time = bigtimearr(an(s))
              sort_corrflux = bigcorrfluxarr(an(s))
              sort_fluxarr = bigfluxarr(an(s))
              sort_fluxerrarr = bigfluxerrarr(an(s))

              ;;print, 'sort time', sort_time
              
              h = histogram(sort_time, OMIN=om, binsize = 1.0, reverse_indices = ri)
              bin_corrflux = findgen(n_elements(h))
              bin_fluxarr = bin_corrflux
              bin_fluxerrarr = bin_corrflux
              bin_time = bin_corrflux
              
              c = 0
              for j = 0L, n_elements(h) - 1 do begin
                 if (ri[j+1] gt ri[j] + 2)  then begin ;require 3 elements in the bin
                    
                    meanclip, sort_corrflux[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
                    bin_corrflux[c] = meanx ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
                    
                    meanclip, sort_fluxarr[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
                    bin_fluxarr[c] = meanx ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
                    
                    idataerr = sort_fluxerrarr[ri[ri[j]:ri[j+1]-1]]
                    bin_fluxerrarr[c] =   sigmax/sqrt(n_elements(idataerr)) ; sqrt(total(idataerr^2))/ (n_elements(idataerr))
                    
                    meanclip, sort_time[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
                    bin_time[c] = meanx ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
                    
                    c = c + 1
                 endif
              endfor
              
              bin_corrflux = bin_corrflux[0:c-1]
              bin_fluxarr = bin_fluxarr[0:c-1]
              bin_fluxerrarr = bin_fluxerrarr[0:c-1]
              bin_time = bin_time[0:c-1]

              if n eq 19 then print, bin_time
              
              print, 'n_elements bin_corrflux', n_elements(bin_corrflux)
;              pb = errorplot(bin_time , bin_corrflux/median(bin_corrflux[0:18]), $
;                             bin_fluxerrarr/median(bin_corrflux), '1s', sym_size = 0.5, ERRORBAR_COLOR = colornames(n),$
;                             sym_filled = 1, ytitle = 'Binned Corrected Flux',color = colornames(n), overplot = pb,$
;                             yrange = [0.95,1.05], xrange = [min(bigtimearr),max(bigtimearr) + 100.], xtickunits = 'months', $ ;
;                             XTICKFORMAT='(C(CDI,1x,CMoA,1x,CYI))', title = 'Ch' + string(ch), xmajor = 7 ) ;,- bigtimearr(0))/60./60./ 24.
              dummy = LABEL_DATE(DATE_FORMAT=['%D-%M-%Y'])
              print, min(bin_time), 'min(bin_time)'
              if n_elements(bin_corrflux) lt 30 then nend = n_elements(bin_corrflux) - 1 else nend = 29
              print,'testing nend and n_elements(bin_corrflux): ', nend, n_elements(bin_corrflux)
              print, 'testing more', n_elements(bin_time), n_elements(bin_corrflux), n_elements(bin_fluxerrarr)
              pb = errorplot(bin_time , bin_corrflux/median(bin_corrflux[0:nend]), $
                             bin_fluxerrarr/median(bin_corrflux), '1s', sym_size = 0.5, ERRORBAR_COLOR = colornames(n),$
                             sym_filled = 1, ytitle = 'Binned Corrected Flux',color = colornames(n), overplot = pb,$
                             yrange = [0.95,1.05], xrange = [2456293, max(bin_time) + 300.],  XTICKUNITS = ['Time'], $ ;xrange = [min(bin_time),max(bin_time) + 300.], 2.45714e+06
                             XTICKFORMAT='LABEL_DATE', title = 'Ch' + string(ch),  xminor = 11, errorbar_capsize = 0.1, /buffer) ;,  xtickinterval = 1, xmajor = 7,- bigtimearr(0))/60./60./ 24.

;;              pb = plot(bin_time , bin_corrflux/median(bin_corrflux[0:nend]), $
;;                        '1s', sym_size = 0.5, ERRORBAR_COLOR = colornames(n),$
;;                        sym_filled = 1, ytitle = 'Binned Corrected Flux Density',color = colornames(n), overplot = pb,$
;;                        yrange = [0.95,1.05], xrange = [2456293, max(bin_time) + 300.],  XTICKUNITS = ['Time'], $ ;xrange = [min(bin_time),max(bin_time) + 300.], 2.45714e+06
;;                        XTICKFORMAT='LABEL_DATE', title = 'Ch' + string(ch),  xminor = 11) ;,  xtickinterval = 1, xmajor = 7,- bigtimearr(0))/60./60./ 24.

               
              pbl = plot(bin_time, intarr(n_elements(bin_time)) + 1.0, overplot = pb)

                            ;;add lines to the plot for May 1 2015 and Sep 17 2015
              JDCNV, 2015, 5, 1,0., Julian
              x = fltarr(5) + Julian
              y = findgen(5)
;;              pbll = plot(x, y, overplot = pb, thick = 2)
;;              t = text(Julian+5., 0.955, 'May 2015',orientation = 50, /data, overplot = pb)

              JDCNV, 2015, 9, 23,0., Julian
              ;;print, 'S19.2 Julian', julian
              x = fltarr(5) + Julian
              y = findgen(5)
;;              pbll = plot(x, y, overplot = pb, thick = 2)
;;              t = text(Julian- 55., 0.960, 'S19.2',orientation = 50,/data, overplot = pb)


              JDCNV, 2015, 10, 2,0., Julian
              ;;print, 'S19.2 Julian', julian
              x = fltarr(5) + Julian
              y = findgen(5)
;;              pbll = plot(x, y, overplot = pb, thick = 2)
;;              t = text(Julian+ 5., 0.955, 'PC162',orientation = 50,/data, overplot = pb)

               JDCNV, 2015, 10, 15,0., Julian
              ;;print, 'S19.2 Julian', julian
              x = fltarr(5) + Julian
              y = findgen(5)
;;              pbll = plot(x, y, overplot = pb, thick = 2)
;;              t = text(Julian+ 5., 0.955, 'PC163',orientation = 50,/data, overplot = pb)
              
              JDCNV, 2015, 10, 29,0., Julian
              ;;print, 'S19.2 Julian', julian
              x = fltarr(5) + Julian
              y = findgen(5)
;;              pbll = plot(x, y, overplot = pb, thick = 2)
;;              t = text(Julian+ 5., 0.955, 'PC164',orientation = 50,/data, overplot = pb)
             
              JDCNV, 2015, 12, 21,0., Julian
              x = fltarr(5) + Julian
              y = findgen(5)
              pbll = plot(x, y, overplot = pb, thick = 2)
              t = text(Julian+ 5., 0.960, 'Dec 21',orientation =50,/data, overplot = pb)
              t = text(Julian+ 10., 0.955, 'turn on',orientation =50,/data, overplot = pb)

           endif   ;if keyword set binning
           
        endif
        
        
     endfor
     
     
  endif
basedir = '/Volumes/IRAC/Calibration/Trending/'
plotname = strcompress(basedir + 'ch' + string(ch) + '_track_binned.png',/remove_all)
if keyword_set(binning) then pb.save, plotname
plotname = strcompress(basedir + 'ch' + string(ch) + '_track_corrected.png',/remove_all)
p.save, plotname
plotname = strcompress(basedir + 'ch' + string(ch) + '_track_flux.png',/remove_all)
pz.save, plotname
print, 'time check', systime(1) - st1


;need these for running this as a cron job so that the script finishes
;remove for actually looking at the plots
;if keyword_set(make_plot) then begin
;   p.close
;   pz.close
;   if keyword_set(binning) then pb.close
;endif

end




