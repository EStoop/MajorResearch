pro plot_hd209458_allch_v2

restore, '/Users/jkrick/irac_warm/hd209458/hd209458_plot_ch1.sav'
;normalize the flux
y = bin_sub+add
y1 = y / mean(y)
x1 = x

;junk = plot(findgen(n_elements(y1)), y1,'1r')
restore, '/Users/jkrick/irac_warm/hd209458/hd209458_plot_ch1_bad.sav'
;normalize the flux
yb = bin_sub+add
yb1 = yb / mean(yb)
xb1=x
;junk = plot(findgen(n_elements(yb1)), yb1,'1b',/overplot)


finalx1 = [x1[0:1800], xb1[1800:2600],x1[1801:9800],xb1[10100:10800],x1[9801:17400],xb1[18400:19000], x1[17401:*]]
finaly1 =[y1[0:1800], yb1[1800:2600],y1[1801:9800],yb1[10100:10800]-0.001,y1[9801:17400],yb1[18400:19000]-0.00005, y1[17401:*]]
 p1 = plot(finalx1, finaly1,  '1.',  sym_filled = 1.,color = 'black',name = 'Ch1',xtitle = 'Orbital Phase', ytitle = ' Normalized Flux', title = 'HD209458', yrange = [0.985, 1.01], xrange = [-0.7,0.7] ) ;, sym_size = 0.1, xrange = [2,4]


;-----
restore, '/Users/jkrick/irac_warm/hd209458/hd209458_plot_ch4.sav'
;normalize the flux
y = bin_sub+add
y4 = y / mean(y)
x4=x

;junk = plot(findgen(n_elements(y4)), y4,'1r')
restore, '/Users/jkrick/irac_warm/hd209458/hd209458_plot_ch4_bad.sav'
;normalize the flux
yb = bin_sub+add
yb4 = yb / mean(yb)
xb4=x
;junk = plot(findgen(n_elements(yb4)), yb4,'1b',/overplot)

finalx4 = [x4[0:180], xb4[180:400],x4[181:*],xb4[3100:*]]
finaly4 =[y4[0:180], yb4[180:400]-0.00005,y4[181:*],yb4[3100:*]-0.0005]

p4 = plot(finalx4, finaly4,  '1.b', sym_filled = 1.,/overplot, name = 'Ch4') ;,  sym_size = 0.1

;---------
restore, '/Users/jkrick/irac_warm/hd209458/hd209458_plot_ch2.sav'
;normalize the flux
y = bin_sub+add
y2 = y / mean(y)
x2 = x

;junk = plot(findgen(n_elements(y2)), y2,'1r')
restore, '/Users/jkrick/irac_warm/hd209458/hd209458_plot_ch2_bad.sav'
;normalize the flux
yb = bin_sub+add
yb2 = yb / mean(yb)
xb2=x
;junk = plot(findgen(n_elements(yb2)), yb2,'1b',/overplot)

finalx2 = [x2[0:300], xb2[300:600],x2[301:3050],xb2[3250:3500],x2[3051:5550],xb2[6000:6200], x2[5501:*]]
finaly2 =[y2[0:300], yb2[300:600]-0.0005,y2[301:3050],yb2[3250:3500]-0.001,y2[3051:5550],yb2[6000:6200]-0.0001, y2[5501:*]]
p2 = plot(finalx2, finaly2,  '1.r', /overplot,sym_filled = 1. );,   sym_size = 0.1

;l = legend(target = [p1,p2,p4], position = [-0.6, 0.990], /data)
t1 = text(-0.6, 0.992, 'Ch1', color = 'black',/data)
t2 = text(-0.6, 0.991, 'Ch2', color = 'red',/data)
t4 = text(-0.6, 0.990, 'Ch4', color = 'blue',/data)
ts = text(-0.6, 0.988, 'S', color = 'green',/data)
ts = text(-0.59, 0.988, 'n', color = 'medium_purple',/data)
ts = text(-0.58, 0.988, 'a', color = 'hot_pink',/data)
ts = text(-0.57, 0.988, 'p', color = 'dodger_blue',/data)

;---------------------------------------------------------
;what about adding the snapshot data.
;--------------------------------

colorarr = ['deep_pink', 'magenta', 'medium_purple', 'hot_pink', 'light_pink', 'rosy_brown', 'chocolate', 'saddle_brown', 'maroon', 'orange_red', 'dark_orange', 'gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green', 'aquamarine', 'teal', 'steel_blue', 'dodger_blue', 'dark_blue', 'indigo','dark_slate_blue', 'blue_violet', 'purple','dim_grey', 'slate_grey', 'dark_slate_grey', 'khaki','black', 'light_cyan', 'lavender', 'peach_puff', 'pale_goldenrod','red']
restore,'/Users/jkrick/irac_warm/hd209458/hd209458_snap.sav'
;restore,'/Users/jkrick/irac_warm/hd209458/hd209458_snap_JI.sav'
;snapshots = snapshots2

 aorname = ['0045188864','0045189120','0045189376','0045189632','0045189888','0045190144','0045190400','0045190656','0045190912','0045191168','0045191424','0045191680','0045191936','0045192192','0045192704','0045195264','0045192960','0045193216','0045193472','0045193984','0045193728','0045195520','0045194240','0045194496','0045194752','0045195008','0045196288','0045195776','0045197312','0045196032','0045196544','0045196800','0045197056','0045197568','0045197824','0045198080','0045192448']
;--------------------------------
;start with the timing
 ;know the period
 period = 3.52474859  ;days
 period_err =0.00000038  ;days

 ;ok and I know when one transit is
 zero_phase = 2455942.56608; 2455928.46708D  ;UTC JD  roughly 01.16.12
 zero_phase = zero_phase -  2400000.5D  ;now in MJD

;transform bmjd from header into phase
;need to worry a bit about transforming JD zero phase to BMJD headers  
;this transformation is 8 min max.

;--------------------------------
;need a normalization factor for all snapshots together
;resistant_mean, snapshots.corrflux2,3,rmean, rsig, num

;print, 'normcorr', rmean, mean(snapshots.corrflux,/nan), median(snapshots.corrflux)
rmean = median(snapshots.corrflux)
;rmean2 = median(snapshots.corrflux2)
;print, 'normalization', rmean, rmean2
;--------------------------------
 bin_level = 360L

  for a = 0, n_elements(aorname) - 1 do begin
     print, 'working on a, color', a, aorname(a), colorarr(a)
   
;     print, 'nflux', n_elements(snapshots[a].corrflux)
;remove outliers
;get the nan's out of there.
     nonan = where(finite(snapshots[a].corrflux) gt 0)
;     print, 'nonan', n_elements(nonan)
                                ;see if it worked
     corrflux = snapshots[a].corrflux(nonan)
  ;   corrflux2 = snapshots[a].corrflux2(nonan)
     corrfluxerr = snapshots[a].corrfluxerr(nonan)
     bmjd = snapshots[a].bmjd(nonan)
                                ;print, 'nflux', n_elements(snapshots[a].corrflux)
     
;try getting rid of flux outliers.
;do some running mean with clipping
     start = 0
;     print, 'nflux', n_elements(corrflux)
     for ni = 100, n_elements(corrflux) -1, 100 do begin
        meanclip,corrflux[start:ni], m, s, subs = subs ;,/verbose
  ;      meanclip,corrflux2[start:ni], m2, s2, subs = subs2 ;,/verbose
                                ; print, 'good', subs+start
                                ;now keep a list of the good ones
        if ni eq 100 then good_ni = subs+start else good_ni = [good_ni, subs+start]
  ;      if ni eq 100 then good_ni2 = subs2+start else good_ni = [good_ni, subs2+start]
        start = ni + 1
     endfor
;     print, 'n good_ni', n_elements(good_ni)
     
                                ;see if it worked
     corrflux = corrflux[good_ni]
  ;   corrflux2 = corrflux2[good_ni]
     corrfluxerr = corrfluxerr[good_ni]
     bmjd = bmjd[good_ni]
;     print, 'nflux', n_elements(corrflux)

   ;must start with binning the corrected snapshot
   ; pick my own binning level 
    if n_elements(corrflux) lt bin_level then nframes = n_elements(corrflux) else nframes = bin_level
    nfits = long(n_elements(corrflux)) / nframes
    print, 'nframes, nfits', nframes, nfits
    ;help, nfits
   
;try binning the fluxes by subarray image
    bin_corrflux= dblarr(nfits)
 ;   bin_corrflux2= dblarr(nfits)
    bin_corrfluxerr= dblarr(nfits)
    bin_bmjd = dblarr(nfits)

    for si = 0L, long(nfits) - 1L do begin
       ;print, 'working on si', si, n_elements(corrflux), si*nframes, si*nframes + (nframes - 1)
       idata = corrflux[si*nframes:si*nframes + (nframes - 1)]
  ;     idata2 = corrflux2[si*nframes:si*nframes + (nframes - 1)]
       idataerr = corrfluxerr[si*nframes:si*nframes + (nframes - 1)]
       bin_corrflux[si] = mean(idata,/nan)
  ;     bin_corrflux2[si] = mean(idata2,/nan)
       bin_corrfluxerr[si] =   sqrt(total(idataerr^2))/ (n_elements(idataerr))
       bin_bmjd[si]= mean(bmjd[si*nframes:si*nframes + (nframes - 1)])
    endfor                      ;for each fits image
;        print, 'BIN_BMJDF, corr', n_elements(bin_bmjd), n_elements(bin_corrflux)

;--------------------------------
  ;now work on phasing the timearray,  and plotting
 
    bmjd = bin_bmjd; snapshots[sa].bmjd
    bmjd_dist = bmjd - zero_phase ; how many UTC away from the transit center
    phase =( bmjd_dist / period )- fix(bmjd_dist/period)

   ;ok, but now I want -0.5 to 0.5, not 0 to 1
   ;need to be careful here because subtracting half a phase will put things off, need something else
    ;print, ' phase',  phase, format = '(A,F0)'    
    pa = where(phase gt 0.5,pacount)
    ;print, 'pa', pa
    if pacount gt 0 then phase[pa] = phase[pa] - 1.0
;    print, ' phase',  phase, format = '(A,F0)'
;    print, 'nphase, corr', n_elements(phase), n_elements(bin_corrflux)

;    normcorr = bin_corrflux/mean(bin_corrflux,/nan)
    if a ne 20 then begin
       ps = errorplot(phase, bin_corrflux /rmean, bin_corrfluxerr / rmean, '1s',/overplot, name = string(aorname(a)),sym_filled = 1.,color = colorarr[a], errorbar_capsize = 0.005, xrange =[-0.7, 0.7], yrange = [0.985, 1.005], sym_size = 0.5)
       
; ps = plot(phase, bin_corrflux2 /rmean2,  '1s',/overplot, name = string(aorname(a)),color = colorarr[a], sym_size = 2);, sym_filled = 1.)
    endif
endfor

end
