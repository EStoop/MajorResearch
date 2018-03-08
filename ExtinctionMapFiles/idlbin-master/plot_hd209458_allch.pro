pro plot_hd209458_allch

restore, '/Users/jkrick/irac_warm/hd209458/hd209458_plot_ch1.sav'
;normalize the flux
y = bin_sub+add
y1 = y / mean(y)
x1 = x
 p1 = plot(x1, y1,  '1.',  color = 'black',xtitle = 'Orbital Phase', ytitle = ' Normalized Flux', title = 'HD209458 snaps with pmap correction', sym_filled = 1.,yrange = [0.985, 1.01], xrange = [-0.7,0.7], name = 'Ch1') ;, sym_size = 0.1, xrange = [2,4]

restore, '/Users/jkrick/irac_warm/hd209458/hd209458_plot_ch4.sav'
;normalize the flux
y = bin_sub+add
y4 = y / mean(y)
x4=x
p4 = plot(x4, y4,  '1.b', sym_filled = 1.,/overplot, name = 'Ch4') ;,  sym_size = 0.1

restore, '/Users/jkrick/irac_warm/hd209458/hd209458_plot_ch2.sav'
;normalize the flux
y = bin_sub+add
y2 = y / mean(y)
x2 = x
p2 = plot(x2, y2,  '1.r', sym_filled = 1.,/overplot, name = 'Ch2');,   sym_size = 0.1

;l = legend(target = [p1,p2,p4], position = [-0.6, 0.990], /data)
t1 = text(-0.6, 0.992, 'Ch1', color = 'black',/data)
t2 = text(-0.6, 0.991, 'Ch2', color = 'red',/data)
t4 = text(-0.6, 0.990, 'Ch4', color = 'blue',/data)
;ts = text(-0.6, 0.988, 'Snap', color = 'green',/data)

;---------------------------------------------------------
;what about adding the snapshot data.
;--------------------------------
 colorarr = [  'SILVER'   , 'SKY_BLUE'  , 'SLATE_BLUE' , 'SLATE_GREY' ,  'PLUM' ,  'SPRING_GREEN'  ,  'STEEL_BLUE'   ,  'green', 'TEAL' ,  'blue' ,  'TOMATO'  ,  'TURQUOISE' ,  'VIOLET'  ,  'YELLOW_GREEN',' PEACH_PUFF','   PERU',   ' POWDER_BLUE', 'dark_blue','dark_cyan','dark_goldenrod','dark_gray','dark_green','dark_khaki','dark_magenta','dark_olive_green','dark_orange','dark_orchid','dark_red','dark_salmon','dark_sea_green','dark_slate_blue','dark_slate_gray','dark_turquoise','dark_violet','deep_pink','deep_sky_blue']

;restore,'/Users/jkrick/irac_warm/hd209458/hd209458_snap.sav'
restore,'/Users/jkrick/irac_warm/hd209458/hd209458_snap_JI.sav'
snapshots = snapshots2

 aorname = ['0045188864','0045189120','0045189376','0045189632','0045189888','0045190144','0045190400','0045190656','0045190912','0045191168','0045191424','0045191680','0045191936','0045192192','0045192704','0045195264','0045192960','0045193216','0045193472','0045193984','0045193728','0045195520','0045194240','0045194496','0045194752','0045195008','0045196288','0045195776','0045197312','0045196032','0045196544','0045196800','0045197056','0045197568','0045197824','0045198080']
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
print, 'normalization', rmean
;--------------------------------
 bin_level = 300L

  for a = 0, n_elements(aorname) - 1 do begin
     print, 'working on a, color', a, aorname(a), colorarr(a)
   
;     print, 'nflux', n_elements(snapshots[a].corrflux)
;remove outliers
;get the nan's out of there.
     nonan = where(finite(snapshots[a].corrflux) gt 0)
;     print, 'nonan', n_elements(nonan)
                                ;see if it worked
     corrflux = snapshots[a].corrflux(nonan)
     corrfluxerr = snapshots[a].corrfluxerr(nonan)
     bmjd = snapshots[a].bmjd(nonan)
                                ;print, 'nflux', n_elements(snapshots[a].corrflux)
     
;try getting rid of flux outliers.
;do some running mean with clipping
     start = 0
;     print, 'nflux', n_elements(corrflux)
     for ni = 100, n_elements(corrflux) -1, 100 do begin
        meanclip,corrflux[start:ni], m, s, subs = subs ;,/verbose
                                ; print, 'good', subs+start
                                ;now keep a list of the good ones
        if ni eq 100 then good_ni = subs+start else good_ni = [good_ni, subs+start]
        start = ni + 1
     endfor
;     print, 'n good_ni', n_elements(good_ni)
     
                                ;see if it worked
     corrflux = corrflux[good_ni]
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
    bin_corrfluxerr= dblarr(nfits)
    bin_bmjd = dblarr(nfits)

    for si = 0L, long(nfits) - 1L do begin
       ;print, 'working on si', si, n_elements(corrflux), si*nframes, si*nframes + (nframes - 1)
       idata = corrflux[si*nframes:si*nframes + (nframes - 1)]
       idataerr = corrfluxerr[si*nframes:si*nframes + (nframes - 1)]
       bin_corrflux[si] = mean(idata,/nan)
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
    ps = errorplot(phase, bin_corrflux /rmean, bin_corrfluxerr / rmean, '1s',/overplot, name = string(aorname(a)),sym_filled = 1.,color = colorarr[a], errorbar_capsize = 0.005)

endfor

end
