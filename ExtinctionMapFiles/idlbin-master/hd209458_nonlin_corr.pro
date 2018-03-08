pro hd209458_nonlin_corr


restore, '/Users/jkrick/irac_warm/hd209458/hd209458_plot_ch2.sav'
;normalize the flux
y = bin_sub+add
y2 = y ;/ mean(y)
x2 = x
;p2 = plot(x2, y2,  '1.r', sym_filled = 1.,xtitle = 'Orbital Phase', ytitle = ' Normalized Flux', title = 'HD209458', name = 'Ch2',yrange = [0.985, 1.01], xrange = [-0.7,0.7]);,   sym_size = 0.1

;l = legend(target = [p1,p2,p4], position = [-0.6, 0.990], /data)
;t2 = text(-0.6, 0.991, 'Ch2', color = 'red',/data)
;ts = text(-0.6, 0.988, 'Snap', color = 'green',/data)

;---------------------------------------------------------
;what about adding the snapshot data.
;--------------------------------
 colorarr = ['deep_pink', 'magenta', 'medium_purple', 'hot_pink', 'light_pink', 'rosy_brown', 'chocolate', 'saddle_brown', 'maroon', 'orange_red', 'dark_orange', 'gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green', 'aquamarine', 'teal', 'steel_blue', 'dodger_blue', 'dark_blue', 'indigo','dark_slate_blue', 'blue_violet', 'purple','dim_grey', 'slate_grey', 'dark_slate_grey', 'khaki','black', 'light_cyan', 'lavender', 'peach_puff', 'pale_goldenrod','red']

restore,'/Users/jkrick/irac_warm/hd209458/hd209458_snap.sav'
aorname = ['0045188864','0045189120','0045189376','0045189632','0045189888','0045190144','0045190400','0045190656','0045190912','0045191168','0045191424','0045191680','0045191936','0045192192','0045192704','0045195264','0045192960','0045193216','0045193472','0045193984','0045193728','0045195520','0045194240','0045194496','0045194752','0045195008','0045196288','0045195776','0045197312','0045196032','0045196544','0045196800','0045197056','0045197568','0045197824','0045198080','0045192448']
;--------------------------------;start with the timing
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
; resistant_mean, snapshots.corrflux,3,rmean, rsig, num

;print, 'normcorr', rmean, mean(snapshots.corrflux,/nan), median(snapshots.corrflux)
 rmean = median(snapshots.corrflux)

;--------------------------------
 bin_level = 300L
 bin_corrfluxarr = fltarr(n_elements(aorname))
 starefluxarr = fltarr(n_elements(aorname))
 bin_pmapcorrarr = fltarr(n_elements(aorname))
 for a = 0,  n_elements(aorname) - 1 do begin
    print, 'working on a, color', a, aorname(a), colorarr(a)
   
;     print, 'nflux', n_elements(snapshots[a].corrflux)
;remove outliers
;get the nan's out of there.
     nonan = where(finite(snapshots[a].corrflux) gt 0)
;     print, 'nonan', n_elements(nonan)
                                ;see if it worked
     flux = snapshots[a].flux(nonan)
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
     flux = flux[good_ni]
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
   
;-------------------------------------------------------
;try binning the fluxes all the way down to one point
    bin_flux = dblarr(nfits)
    bin_corrflux= dblarr(nfits)
    bin_corrfluxerr= dblarr(nfits)
    bin_bmjd = dblarr(nfits)

    for si = 0L, long(nfits) - 1L do begin
       ;print, 'working on si', si, n_elements(corrflux), si*nframes, si*nframes + (nframes - 1)
       idata = corrflux[si*nframes:si*nframes + (nframes - 1)]
       idataerr = corrfluxerr[si*nframes:si*nframes + (nframes - 1)]
       iflux = flux[si*nframes:si*nframes + (nframes - 1)]
       bin_corrflux[si] = mean(idata,/nan)
       bin_flux[si] = mean(iflux,/nan)
       bin_corrfluxerr[si] =   sqrt(total(idataerr^2))/ (n_elements(idataerr))
       bin_bmjd[si]= mean(bmjd[si*nframes:si*nframes + (nframes - 1)])
    endfor                      ;for each fits image
;        print, 'BIN_BMJDF, corr', n_elements(bin_bmjd), n_elements(bin_corrflux)

    ; the pmap correction is the division of corrflux / flux
    bin_pmapcorrarr[a] = bin_flux / bin_corrflux
    bin_corrfluxarr[a] = bin_corrflux
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


    ;-------------------------------
    ;what is the staring mode flux at each of those binned phases?
    pg = fltarr(n_elements(x2))
    pgcount = 0
    for g = 0, n_elements(x2) -1 do begin
       ;print, g, x2(g), phase+ 0.01, phase - 0.01, format = '(I, F0,F0,F0)'
      ; if x2(g) lt (phase + 0.01)  then print, 'found one', g, x2(g), phase
       if x2(g) gt (phase - 0.005)  and x2(g) lt (phase + 0.005) then begin
          ;print, 'found one', g, x2(g), phase
          pg(pgcount) = g
          pgcount = pgcount + 1
       endif

    endfor

       pg = pg[0:pgcount - 1]
       print, 'pgcount', pgcount
    ;I don't know why this isn't working

    ;pg = where(x2 gt (phase - 0.005)  and x2 lt (phase + 0.005), pgcount)
    ;print, 'pg', pgcount, pg
    starefluxarr[a] = mean(y2(pg)) ; only want one value for now



;    ps = errorplot(phase, bin_corrflux /rmean, bin_corrfluxerr / rmean, '1s',/overplot, name = string(aorname(a)),sym_filled = 1.,color = colorarr[a], errorbar_capsize = 0.005)

endfor

;check
;print, 'stare', starefluxarr
;print, 'bin_corrflux', bin_corrfluxarr

;----------------------------------------'
;plotting
;make each point a different color based on all the other stuff so I can weed out the eclipses/transits by eye.
ps = plot(bin_pmapcorrarr, bin_corrfluxarr / starefluxarr, '1*', xtitle = 'pmap correction', ytitle = 'snap / stare', title = 'HD209458',nodata = 1)

for a = 0, n_elements(aorname) -1 do begin
   ps1 = text(bin_pmapcorrarr[a], bin_corrfluxarr[a] / starefluxarr[a], '$\diamondsuit$', alignment = 0.5, vertical_alignment = 0.5,color = colorarr[a], /data, font_size = 25)
   if a eq 4 then tps = text(bin_pmapcorrarr[a], bin_corrfluxarr[a] / starefluxarr[a], '$\Chi$', alignment = 0.5, vertical_alignment = 0.5, /data)
   if a eq 18 then tps = text(bin_pmapcorrarr[a], bin_corrfluxarr[a] / starefluxarr[a], '$\Chi$', alignment = 0.5, vertical_alignment = 0.5, /data)
   if a eq 13 then tps = text(bin_pmapcorrarr[a], bin_corrfluxarr[a] / starefluxarr[a], '$\Chi$', alignment = 0.5, vertical_alignment = 0.5, /data)
   if a eq 21 then tps = text(bin_pmapcorrarr[a], bin_corrfluxarr[a] / starefluxarr[a], '$\Chi$', alignment = 0.5, vertical_alignment = 0.5, /data)
   if a eq 22 then tps = text(bin_pmapcorrarr[a], bin_corrfluxarr[a] / starefluxarr[a], '$\Chi$', alignment = 0.5, vertical_alignment = 0.5, /data)

endfor


end
