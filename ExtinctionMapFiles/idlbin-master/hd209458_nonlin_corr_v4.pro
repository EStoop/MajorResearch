pro hd209458_nonlin_corr


restore, '/Users/jkrick/irac_warm/hd209458/hd209458_plot_ch2.sav'
;normalize the flux
y = bin_sub+add
y2 = y ;/ mean(y)
x2 = x

;---------------------------------------------------------
; the snapshot data.
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
 bin_corrfluxarr = fltarr(n_elements(aorname))
 bin_fluxarr = fltarr(n_elements(aorname))
 starefluxarr = fltarr(n_elements(aorname))
 bin_pmapcorrarr = fltarr(n_elements(aorname))

 for a = 0,  n_elements(aorname) - 1 do begin
    print, 'working on a, color', a, aorname(a), colorarr(a)
   
;remove outliers
;get the nan's out of there.
     nonan = where(finite(snapshots[a].corrflux) gt 0)
     flux = snapshots[a].flux(nonan)
     corrflux = snapshots[a].corrflux(nonan)
     corrfluxerr = snapshots[a].corrfluxerr(nonan)
     bmjd = snapshots[a].bmjd(nonan)
     
;try getting rid of flux outliers.
;do some running mean with clipping
     start = 0
     print, 'nflux', n_elements(corrflux)
     for ni = 100, n_elements(corrflux) -1, 100 do begin
        meanclip,corrflux[start:ni], m, s, subs = subs ;,/verbose
                                ; print, 'good', subs+start
                                ;now keep a list of the good ones
        if ni eq 100 then good_ni = subs+start else good_ni = [good_ni, subs+start]
        start = ni + 1
     endfor
     
     flux = flux[good_ni]
     corrflux = corrflux[good_ni]
     corrfluxerr = corrfluxerr[good_ni]
     bmjd = bmjd[good_ni]
     print, 'nflux', n_elements(corrflux)
   
;-------------------------------------------------------
;try binning the fluxes all the way down to one point
     bin_corrfluxarr[a] = mean(corrflux,/nan)
     bin_fluxarr[a] = mean(flux,/nan)
     bin_bmjd[a] = mean(bmjd,/nan)
    ; the pmap correction is the division of corrflux / flux
     bin_pmapcorrarr[a] = bin_fluxarr[a] / bin_corrfluxarr[a]

     print, 'n binnedcorrflux', n_elements(bin_corrfluxarr[a])
;--------------------------------
  ;now work on phasing the timearray,  and plotting
 
    bmjd = bin_bmjd[a]; snapshots[sa].bmjd
    bmjd_dist = bmjd - zero_phase ; how many UTC away from the transit center
    phase =( bmjd_dist / period )- fix(bmjd_dist/period)
    print, 'n phase', n_elements(phase)
   ;ok, but now I want -0.5 to 0.5, not 0 to 1
   ;need to be careful here because subtracting half a phase will put things off, need something else
    ;print, ' phase',  phase, format = '(A,F0)'    
 ;   pa = where(phase gt 0.5,pacount)
    ;print, 'pa', pa
 ;   if pacount gt 0 then phase[pa] = phase[pa] - 1.0

    if phase gt 0.5 then phase = phase - 1.0
    print, ' phase',  phase, format = '(A,F0)'
    print, 'nphase, corr', n_elements(phase), n_elements(bin_corrfluxarr[a])


    ;-------------------------------
    ;what is the staring mode flux at each of those binned phases?
    aw = where(x2 gt phase -0.001 and x2 lt phase + 0.001, awcount)
    print, 'x2(aw)', n_elements(aw), x2(aw)

    if awcount gt 1 then begin
       starefluxarr[a] = mean(y2(aw)) ; only want one value for now
    endif else begin ; there is no staring mode comparison becuase the phase is on an eclipse or transit
       starefluxarr[a] = 10000.
    endelse


;    ps = errorplot(phase, bin_corrflux /rmean, bin_corrfluxerr / rmean, '1s',/overplot, name = string(aorname(a)),sym_filled = 1.,color = colorarr[a], errorbar_capsize = 0.005)

endfor

;check
print, 'stare', starefluxarr
print, 'bin_corrflux', bin_corrfluxarr

;----------------------------------------'
;plotting
;make each point a different color based on all the other stuff so I can weed out the eclipses/transits by eye.
ps = plot(bin_pmapcorrarr, bin_corrfluxarr / starefluxarr, '1*', xtitle = 'pmap correction', ytitle = 'snap / stare', title = 'HD209458', yrange = [1.03, 1.05],nodata = 1)

for a = 0, n_elements(aorname) -1 do begin
   ps1 = text(bin_pmapcorrarr[a], bin_corrfluxarr[a] / starefluxarr[a], '$\diamondsuit$', alignment = 0.5, vertical_alignment = 0.5,color = colorarr[a], /data, font_size = 25)
   if a eq 4 then tps = text(bin_pmapcorrarr[a], bin_corrfluxarr[a] / starefluxarr[a], '$\Chi$', alignment = 0.5, vertical_alignment = 0.5, /data)
   if a eq 18 then tps = text(bin_pmapcorrarr[a], bin_corrfluxarr[a] / starefluxarr[a], '$\Chi$', alignment = 0.5, vertical_alignment = 0.5, /data)
   if a eq 13 then tps = text(bin_pmapcorrarr[a], bin_corrfluxarr[a] / starefluxarr[a], '$\Chi$', alignment = 0.5, vertical_alignment = 0.5, /data)
   if a eq 21 then tps = text(bin_pmapcorrarr[a], bin_corrfluxarr[a] / starefluxarr[a], '$\Chi$', alignment = 0.5, vertical_alignment = 0.5, /data)
   if a eq 22 then tps = text(bin_pmapcorrarr[a], bin_corrfluxarr[a] / starefluxarr[a], '$\Chi$', alignment = 0.5, vertical_alignment = 0.5, /data)

endfor


end

   ; pg = fltarr(n_elements(x2))
   ; pgcount = 0
    ;help, x2
    ;help, phase
   ; for g = 0, n_elements(x2) -1 do begin
       ;print, g, x2(g), phase+ 0.01, phase - 0.01, format = '(I, F0,F0,F0)'
                                ;            if x2(g) lt (phase + 0.01)  then print, 'found one', g, x2(g), phase
   ;    aw = where(x2(g) gt (phase - 0.002)  and x2(g) lt (phase + 0.002) ,awcount)

   ;    print, g, awcount
      ; if x2(g) gt (phase - 0.005)  and x2(g) lt (phase + 0.005) then begin
          ;print, 'found one', g, x2(g), phase
      ;    pg(pgcount) = g
      ;    pgcount = pgcount + 1
      ; endif
 ;  endfor

      ; pg = pg[0:pgcount - 1]
       ;print, 'pgcount', pgcount
    ;I don't know why this isn't working

    ;pg = where(x2 gt (phase - 0.005)  and x2 lt (phase + 0.005), pgcount)
    ;print, 'pg', pgcount, pg
    ;starefluxarr[a] = mean(y2(pg)) ; only want one value for now
 
