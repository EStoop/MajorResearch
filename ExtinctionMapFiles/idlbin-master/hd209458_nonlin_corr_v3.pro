pro hd209458_nonlin_corr_v3


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
;--------------------------------

;--------------------------------
;need a normalization factor for all snapshots together
; resistant_mean, snapshots.corrflux,3,rmean, rsig, num

;print, 'normcorr', rmean, mean(snapshots.corrflux,/nan), median(snapshots.corrflux)
 rmean = median(snapshots.corrflux)

;corrflux = dblarr(n_elements(snapshots.corrflux)
;flux = dblarr(n_elements(snapshots.corrflux)
;bmjd = dblarr(n_elements(snapshots.corrflux)
for a = 0, n_elements(aorname) -1 do begin
   if a eq 0 then begin
      corrflux = snapshots[a].corrflux 
      corrfluxerr = snapshots[a].corrfluxerr
      flux = snapshots[a].flux 
      bmjd = snapshots[a].bmjd
   endif else begin
      corrflux = [corrflux, snapshots[a].corrflux]
      corrfluxerr = [corrflux, snapshots[a].corrfluxerr]
      flux = [flux, snapshots[a].flux]
      bmjd = [bmjd, snapshots[a].bmjd]
   endelse
endfor

;--------------------------------
 ;remove outliers
;get the nan's out of there.
 nonan = where(finite(corrflux) gt 0)
;     print, 'nonan', n_elements(nonan)
                                ;see if it worked
 flux = flux(nonan)
 corrflux = corrflux(nonan)
 corrfluxerr = corrfluxerr(nonan)
 bmjd = bmjd(nonan)
                                ;print, 'nflux', n_elements(snapshots[a].corrflux)
 
;try getting rid of flux outliers.
;do some running mean with clipping
; start = double(0)
; help, start
;;     print, 'nflux', n_elements(corrflux)
; for ni = double(100), n_elements(corrflux) -1, double(100) do begin
;    if ni eq 100 then help, ni
;;print, 'ni, start', ni, start, n_elements(corrflux)
;   meanclip,corrflux[start:ni], m, s, subs = subs ;,/verbose
;;now keep a list of the good ones
;    if ni eq 100 then good_ni = subs+start else good_ni = [good_ni, subs+start]
;    start = ni + 1
; endfor
 
;see if it worked
; flux = flux[good_ni]
; corrflux = corrflux[good_ni]
; corrfluxerr = corrfluxerr[good_ni]
; bmjd = bmjd[good_ni]
 
; the pmap correction is the division of corrflux / flux
 pmapcorrarr = flux / corrflux
 
;--------------------------------
;now work on phasing the timearray,  and plotting
 ;start with the timing
 ;know the period
 period = 3.52474859            ;days
 period_err =0.00000038         ;days
 
;ok and I know when one transit is
 zero_phase = 2455942.56608     ; 2455928.46708D  ;UTC JD  roughly 01.16.12
 zero_phase = zero_phase -  2400000.5D ;now in MJD
 
;transform bmjd from header into phase
;need to worry a bit about transforming JD zero phase to BMJD headers  
;this transformation is 8 min max.
 
 bmjd_dist = bmjd - zero_phase  ; how many UTC away from the transit center
 phase =( bmjd_dist / period )- fix(bmjd_dist/period)
 
 ;ok, but now I want -0.5 to 0.5, not 0 to 1
 ;need to be careful here because subtracting half a phase will put things off, need something else
 ;print, ' phase',  phase, format = '(A,F0)'    
 pa = where(phase gt 0.5,pacount)
 if pacount gt 0 then phase[pa] = phase[pa] - 1.0

  ;-------------------------------
 ;what is the staring mode flux at each of those  phases?

 starefluxarr =dblarr(n_elements(phase)) 
 print, 'starting for loop'
 for s = 0, n_elements(phase) -1 do begin
;    print, 'working on', s, n_elements(phase)
    pg = fltarr(n_elements(x2))
    pgcount = 0

    aw = where(x2 gt  (phase(s) - 0.002)  and x2 lt (phase(s) + 0.002), awcount)
    starefluxarr(s) = mean(y2(aw)) 

;    for g = 0, n_elements(x2) -1 do begin
                                ;print, g, x2(g), phase+ 0.01, phase - 0.01, format = '(I, F0,F0,F0)'
                                ; if x2(g) lt (phase + 0.01)  then print, 'found one', g, x2(g), phase
;       if x2(g) gt (phase(s) - 0.002)  and x2(g) lt (phase(s) + 0.002) then begin
                                ;print, 'found one', g, x2(g), phase
;          pg(pgcount) = g
;          pgcount = pgcount + 1
;       endif
       
;    endfor
    
   ; pg = pg[0:pgcount - 1]
                                ;print, 'pgcount', pgcount
                                ;I don't know why this isn't working
    
                                ;pg = where(x2 gt (phase - 0.005)  and x2 lt (phase + 0.005), pgcount)
                                ;print, 'pg', pgcount, pg
    
 endfor
 
 
;    ps = errorplot(phase, bin_corrflux /rmean, bin_corrfluxerr / rmean, '1s',/overplot, name = string(aorname(a)),sym_filled = 1.,color = colorarr[a], errorbar_capsize = 0.005)
 

;check
print, 'nelements check', n_elements(starefluxarr), n_elements(corrflux), n_elements(pmapcorrarr)
save, /all, filename= '/Users/jkrick/irac_warm/hd209458/hd209458_nonlincorr.sav'

;----------------------------------------'
;plotting

colorfull = dblarr(n_elements(corrflux))
for a = 0, n_elements(colorfull) - 1 do begin
   if a / 121797 le 1  and a/121797 gt 0 then colorfull[a] = colorarr[0]
   if a / 121797 le 2  and a/121797 gt 1 then colorfull[a] = colorarr[1]
   if a / 121797 le 3  and a/121797 gt 2 then colorfull[a] = colorarr[2]
   if a / 121797 le 4  and a/121797 gt 3 then colorfull[a] = colorarr[3]
   if a / 121797 le 5  and a/121797 gt 4 then colorfull[a] = colorarr[4]
   if a / 121797 le 6  and a/121797 gt 5 then colorfull[a] = colorarr[5]
   if a / 121797 le 7  and a/121797 gt 6 then colorfull[a] = colorarr[6]
   if a / 121797 le 8  and a/121797 gt 7 then colorfull[a] = colorarr[7]
   if a / 121797 le 9  and a/121797 gt 8 then colorfull[a] = colorarr[8]
   if a / 121797 le 10  and a/121797 gt 9 then colorfull[a] = colorarr[9]
   if a / 121797 le 11  and a/121797 gt 10 then colorfull[a] = colorarr[10]
   if a / 121797 le 12  and a/121797 gt 11 then colorfull[a] = colorarr[11]
   if a / 121797 le 13  and a/121797 gt 12 then colorfull[a] = colorarr[12]
   if a / 121797 le 14  and a/121797 gt 13 then colorfull[a] = colorarr[13]
   if a / 121797 le 15  and a/121797 gt 14 then colorfull[a] = colorarr[14]
   if a / 121797 le 16  and a/121797 gt 15 then colorfull[a] = colorarr[15]
   if a / 121797 le 17  and a/121797 gt 16 then colorfull[a] = colorarr[16]
   if a / 121797 le 18  and a/121797 gt 17 then colorfull[a] = colorarr[17]
   if a / 121797 le 19  and a/121797 gt 18 then colorfull[a] = colorarr[18]
   if a / 121797 le 20  and a/121797 gt 19 then colorfull[a] = colorarr[19]
   if a / 121797 le 21  and a/121797 gt 20 then colorfull[a] = colorarr[20]
   if a / 121797 le 22  and a/121797 gt 21 then colorfull[a] = colorarr[21]
   if a / 121797 le 23  and a/121797 gt 22 then colorfull[a] = colorarr[22]
   if a / 121797 le 24  and a/121797 gt 23 then colorfull[a] = colorarr[23]
   if a / 121797 le 25  and a/121797 gt 24 then colorfull[a] = colorarr[24]
   if a / 121797 le 26  and a/121797 gt 25 then colorfull[a] = colorarr[25]
   if a / 121797 le 27  and a/121797 gt 26 then colorfull[a] = colorarr[26]
   if a / 121797 le 28  and a/121797 gt 27 then colorfull[a] = colorarr[27]
   if a / 121797 le 29  and a/121797 gt 28 then colorfull[a] = colorarr[28]
   if a / 121797 le 30  and a/121797 gt 29 then colorfull[a] = colorarr[29]
   if a / 121797 le 31  and a/121797 gt 30 then colorfull[a] = colorarr[30]
   if a / 121797 le 32  and a/121797 gt 31 then colorfull[a] = colorarr[31]
   if a / 121797 le 33  and a/121797 gt 32 then colorfull[a] = colorarr[32]
   if a / 121797 le 34  and a/121797 gt 33 then colorfull[a] = colorarr[33]
   if a / 121797 le 35  and a/121797 gt 34 then colorfull[a] = colorarr[34]
   if a / 121797 le 36  and a/121797 gt 35 then colorfull[a] = colorarr[35]
   if a / 121797 le 37  and a/121797 gt 36 then colorfull[a] = colorarr[36]
 endfor

;print, colorfull[0],colorfull[122797],colorfull[353900],colorfull[100000],colorfull[200000],colorfull[300000]
;make each point a different color based on all the other stuff so I can weed out the eclipses/transits by eye.
ps = plot(pmapcorrarr, corrflux / starefluxarr, '1*', xtitle = 'pmap correction', ytitle = 'snap / stare', title = 'HD209458', yrange = [1.01, 1.06], sym_size=0.8);,nodata = 1);nodata = 1

;for a = 0, n_elements(pmapcorrarr) -1 do begin
;   ps1 = text(pmapcorrarr[a], corrflux[a] / starefluxarr[a], '$\diamondsuit$', alignment = 0.5, vertical_alignment = 0.5,color = colorfull[a], /data, font_size = 2)

;endfor


end
