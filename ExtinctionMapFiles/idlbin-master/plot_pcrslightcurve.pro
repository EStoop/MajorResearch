pro plot_pcrslightcurve


;ps_start, filename= '/Users/jkrick/nutella/spitzer/planets_proposal/pixnoise/centroids.ps'
;ps_open,filename= '/Users/jkrick/nutella/spitzer/planets_proposal/pixnoise/centroids.ps',/color,/landscape
  ifps, '/Users/jkrick/irac_warm/pcrsstare/centroids.ps',19.5,2.5,/color,/post,/encap

;

  restore, filename='/Users/jkrick/irac_warm/pcrsstare/pixnoise.sav'
  
;  a = where (calstar.aor eq 39161344 and calstar.ch eq 2 and strmid(calstar.filename, 42, 2) ne '12' and strmid(calstar.filename, 42, 2) ne '00', count)
;print,calstar[a].bcdycen
;plot, calstar[a].bcdxcen[1:*] -15.0, calstar[a].bcdycen[1:*] - 15.0, xrange = [-0.8,0.8], yrange = [-0.6,0.8], psym = 1, xstyle = 1, ystyle = 1
  

reqkey = [ 0042023424, 0042023680, 0042023936, 0042024192, 0042024448, 0042024704, 0042024960, 0042025216, 0042025472, 0042025728, 0042025984, 0042026240, 0042026496, 0042026752, 0042027008, 0042027264, 0042027520, 0042027776, 0042028032, 0042028288]
  

  xrange=[-.75,0.75]
  yrange=[-.6,.75]
  x0 = 15.0
  y0 = 15.0
;  allreq = reqkey_set
;  iuniq_req = UNIQ(allreq,BSORT(allreq))
;  sort_bytime = BSORT(time[iuniq_req])
;  reqkey = allreq[iuniq_req[sort_bytime]]    ;;; Sorted by time
  nreq = n_elements(reqkey)
  cred = FSC_COLOR('Indian Red',!d.table_size-3)
  cgrey = FSC_COLOR('Medium Gray', !d.table_size-2)
  multiplot,[nreq,1],/initialize
;  dx = calstar[a].bcdxcen-x0
;  dy = calstar[a].bcdxcen-y0
  xdiff = xrange[1]-xrange[0]
  ydiff = yrange[1]-yrange[0]
  FOR i = 0,nreq-1 do begin
     ;iuse = where(reqkey_set EQ reqkey[i],nuse)
       a = where (calstar.aor eq reqkey(i) and calstar.ch eq 2 and strmid(calstar.filename, 33, 2) ne '12' and strmid(calstar.filename, 33, 2) ne '00', count)
      ; a = where (calstar.aor eq reqkey(i) and calstar.ch eq 2 , count)

     multiplot
    ; if nuse NE 0 THEN BEGIN
        xtn = !x.tickname
        ytn = !y.tickname
        !x.tickname = ' ' 
        !y.tickname = ' '
        plot,xrange-0.5,yrange-0.5,xstyle=1,ystyle=1,/nodata,xticklen=1.0,yticklen=1.0,ytickinterval=1.0,xtickinterval=1.0,$
             color=cgrey,xtickname=' ',ytickname=' '
        !x.tickname = xtn
        !y.tickname = ytn
        IF i EQ 0 THEN !y.title='!4D!3Y (pix)' ELSE !y.title=' '
        plot,calstar[a].bcdxcen - x0, calstar[a].bcdycen - y0,ystyle=1,psym=1,xcharsize=0.4,ycharsize=0.5,xstyle=1,xrange=xrange,yrange=yrange,   /noerase
        avg_x = MEAN(calstar[a].bcdxcen - x0)
        avg_y = MEAN(calstar[a].bcdycen - y0)
        xyouts,avg_x,avg_y-0.2,'('+STRING(avg_x,format='(f5.2)')+','+STRING(avg_y,format='(f5.2)')+')',alignment=0.5,charsize=0.6
        !y.title= ' '
        plots,0,0,psym=symcat(7),color=cred

        ;testing how this plotting works
        plot, calstar.sclktime, calstar.bcdflux, psym = 3
    ; ENDIF
  ENDFOR
xyouts,0.5,0.2,'!4D!3X from intended position (pix)',alignment=0.5,/Norm,charsize=0.8
multiplot, /reset

;ps_close, /noprint,/noid

;ps_end
ENDPS, '/Users/jkrick/irac_warm/pcrsstare/centroids.ps',/pdf

end
