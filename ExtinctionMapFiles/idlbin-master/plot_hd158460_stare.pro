pro plot_hd158460_stare
  
restore,  '/Users/jkrick/irac_warm/snapshots/HD158460_stare_v2.sav'

  aorname = ['r42506496','r42051584']
  
  for a = 0, n_elements(aorname) - 1 do begin
 
     if a eq 0 then begin
        xarr = xarr_0
        yarr = yarr_0
        timearr = timearr_0
        st = plot(timearr, yarr,'6r1.', yrange = [14.5,15.5],xtitle = 'Time(hours)',ytitle = 'Y pix',  title = 'HD158460_stare') ;
        t = text(7, 14.85, aorname(a), color = 'red', /data, /current)
       st2 = plot(timearr, xarr,'6r1.',yrange = [14.5,15.5], xtitle = 'Time(hours)',ytitle = 'X pix',  title = 'HD158460_stare') ;
                                ;st3 = errorplot(timearr, fluxarr,fluxerrarr,'6r1.', xtitle = 'Time(hours)',ytitle = 'Flux',  title = 'HD158460_stare') ;
                                ; st3 = plot(timearr, fluxarr,'6r1.', yrange = [4.0, 4.3],xtitle = 'Time(hours)',ytitle = 'Flux',  title = 'HD158460_stare') ;
                                ; st4 = plot(timearr, corrfluxarr,'6r1.', yrange = [4.0, 4.15],xtitle = 'Time(hours)',ytitle = 'Corrected Flux',  title = 'HD158460_stare') ;
        t = text(7, 15.1, aorname(a), color = 'red', /data, /current)
     endif
     
     if a gt 0 then begin
        xarr = xarr_1
        yarr = yarr_1
        timearr = timearr_1
        st.Select
        st = plot( timearr, yarr, '6b1.',/overplot,/current)
        t = text(7, 15.2, aorname(a), color = 'blue', /data, /current)
       st2.Select
        st2 = plot( timearr, xarr- 0.2, '6b1.',/overplot, /current)
        t = text(7, 14.7, aorname(a), color = 'blue', /data, /current)

      ;  st3.Select
                                ;st3 = errorplot(timearr, fluxarr-0.2,fluxerrarr,'6r1.',/overplot, /current)
     ;   st3 = plot(timearr, fluxarr-0.02,'6b1.',/overplot, /current)
    ;    st4.Select
                                ;st4 = errorplot(timearr, fluxarr-0.2,fluxerrarr,'6r1.',/overplot, /current)
     ;   st4 = plot(timearr, corrfluxarr-0.02,'6b1.',/overplot, /current)
        
     endif
     
   
  endfor                        ;for each AOR
  
   
; save these files
  st.Save , '/Users/jkrick/irac_warm/snapshots/yvstime_stare.png'
  st2.Save ,  '/Users/jkrick/irac_warm/snapshots/xvstime_stare.png'
 ; st3.Save, '/Users/jkrick/irac_warm/snapshots/fluxvstime_stare.png'
 ; st4.Save, '/Users/jkrick/irac_warm/snapshots/corrfluxvstime_stare.png'
 ; xy.Save ,  '/Users/jkrick/irac_warm/snapshots/xy_stare.png'


end
