pro irac_completeness
ps_open, filename='/Users/jkrick/spitzer/irac/completeness.ps',/portrait,/color

redcolor = FSC_COLOR("Red", !D.Table_Size-2)

!p.multi=[0,2,2]
!X.thick=3
!Y.thick=3

restore, '/Users/jkrick/idlbin/objectnew.sav'
name = 'irac1mag'


for ch = 0, 3 do begin
   if ch eq 0 then begin
      a = where(objectnew.irac1flux gt 0 and objectnew.irac1flux lt 90)
      plothist, objectnew[a].irac1flux, xhist, yhist, bin = 0.01, xrange=[1E-2,100], /xlog, $
                thick=3, charthick=3, xtitle='micro Jy', ytitle='ch 1 number'
;only want to fit y from x =0.5 - 1.1
      front = where(xhist gt 0.5 and xhist lt 0.51)
      back = where(xhist gt 1.1 and xhist lt 1.11)
   endif
  if ch eq 1 then begin
      a = where(objectnew.irac2flux gt 0 and objectnew.irac2flux lt 90)
      plothist, objectnew[a].irac2flux, xhist, yhist, bin = 0.01, xrange=[1E-2,100], /xlog, $
                thick=3, charthick=3, xtitle='micro Jy', ytitle='ch 2 number'
      front = where(xhist gt 0.4 and xhist lt 0.41)
      back = where(xhist gt 1.1 and xhist lt 1.11)
   endif
  if ch eq 2 then begin
      a = where(objectnew.irac3flux gt 0 and objectnew.irac3flux lt 90)
      plothist, objectnew[a].irac3flux, xhist, yhist, bin = 0.01, xrange=[1E-2,100], /xlog, $
                thick=3, charthick=3, xtitle='micro Jy', ytitle='ch 3 number'
      front = where(xhist gt 0.5 and xhist lt 0.51)
      back = where(xhist gt 1.1 and xhist lt 1.11)
   endif
  if ch eq 3 then begin
      a = where(objectnew.irac4flux gt 0 and objectnew.irac4flux lt 90)
      plothist, objectnew[a].irac4flux, xhist, yhist, bin = 0.01, xrange=[1E-2,100], /xlog, $
                thick=3, charthick=3, xtitle='micro Jy', ytitle='ch 4 number'
      front = where(xhist gt 0.8 and xhist lt 0.81)
      back = where(xhist gt 1.2 and xhist lt 1.21)
   endif



;new xrange and yrange
   limitxhist = xhist[front:back]
   limityhist = yhist[front:back]
   start = [500.0,-10000.0]

;what do I think the errors are? 
;this is just number counts, make all equal errors.
   noise = fltarr(n_elements(limityhist))
   noise = noise + 1.0

   result= MPFITFUN('linear',limitxhist,limityhist, noise, start, /quiet)  
   yfit = result(0)*xhist + result(1)
   oplot, xhist, yfit, thick = 3, color = redcolor, linestyle = 0

;where does the real histogram equal 95% of the fitted line?
   dist = (yfit)  - (yhist)
   
   ratio = where(yhist gt  0.95*(yfit) )
   print, "ch", ch + 1, '  95% completeness = ', xhist(ratio(0))
endfor

ps_close, /noprint,/noid

end
