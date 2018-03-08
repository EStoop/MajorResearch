pro eclipse_test
restore, '/Users/jkrick/irac_warm/hd209458/hd209458_ch4.sav'

  parinfo = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0.D]}, n_elements(pa0))
                               

  afargs = {FLUX:y, DX:adx, DY:ady, T:x, ERR:yerr}
  pa = mpfit('eclipse_func', pa0, FUNCTARGS=afargs, PERROR=spa, BESTNORM=achi, DOF=adof, COVAR = COV, status = status, errmsg = errmsg, parinfo = parinfo) ;,/quiet)
  print, 'status', status
  print, errmsg
  achi = achi / adof
  print, 'reduced chi squared', achi
  if keyword_set(sine) then begin  ;use a fitting function with a sine curve or not
     model_fit = pa[0] * (1+ pa[1] * adx + pa[2] * ady + pa[3] * adx * ady + pa[4] * adx * adx + pa[5] * ady * ady)+ pa[6]*sin(x/pa[7] + pa[8]) 
  endif else begin
     model_fit = pa[0] * (1+ pa[1] * adx + pa[2] * ady + pa[3] * adx * ady + pa[4] * adx * adx + pa[5] * ady * ady)
  endelse


end


function eclipse_func,x,p
; p is parameters, f0, d, t1, t3, dt
; t is time
; Start of ingress is t1, end of ingress is t2, start of egress is t3, end of egress is t4
; assume ingress and egress have same durations, dt 
   n = n_elements(x)
   model = fltarr(n) + p[0]

   t2 = p[2] + p[4]
   t4 = p[3] + p[4]

; Out of eclipse
;    ptr = where(x le p[2], count)
;    if (count gt 0) then model[ptr] = p[0]

;    ptr = where(x gt t4, count)
;    if (count gt 0) then model[ptr] = p[0]

; Beginning eclipse
   ptr = where(x gt p[2] and x le t2, count)
   if (count gt 0) then $
       model[ptr] = model[ptr] - p[1] * (1.0 - (t2 - x[ptr]) / p[4])

; In eclipse
   ptr = where(x gt t2 and x le p[3], count)
   if (count gt 0) then model[ptr] = model[ptr] - p[1]

; Leaving eclipse
   ptr = where(x gt p[3] and x le t4, count)
   if (count gt 0) then $
       model[ptr] = model[ptr] + p[1] * ((x[ptr] - p[3]) / p[4] - 1.0)

return, model
end
