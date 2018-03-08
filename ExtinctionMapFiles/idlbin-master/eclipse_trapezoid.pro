function eclipse_trapezoid, p, FLUXfit=fluxfit, Xfit=xfit, yERRfit=yerrfit
; p is parameters, f0, d, t1, t3, dt
; t is time
; Start of ingress is t1, end of ingress is t2, start of egress is t3, end of egress is t4
; assume ingress and egress have same durations, dt 
   n = n_elements(xfit)
   model = fltarr(n) + p[0]

   t2 = p[2] + p[4]
   t4 = p[3] + p[4]

; Out of eclipse
    ptr = where(xfit le p[2], count)
    if (count gt 0) then model[ptr] = p[0]

    ptr = where(xfit gt t4, count)
    if (count gt 0) then model[ptr] = p[0]

; Beginning eclipse
   ptr = where(xfit gt p[2] and xfit le t2, count)
   if (count gt 0) then $
       model[ptr] = model[ptr] - p[1] * (1.0 - (t2 - xfit[ptr]) / p[4])

; In eclipse
   ptr = where(xfit gt t2 and xfit le p[3], count)
   if (count gt 0) then model[ptr] = model[ptr] - p[1]

; Leaving eclipse
   ptr = where(xfit gt p[3] and xfit le t4, count)
   if (count gt 0) then $
       model[ptr] = model[ptr] + p[1] * ((xfit[ptr] - p[3]) / p[4] - 1.0)

 model = (fluxfit - model) / yerrfit

return, model
end
