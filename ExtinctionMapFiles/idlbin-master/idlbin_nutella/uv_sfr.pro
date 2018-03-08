function uv_sfr

common sharevariables

;-----------------------------------------------------------
;correct u and g mags for Galactic extinction
;using Schlegel average of 0.042 = E(B-V) in the area of the clusters.
;using Rv=3.1 to convert to Au and Ag (schlegel app. B)

Au = 0.042*5.155
Ag = 0.042*3.793

objectnew[member].umag = objectnew[member].umag - Au
objectnew[member].gmag = objectnew[member].gmag - Ag

;-----------------------------------------------------------
;convert u and g mags into fluxes
fu = 10^((objectnew[member].umag + 48.6)/(-2.5))
fg = 10^((objectnew[member].gmag + 48.6)/(-2.5))

;print, 'fu', fu

fuerr_low = 10^(((objectnew[member].umag - objectnew[member].umagerr) + 48.6)/(-2.5))
fuerr_up= 10^(((objectnew[member].umag + objectnew[member].umagerr) + 48.6)/(-2.5))
fuerr_up = abs(fu - fuerr_up)
fuerr_low = abs(fuerr_low - fu)
fuerr = (fuerr_up + fuerr_low) / 2

fgerr_low = 10^(((objectnew[member].gmag - objectnew[member].gmagerr) + 48.6)/(-2.5))
fgerr_up = 10^(((objectnew[member].gmag + objectnew[member].gmagerr) + 48.6)/(-2.5))
fgerr_up = abs(fg - fgerr_up)
fgerr_low = abs(fgerr_low - fg)
fgerr = (fgerr_up + fgerr_low) / 2


a = where(fgerr lt 0.5*fg and fuerr lt 0.5*fu)
print,'number with errors lt 0.3 flux ', n_elements(a)


;-----------------------------------------------------------
;derive Beta from u - g
;fit power law slope to u and g
;f(lambda) = c*lambda^(beta)
lambda = [3400,5200]
start = [1.0E-35,2.0]
beta = fltarr(n_elements(member))
const = fltarr(n_elements(member))
for i = 0, n_elements(fu) - 1 do begin
   flux = [fu(i), fg(i)]
   fluxerr = [fuerr(i), fgerr(i)]
   result= MPFITFUN('powerlaw',lambda,flux, fluxerr, start, bestnorm = bn, perror = perr,/quiet) ;ICL
   if objectnew[member[i]].umag lt 90 and objectnew[member[i]].umag gt 0 and objectnew[member[i]].gmag lt 90 $
      and objectnew[member[i]].gmag gt 0 then begin
      beta[i] = result(1)
      const[i] = result(0)
   endif else begin
      beta[i] = 99
      const[i] = 0
   endelse

;  print, 'test beta', beta(i), objectnew[member[i]].umag, objectnew[member[i]].gmag, fu(i), fg(i), const[i], const(i) * 3400^(beta[i]),  const(i) * 5200^(beta[i])

endfor

plothist, beta
;-----------------------------------------------------------
;convert beta into an extinction at 1600 
A1600 = 4.43 + 1.99*beta  ;in mags
fu_corr = fltarr(n_elements(member))
;correct fu for A1600
for i = 0, n_elements(member) - 1 do begin
   if objectnew[member[i]].umag lt 90 and objectnew[member[i]].umag gt 0 and objectnew[member[i]].gmag lt 90 $
      and objectnew[member[i]].gmag gt 0 and beta[i] lt 4. then begin
       fu_corr[i] = fu[i]*10^(0.4*A1600[i])
       
   endif else begin
      fu_corr[i] = 0
   endelse
endfor


;-----------------------------------------------------------
;convert fu into Lu
z = 1.
Lu = lonarr(n_elements(member))

Lu = (((fu_corr*4*!Pi*luminositydistance^2)/(1 + z)) * ((3.0856D24)^2) ) ; include conversion from cm to Mpc 
Lu_solar = Lu / 3.839D33 ;and erg/s to solar lum.
;print, "Lu_solar", Lu_solar

;convert Lu into sfr
sfr_uv = Lu * 1.4E-28   ;kennicutt

plothist, sfr_uv
print, sfr_uv


;for i = 0, n_elements(member) - 1 do print, format='(I10, F10.2, F10.2, E10.2, E10.2, E10.2,E10.2,F10.2, F10.2, E10.2, E10.2, F10.2 )', member(i), objectnew[member(i)].umag, objectnew[member(i)].gmag, fu(i), fg(i), fuerr(i), fgerr(i), beta(i), pcerror(i), lu(i), lu_solar(i), sfr_uv(i)
;-----------------------------------------------------------

return, 0

end

