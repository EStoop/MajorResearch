function get_nflux

restore, '/Users/jkrick/nutella/idlbin/objectnew.sav'
nflux = fltarr(n_elements(objectnew.ra))

for i = 0L, n_elements(objectnew.ra) - 1 do begin

if objectnew[i].umag gt 0 and objectnew[i].umag ne 99  then nflux(i) = nflux(i) + 1 ;and objectnew[i].umagerr lt 1.0 
if objectnew[i].gmag gt 0 and objectnew[i].gmag ne 99   then nflux(i) = nflux(i) + 1 ;and objectnew[i].umagerr lt 1.0 
if objectnew[i].rmag gt 0 and objectnew[i].rmag ne 99  then nflux(i) = nflux(i) + 1 ; and objectnew[i].umagerr lt 1.0 
if objectnew[i].imag gt 0 and objectnew[i].imag ne 99   then nflux(i) = nflux(i) + 1 ;and objectnew[i].umagerr lt 1.0
if objectnew[i].acsmag gt 0 and objectnew[i].acsmag ne 99 then nflux(i)=nflux(i)+1
if objectnew[i].zmagbest gt 0 and objectnew[i].zmagbest ne 99 then nflux(i)=nflux(i)+1
if objectnew[i].flamjmag gt 0 and objectnew[i].flamjmag ne 99 then nflux(i) = nflux(i) + 1
if objectnew[i].wircjmag gt 0 and objectnew[i].wircjmag ne 99 then nflux(i) = nflux(i) + 1
if objectnew[i].wirchmag gt 0 and objectnew[i].wirchmag ne 99 then nflux(i) = nflux(i) + 1
if objectnew[i].wirckmag gt 0 and objectnew[i].wirckmag ne 99 then nflux(i) = nflux(i) + 1
if objectnew[i].tmassjmag gt 0 and objectnew[i].tmassjmag ne 99 then nflux(i) = nflux(i) + 1
if objectnew[i].tmasshmag gt 0 and objectnew[i].tmasshmag ne 99 then nflux(i) = nflux(i) + 1
if objectnew[i].tmasskmag gt 0 and objectnew[i].tmasskmag ne 99 then nflux(i) = nflux(i) + 1
if objectnew[i].irac1mag gt 0 and objectnew[i].irac1mag ne 99 then nflux(i) = nflux(i) + 1
if objectnew[i].irac2mag gt 0 and objectnew[i].irac2mag ne 99 then nflux(i) = nflux(i) + 1
if objectnew[i].irac3mag gt 0 and objectnew[i].irac3mag ne 99 and objectnew[i].irac3flux gt 2.* objectnew[i].irac3fluxerr then nflux(i) = nflux(i) + 1
if objectnew[i].irac4mag gt 0 and objectnew[i].irac4mag ne 99 and objectnew[i].irac4flux gt 2.* objectnew[i].irac4fluxerr  then nflux(i) = nflux(i) + 1
if objectnew[i].mips24mag gt 0 and objectnew[i].mips24mag ne 99 then nflux(i) = nflux(i) + 1
endfor

return, nflux
end

