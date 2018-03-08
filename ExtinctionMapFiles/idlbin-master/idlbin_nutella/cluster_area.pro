pro cluster_area

restore, '/Users/jkrick/idlbin/object.old.sav'
m=-0.12
b = 5.1
totfaintcount = 0
totbrightcount = 0
canddec = 69.045792;69.04481 
candra = 264.69266;264.68160 

openw, outlunbright, '/Users/jkrick/nep/clusters/bright.reg', /get_lun
printf, outlunbright, 'fk5'
openw, outlunfaint, '/Users/jkrick/nep/clusters/faint.reg', /get_lun
printf, outlunfaint, 'fk5'

;want a half annulus between 500 and 750 kpc, 0.017degr and 0.026degr
 good = where(sphdist(objectnew.ra, objectnew.dec, candra, canddec, /degrees) lt 0.026 and  $
              sphdist(objectnew.ra, objectnew.dec, candra, canddec, /degrees) gt 0.017  $
              and objectnew.ra ge candra $
              and objectnew.irac1mag gt 0 and objectnew.irac1mag lt 90 and objectnew.acsmag gt 0 $
              and objectnew.acsmag lt 90) 


print,"number in half annulus",  n_elements(good)
memberbrightarr = fltarr(n_elements(good))
memberfaintarr = fltarr(n_elements(good))

;want members in that annulus
for count = 0, n_elements(good) - 1 do begin
   limithigh = (m*objectnew[good[count]].irac1mag) + b +0.5    ; 2*sig1;0.5 
   limitlow = (m*objectnew[good[count]].irac1mag) + b - 0.5    ;2*sig1;0.5 
   printf, outlunfaint, 'circle( ', objectnew[good[count]].ra, objectnew[good[count]].dec, ' 3")'

   if (objectnew[good[count]].acsmag - objectnew[good[count]].irac1mag LT limithigh) and $
      (objectnew[good[count]].acsmag - objectnew[good[count]].irac1mag GT limitlow) and $
      (objectnew[good[count]].irac1mag lt 24.6)  and (objectnew[good[count]].irac1mag gt 22.6) then begin
      memberfaintarr(totfaintcount) = good[count]
      totfaintcount = totfaintcount + 1
   endif
   if (objectnew[good[count]].acsmag - objectnew[good[count]].irac1mag LT limithigh) and $
      (objectnew[good[count]].acsmag - objectnew[good[count]].irac1mag GT limitlow) and $
      (objectnew[good[count]].irac1mag le 22.6)  and (objectnew[good[count]].irac1mag gt 20.6) then begin
      memberbrightarr(totbrightcount) = good[count]
      totbrightcount = totbrightcount + 1
      printf, outlunbright, 'circle( ', objectnew[good[count]].ra, objectnew[good[count]].dec, ' 3")'

   endif

endfor

memberbrightarr = memberbrightarr[0:totbrightcount-1]   
print, 'number of bright members in half annnulus', totbrightcount
memberfaintarr = memberfaintarr[0:totfaintcount-1]   
print, 'number of faint membersin half annnulus', totfaintcount


;---------------------------------
;now count eveerything less than 0.017 degrees
totfaintcount = 0
totbrightcount = 0

good = where(sphdist(objectnew.ra, objectnew.dec, candra, canddec, /degrees) le 0.017 $
              and objectnew.irac1mag gt 0 and objectnew.irac1mag lt 90 and objectnew.acsmag gt 0 $
              and objectnew.acsmag lt 90) 

print,"number in half annulus",  n_elements(good)
memberbrightarr = fltarr(n_elements(good))
memberfaintarr = fltarr(n_elements(good))

;want members in that annulus
for count = 0, n_elements(good) - 1 do begin
   limithigh = (m*objectnew[good[count]].irac1mag) + b +0.5    ; 2*sig1;0.5 
   limitlow = (m*objectnew[good[count]].irac1mag) + b - 0.5    ;2*sig1;0.5 

   if (objectnew[good[count]].acsmag - objectnew[good[count]].irac1mag LT limithigh) and $
      (objectnew[good[count]].acsmag - objectnew[good[count]].irac1mag GT limitlow) and $
      (objectnew[good[count]].irac1mag lt 24.6)  and (objectnew[good[count]].irac1mag gt 22.6) then begin
      memberfaintarr(totfaintcount) = good[count]
      totfaintcount = totfaintcount + 1
   endif
   if (objectnew[good[count]].acsmag - objectnew[good[count]].irac1mag LT limithigh) and $
      (objectnew[good[count]].acsmag - objectnew[good[count]].irac1mag GT limitlow) and $
      (objectnew[good[count]].irac1mag le 22.6)  and (objectnew[good[count]].irac1mag gt 20.6) then begin
      memberbrightarr(totbrightcount) = good[count]
      totbrightcount = totbrightcount + 1
   endif

endfor

memberbrightarr = memberbrightarr[0:totbrightcount-1]   
print, 'number of bright members inside 500', totbrightcount
memberfaintarr = memberfaintarr[0:totfaintcount-1]   
print, 'number of faint members inside 500', totfaintcount

close, outlunbright
free_lun, outlunbright
close, outlunfaint
free_lun, outlunfaint


end

