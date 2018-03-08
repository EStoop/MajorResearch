pro addakari

restore, '/Users/jkrick/idlbin/objectnew.sav'

;read in akari data  
readcol, '/Users/jkrick/akari/cat_2.txt', hh, mm, ss, dd, dm, ds, flux11a, err11a, flux15a, err15a, flux18a, err18a,format="A"

;convert h,m,s into degrees
ra = 15.* tenv(hh, mm, ss)
dec = tenv(dd, dm, ds)

;places where flux is 0, really means undetected which is -99 in our
;catalog

a = where(flux11a eq 0.0)
flux11a(a) = -99.0
err11a(a) = -99.0

b = where(flux15a eq 0.0)
flux15a(b) = -99.0
err15a(b) = -99.0

c = where(flux18a eq 0.0)
flux18a(c) = -99.0
err18a(c) = -99.0

;add to catalog
akari  = { akarira:0D, akaridec:0D, flux11:0D, err11:0D, flux15:0D, err15:0D,flux18:0D, err18:0D, akarimatch:0D, akarimatchdist:0D}
ab = replicate(akari, n_elements(objectnew.ra))
objectnew = struct_addtags(objectnew, ab)

objectnew.akarimatch = -999
objectnew.akarimatchdist = 999

;do some matching
m = n_elements(ra)
mmatch = fltarr(m)
mmatch[*] = -999

print, 'starting matching'
for q=0,m -1 do begin
   dist=sphdist( ra(q), dec(q),objectnew.ra,objectnew.dec,/degrees)
   sep=min(dist,ind)
   if (sep LT 0.00055)  then begin  ;3"
      if objectnew[ind].akarimatch ge 0 then begin   ;if there was a previous match
         if sep lt objectnew[ind].akarimatchdist then begin   ;new match is a closer better match
           ;how do I put the old match into the unmatched list?  need to know its 'q'
            ;aka which number it is in mips24object
            mmatch[objectnew[ind].akarimatch] = -999
            ;put the new match into the matched pile
            mmatch[q]=ind
            objectnew[ind].akarira = ra(q)
            objectnew[ind].akaridec = dec(q)
            objectnew[ind].flux11 = flux11a(q)
            objectnew[ind].err11 = err11a(q)
            objectnew[ind].flux15 = flux15a(q)
            objectnew[ind].err15 = err15a(q)
            objectnew[ind].flux18 = flux18a(q)
            objectnew[ind].err18 = err18a(q)
            objectnew[ind].akarimatch = q
            objectnew[ind].akarimatchdist = sep
         endif
      endif else begin          ;if there was no previous match
         mmatch[q]=ind
         objectnew[ind].akarira = ra(q)
         objectnew[ind].akaridec = dec(q)
         objectnew[ind].flux11 = flux11a(q)
         objectnew[ind].err11 = err11a(q)
         objectnew[ind].flux15 = flux15a(q)
         objectnew[ind].err15 = err15a(q)
         objectnew[ind].flux18 = flux18a(q)
         objectnew[ind].err18 = err18a(q)
         objectnew[ind].akarimatch = q
         objectnew[ind].akarimatchdist = sep
       endelse
   endif 
endfor
matched=where(mmatch GE 0)
nonmatched = where(mmatch lt 0)

print, "matched", n_elements(where(mmatch GT 0))
print, "nonmatched", n_elements(where(mmatch lt 0))

;print out any that do not have a match


openw, outlunred, '/Users/jkrick/akari/nonmatched.reg', /get_lun
printf, outlunred, 'fk5'
for num= 0,  n_elements(nonmatched) -1. do  printf, outlunred, 'circle( ',ra(nonmatched[num]), dec(nonmatched[num]) , ' 3")'
close, outlunred
free_lun, outlunred

openw, outlunred, '/Users/jkrick/akari/matched.reg', /get_lun
printf, outlunred, 'fk5'

for num= 0,  n_elements(matched) -1. do  begin
   printf, outlunred, 'circle( ',ra(matched[num]), dec(matched[num]) , ' 3") # color=red'
endfor

close, outlunred
free_lun, outlunred


save, objectnew, filename='/Users/jkrick/idlbin/objectnew_akari.sav'

end
