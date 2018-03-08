pro add_mips70

restore, '/Users/jkrick/idlbin/objectnew.sav'

;read in akari data  
readcol, '/Users/jkrick/spitzer/mips/mips70/mosaic/apex/commandline/mosaic_extract.tbl',id, ra, dec, x, y, flux, ap5, snr, chi2dof,format="A"

;mips70  = { mips70ra:0D, mips70dec:0D, mips70xcenter:0D, mips70ycenter:0D,mips70flux:0D, mips70apflux:0D, mips70snr:0D, mips70chi2dof:0D,mips70match:0D, mips70matchdist:0D}

;add to catalog
newmips70 = {  mips70apflux:0D, mips70snr:0D, mips70chi2dof:0D}
ab = replicate(newmips70, n_elements(objectnew.ra))
objectnew = struct_addtags(objectnew, ab)

print, 'clearing matches'
;clear the matches that are already there.
objectnew.mips70xcenter=0
objectnew.mips70ycenter=0
objectnew.mips70ra=0
objectnew.mips70dec=0
objectnew.mips70flux=-99
objectnew.mips70fluxerr=-99
objectnew.mips70mag=99
objectnew.mips70magerr=99
objectnew.mips70bckgrnd=0
objectnew.mips70isoarea=0
objectnew.mips70fwhm=0
objectnew.mips70flags=0
objectnew.mips70nndist=0
objectnew.mips70match = -999
objectnew.mips70matchdist = 999


;do some matching
m = n_elements(ra)
mmatch = fltarr(m)
mmatch[*] = -999

print, 'starting matching'
for q=0,m -1 do begin
   dist=sphdist( ra(q), dec(q),objectnew.ra,objectnew.dec,/degrees)
   sep=min(dist,ind)

   separation = dist[reverse(sort(dist))]
   index = reverse(sort(dist))
   for test = long(0), n_elements(dist) - 1 do begin
      if objectnew[index[test]].mips24flux gt 0 and separation[index[test]] lt 0.00083 then begin;not acceptable, need to find a new better match
        ind = index[test]
     endif
   endfor

 print, 'found a match', objectnew[ind].ra,objectnew[ind].dec,ra(q), dec(q), objectnew[ind].mips24flux      

;   if (sep LT 0.00083)  then begin  ;3"

      if objectnew[ind].mips70match ge 0 then begin   ;if there was a previous match
         if sep lt objectnew[ind].mips70matchdist then begin   ;new match is a closer better match
           ;how do I put the old match into the unmatched list?  need to know its 'q'
            ;aka which number it is in mips24object
            mmatch[objectnew[ind].mips70match] = -999
            ;put the new match into the matched pile
            mmatch[q]=ind
            objectnew[ind].mips70ra = ra(q)
            objectnew[ind].mips70dec = dec(q)
            objectnew[ind].mips70flux = flux(q)
            objectnew[ind].mips70apflux = ap5(q)
            objectnew[ind].mips70xcenter = x(q)
            objectnew[ind].mips70ycenter = y(q)
            objectnew[ind].mips70snr = snr(q)
            objectnew[ind].mips70chi2dof = chi2dof(q)
            objectnew[ind].mips70match = q
            objectnew[ind].mips70matchdist = sep
         endif
      endif else begin          ;if there was no previous match
         mmatch[q]=ind
         objectnew[ind].mips70ra = ra(q)
         objectnew[ind].mips70dec = dec(q)
         objectnew[ind].mips70flux = flux(q)
         objectnew[ind].mips70apflux = ap5(q)
         objectnew[ind].mips70xcenter = x(q)
         objectnew[ind].mips70ycenter = y(q)
         objectnew[ind].mips70snr = snr(q)
         objectnew[ind].mips70chi2dof = chi2dof(q)
         objectnew[ind].mips70match = q
         objectnew[ind].mips70matchdist = sep
      endelse
;   endif 
endfor
matched=where(mmatch GE 0)
nonmatched = where(mmatch lt 0)

print, "matched", n_elements(where(mmatch GT 0))
print, "nonmatched", n_elements(where(mmatch lt 0))

;print out any that do not have a match


openw, outlunred, '/Users/jkrick/spitzer/mips/mips70/mosaic/apex/commandline/nonmatched.reg', /get_lun
printf, outlunred, 'fk5'
for num= 0,  n_elements(nonmatched) -1. do  printf, outlunred, 'circle( ',ra(nonmatched[num]), dec(nonmatched[num]) , ' 10")'
close, outlunred
free_lun, outlunred

openw, outlunred, '/Users/jkrick/spitzer/mips/mips70/mosaic/apex/commandline/matched.reg', /get_lun
printf, outlunred, 'fk5'

;for num= 0,  n_elements(matched) -1. do  begin
;   printf, outlunred, 'circle( ',ra(matched[num]), dec(matched[num]) , ' 9") # color=blue'
;endfor

;try printing out the object ra and dec of the matched sources, not
;that from the mips70 catalog

m = where(objectnew.mips70flux gt 0)
for num= 0,  n_elements(m) -1. do  begin
   printf, outlunred, 'circle( ',objectnew[m[num]].ra, objectnew[m[num]].dec , ' 9") # color=blue'
endfor


close, outlunred
free_lun, outlunred


;save, objectnew, filename='/Users/jkrick/idlbin/objectnew_akari.sav'

end
