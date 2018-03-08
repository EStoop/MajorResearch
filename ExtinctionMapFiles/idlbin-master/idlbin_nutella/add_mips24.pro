pro addmips24

numoobjects = 10000

mips24object = replicate({mips24ob, mips24xcenter:0D, mips24ycenter:0D, mips24ra:0D, mips24dec:0D, mips24flux:0D, mips24fluxerr:0D, mips24mag:0D, mips24magerr:0D, mips24bckgrnd:0D, iter:0D,sharp:0D, mips24match:0D, mips24nndist:0D, mips24matchdist:0D},numoobjects)


openr, lunw, "/Users/jkrick/Spitzer/mips/mips24/dao/newastrometry/mips24.ra.phot", /get_lun
w = 0
WHILE (NOT EOF(lunw)) DO BEGIN
  READF, lunw, num, ra, dec, xcen,ycen, mag, magerr, back, niter, sh, junk, junk
   mips24object[w] = {mips24ob, xcen, ycen, ra,dec, 1.432*(10^((25-mag)/2.5)),  10^((25-(mag -magerr))/2.5) -(10^((25-mag)/2.5)) , mag,magerr,back, niter,sh,0.,0.,0.}
   w = w +1
endwhile
mips24object = mips24object[0:w - 1]
print, "there are ",w," mips24 objects"
close, lunw
free_lun, lunw

;determine nndist
for nn = 0l, n_elements(mips24object.mips24ra) - 1 do begin
   a = where(sqrt( (mips24object[nn].mips24ra - mips24object.mips24ra)^2 + (mips24object[nn].mips24dec - mips24object.mips24dec)^2 ) le 0.0008 )
   mips24object[nn].mips24nndist = n_elements(a)
endfor

restore, '/Users/jkrick/idlbin/objectnew.sav'

;clear the matches that are already there.

objectnew.mips24xcenter = 0 
objectnew.mips24ycenter = 0 
objectnew.mips24ra = 0 
objectnew.mips24dec = 0 
objectnew.mips24flux = -1 
objectnew.mips24fluxerr = -1 
objectnew.mips24mag = 99.
objectnew.mips24magerr = 99.
objectnew.mips24bckgrnd = 0 
objectnew.iter = 0
objectnew.sharp = 0 
objectnew.mips24nndist = 0
objectnew.mips24match = -999
objectnew.mips24matchdist = 999

;do some matching
m = n_elements(mips24object.mips24ra)
mmatch = fltarr(m)
mmatch[*] = -999
for q=0,m - 1 do begin
   dist=sphdist( mips24object[q].mips24ra, mips24object[q].mips24dec,objectnew.ra,objectnew.dec,/degrees)
   sep=min(dist,ind)
;   print, 'q', mips24object[q].mips24ra, mips24object[q].mips24dec, sep
   if (sep LT 0.0008)  then begin  ;4"
      if objectnew[ind].mips24match ge 0 then begin   ;if there was a previous match
         if sep lt objectnew[ind].mips24matchdist then begin   ;new match is a closer better match
           ;how do I put the old match into the unmatched list?  need to know its 'q'
            ;aka which number it is in mips24object
            mmatch[objectnew[ind].mips24match] = -999
            ;put the new match into the matched pile
            mmatch[q]=ind
            objectnew[ind].mips24ra = mips24object[q].mips24ra
            objectnew[ind].mips24dec = mips24object[q].mips24dec
            objectnew[ind].mips24xcenter = mips24object[q].mips24xcenter 
            objectnew[ind].mips24ycenter = mips24object[q].mips24ycenter 
            objectnew[ind].mips24flux = mips24object[q].mips24flux 
            objectnew[ind].mips24fluxerr = mips24object[q].mips24fluxerr 
            objectnew[ind].mips24mag = mips24object[q].mips24mag 
            objectnew[ind].mips24magerr = mips24object[q].mips24magerr
            objectnew[ind].mips24bckgrnd = mips24object[q].mips24bckgrnd 
            objectnew[ind].iter = mips24object[q].iter
            objectnew[ind].sharp = mips24object[q].sharp 
            objectnew[ind].mips24nndist = mips24object[q].mips24nndist
            objectnew[ind].mips24match = q
            objectnew[ind].mips24matchdist = sep

         endif
      endif else begin          ;if there was no previous match
         mmatch[q]=ind
         objectnew[ind].mips24ra = mips24object[q].mips24ra
         objectnew[ind].mips24dec = mips24object[q].mips24dec
         objectnew[ind].mips24xcenter = mips24object[q].mips24xcenter 
         objectnew[ind].mips24ycenter = mips24object[q].mips24ycenter 
         objectnew[ind].mips24flux = mips24object[q].mips24flux 
         objectnew[ind].mips24fluxerr = mips24object[q].mips24fluxerr 
         objectnew[ind].mips24mag = mips24object[q].mips24mag 
         objectnew[ind].mips24magerr = mips24object[q].mips24magerr
         objectnew[ind].mips24bckgrnd = mips24object[q].mips24bckgrnd 
         objectnew[ind].iter = mips24object[q].iter
         objectnew[ind].sharp = mips24object[q].sharp 
         objectnew[ind].mips24nndist = mips24object[q].mips24nndist
         objectnew[ind].mips24match = q
         objectnew[ind].mips24matchdist = sep

      endelse
   endif 
endfor
matched=where(mmatch GE 0)
nonmatched = where(mmatch lt 0)

print, "matched", n_elements(where(mmatch GT 0))
print, "nonmatched", n_elements(where(mmatch lt 0))

;print out any that do not have a match
openw, outlunred, '/Users/jkrick/spitzer/mips/mips24/dao/newastrometry/nonmatch.reg', /get_lun
printf, outlunred, 'fk5'

for num= 0,  n_elements(nonmatched) -1. do begin
   
   printf, outlunred, 'circle( ', mips24object[nonmatched[num]].mips24ra,mips24object[nonmatched[num]].mips24dec, ' 3")'
endfor
close, outlunred
free_lun, outlunred

;these are nonmatches becuase they are noise or blends.



save, objectnew, filename='/Users/jkrick/idlbin/objectnew.test.sav'

end
