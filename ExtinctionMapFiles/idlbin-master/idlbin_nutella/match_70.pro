pro match_70
restore, '/Users/jkrick/idlbin/object.sav'


p = 0
mips70object = replicate({mips70ob, mips70xcenter:0D, mips70ycenter:0D, mips70ra:0D, mips70dec:0D, mips70flux:0D, mips70fluxerr:0D, mips70mag:0D, mips70magerr:0D, mips70bckgrnd:0D, mips70isoarea:0D,mips70fwhm:0D, mips70flags:0D},500)

openr, luni, "/Users/jkrick/spitzer/mips/mips70/mosaic/Combine/mosaic.asc", /get_lun, error=err
if (err ne 0.) then print, "file didn't open", ﻿!ERROR_STATE.MSG  

WHILE (NOT EOF(luni)) DO BEGIN

   READF, luni, num,ra, dec, x, y, fluxiso, fluxerriso,fluxauto, fluxerrauto,kron,back,isoarea,a, b, theta, flag, fwhm,elon, ellip
; I do not know what the magnitudes or the correct fluxes are.  there is an area issue.

   mips70object[p] ={mips70ob, x,y, ra, dec, fluxauto*23.5045, fluxerrauto*23.5045, fluxauto, fluxerrauto,back, isoarea, fwhm, flag}

   p = p + 1
   
 endwhile
mips70object = mips70object[0:p-1]
print, "there are ",p," mips70 objects"



;---------------------------------------------------------

;match mips70
; create initial arrays
m=n_elements(mips70object.mips70ra)
ir=n_elements(object.ra)

irmatch=fltarr(ir)
mmatch=fltarr(m)
irmatch[*]=-999
mmatch[*]=-999


print,'Matching mips70 to object'
print,"Starting at "+systime()
dist=irmatch
dist[*]=0

for q=0,m-1 do begin

   dist=sphdist( mips70object[q].mips70ra, mips70object[q].mips70dec,object.ra,object.dec,/degrees)
   sep=min(dist,ind)
   if (sep LT 0.0008)  then begin
      mmatch[q]=ind
      object[ind].mips70ra = mips70object[q].mips70ra
      object[ind].mips70dec = mips70object[q].mips70dec
      object[ind].mips70xcenter = mips70object[q].mips70xcenter
      object[ind].mips70ycenter = mips70object[q].mips70ycenter
      object[ind].mips70mag = mips70object[q].mips70mag
      object[ind].mips70magerr = mips70object[q].mips70magerr
      object[ind].mips70fwhm = mips70object[q].mips70fwhm
      object[ind].mips70isoarea = mips70object[q].mips70isoarea
      object[ind].mips70flags = mips70object[q].mips70flags
      object[ind].mips70bckgrnd = mips70object[q].mips70bckgrnd
      object[ind].mips70flux = mips70object[q].mips70flux
      object[ind].mips70fluxerr = mips70object[q].mips70fluxerr
   endif 
endfor
print, "matched", n_elements(where(mmatch GT 0))
print, "nonmatched", n_elements(where(mmatch lt 0))

print,"Finished at "+systime()


plothist, object.mips70mag, xhist, yhist, /noprint, yrange=[0,10]

save, object, filename='/Users/jkrick/idlbin/object.sav'

end
