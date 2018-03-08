pro add_specz

restore, '/Users/jkrick/nutella/idlbin/objectnew_swirc.sav'
;read in the excell spreadsheet of spec z's and photz's from the
;catalog.  

readcol, '/Users/jkrick/nutella/palomar/cosmic/01june2008/spectra/specz.txt', masknum, catnum, photz, xccz, xcczerr,	zxc,	xcr,	xcgrade,	emcz,	emczerr,	zem,	emgrade,	whichmethod, format="(I1,L5,F10.3, F10.2,F10.2, F10.3, F10.2, I1, F10.2,F10.2,F10.3,I1,A)"

speczobject = replicate({speczob, specz:0D, speczerr:0D, speczmethod:"A"},n_elements(objectnew.ra))

objectnew = struct_addtags(objectnew, speczobject)

for i = 0, n_elements(catnum) - 1 do begin
   if whichmethod(i) eq 'em' and emgrade(i) gt 1 then begin
      objectnew[catnum(i)].specz =zem(i)
      objectnew[catnum(i)].speczerr = emczerr(i)
      objectnew[catnum(i)].speczmethod = 'em'
   endif
   if whichmethod(i) eq 'xc' and xcgrade(i) gt 1 then begin
      objectnew[catnum(i)].specz =zxc(i)
      objectnew[catnum(i)].speczerr = xcczerr(i)
      objectnew[catnum(i)].speczmethod = 'xc'
   endif

   
endfor

save, objectnew, filename='/Users/jkrick/nutella/idlbin/objectnew.sav'

end

