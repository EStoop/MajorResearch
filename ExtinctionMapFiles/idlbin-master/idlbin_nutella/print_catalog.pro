pro print_catalog
  ;;for grad student Helen Davidge
  ;; want to print out ra, dec, optical,IRAC, MIPS, Akari photometry


  restore, '/Users/jkrick/idlbin/idlbin_nutella/objectnew_akari.sav'
  openw, outlun, '/Users/jkrick/external/nep/akari/davidge_catalog.txt', /get_lun

  printf, outlun, "num, zphot, ra,  dec, umag,  umagerr, gmag, gmagerr, rmag, rmagerr, acsmag, acsmagerr, zmagbest, zmagerrbest,  wircjmag,  wircjmagerr, wirchmag, wirchmagerr, wirckmag, wirckmagerr, irac1flux,irac1fluxerr,irac2flux,irac2fluxerr, irac3flux, irac3fluxerr, irac4flux, irac4fluxerr, Akari flux11, Akari err11, Akari flux15,  Akari err15, Akari flux18, Akari err18,mips24flux, mips24fluxerr, mips70flux, mips70fluxerr"

  for ind = 0, n_elements(objectnew) - 1 do begin
     printf, outlun, format='(I10,F10.2,F10.5,F10.5,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2, F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)',$
             ind, objectnew[ind].zphot, objectnew[ind].ra,  objectnew[ind].dec, objectnew[ind].umag,  objectnew[ind].umagerr,$
             objectnew[ind].gmag, objectnew[ind].gmagerr,$
             objectnew[ind].rmag, objectnew[ind].rmagerr, objectnew[ind].acsmag, objectnew[ind].acsmagerr,$
             objectnew[ind].zmagbest, objectnew[ind].zmagerrbest,  objectnew[ind].wircjmag,  objectnew[ind].wircjmagerr,  $
             objectnew[ind].wirchmag, $
             objectnew[ind].wirchmagerr, objectnew[ind].wirckmag, objectnew[ind].wirckmagerr,$
             objectnew[ind].irac1flux,objectnew[ind].irac1fluxerr,objectnew[ind].irac2flux,objectnew[ind].irac2fluxerr,$
             objectnew[ind].irac3flux, objectnew[ind].irac3fluxerr, objectnew[ind].irac4flux, objectnew[ind].irac4fluxerr,$
             objectnew[ind].flux11, objectnew[ind].err11, objectnew[ind].flux15, $
             objectnew[ind].err15, objectnew[ind].flux18, objectnew[ind].err18,$
             objectnew[ind].mips24flux, objectnew[ind].mips24fluxerr,$
             objectnew[ind].mips70flux, objectnew[ind].mips70fluxerr
  endfor

  
  free_lun, outlun
end


