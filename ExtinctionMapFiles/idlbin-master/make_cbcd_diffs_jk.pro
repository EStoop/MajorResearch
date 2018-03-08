pro make_cbcd_diffs_jk, aor, channel
;  print, N_params()
;  if (N_params() ne 2) then begin
;    print, 'Syntax == make_cbcd_diffs, reqkey, channel'
;    return
;  endif
  
  dir = '~/external/S19.2repro/'
  
  ;;;;Set the IMASK bits to keep
  bits_to_keep = [2,3,4,5,6,7,8,9,10,11,12,13,14]
  bitmask = FIX(TOTAL(2L^bits_to_keep)) 
  
  saor = strn(aor)
  sch = strn(channel)
  sha_loc = dir + 'S19.1SHA/r' + saor + '/ch' + sch + '/bcd/'
  repro_loc = dir + '00' + saor
;  print, 'Processing request key ' + saor, sha_loc, repro_loc

  bcdrepro_files = file_search(repro_loc + '/IRAC.'+ sch + '.00' + saor + '*.bcd_fp.fits', COUNT=nfiles)
  bcdsha_files = file_search(sha_loc + '/SPITZER_I' + sch + '_' + saor + '*_bcd.fits', COUNT=ccount)
  imask_files = file_search(sha_loc + '/SPITZER_I' + sch + '_' + saor + '*_bimsk.fits', COUNT=icount)
 
  imask_2_files = file_dirname(imask_files,/mark) + file_basename(imask_files,'bimsk.fits') + 'bimsk2.fits'
  diff_files = file_dirname(bcdsha_files,/mark) + file_basename(bcdsha_files,'bcd.fits') + 'diff.fits'
  
  if (ccount eq nfiles AND icount EQ nfiles) then begin
    for j = 0L, nfiles-1 do begin
;      print_pct, j, nfiles
      bcdrepro = DOUBLE(readfits(bcdrepro_files[j], bcdreproh, /SILENT))
      bcdsha = DOUBLE(readfits(bcdsha_files[j],bcdshah,/SILENT))
      imask = readfits(imask_files[j],imskh,/SILENT) AND bitmask
      sxaddhist,['MAKE_CBCD_DIFFS.pro Ran '+systime(),'Keeping bits '+STRCOMPRESS(STRJOIN(bits_to_keep,','),/remove)],imskh
      sxaddhist,['MAKE_CBCD_DIFFS.pro Ran '+systime(),'Difference between '+file_basename(bcdrepro_files[j]),' and '+file_basename(bcdsha_files[j])],bcdreproh
      sxaddpar,cbcdh,'BITPIX',-64
      print, 'writing', diff_files[j]
      writefits,imask_2_files[j],imask,imskh
      writefits,diff_files[j],bcdsha-bcdrepro,bcdreproh
    endfor
  endif
  
return
end
