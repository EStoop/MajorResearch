pro make_cbcd_diffs, aor, channel
  if (N_params() ne 2) then begin
    print, 'Syntax == make_cbcd_diffs, reqkey, channel'
    return
  endif
  
  dir = '/Users/jkrick/iracdata/FINAL_CHECK/IRAC012200/'
  ndir = dir + 'S18.18v2/'
  ;odir = dir + 'S18.7/'
  
  ;;;;Set the IMASK bits to keep
  bits_to_keep = [2,3,4,5,6,7,8,9,10,11,12,13,14]
  bitmask = FIX(TOTAL(2L^bits_to_keep)) 
  
  saor = strn(aor)
  sch = strn(channel)
  suffix = 'r' + saor + '/ch' + sch + '/bcd/'
  sstr = ndir + suffix
  print, 'Processing request key ' + saor

  bcd_files = file_search(sstr + '/SPITZER_I' + sch + '_' + saor + '_*_bcd.fits', COUNT=nfiles)
  cbcd_files = file_search(sstr + '/SPITZER_I' + sch + '_' + saor + '_*_cbcd.fits', COUNT=ccount)
  imask_files = file_search(sstr + '/SPITZER_I' + sch + '_' + saor + '_*_bimsk.fits', COUNT=icount)
  
  imask_2_files = file_dirname(imask_files,/mark) + file_basename(imask_files,'bimsk.fits') + 'bimsk2.fits'
  diff_files = file_dirname(cbcd_files,/mark) + file_basename(cbcd_files,'cbcd.fits') + 'cdiff.fits'
  
  if (ccount eq nfiles AND icount EQ nfiles) then begin
    for j = 0L, nfiles-1 do begin
      print_pct, j, nfiles
      bcd = DOUBLE(readfits(bcd_files[j], bcdh, /SILENT))
      cbcd = DOUBLE(readfits(cbcd_files[j],cbcdh,/SILENT))
      imask = readfits(imask_files[j],imskh,/SILENT) AND bitmask
      sxaddhist,['MAKE_CBCD_DIFFS.pro Ran '+systime(),'Keeping bits '+STRCOMPRESS(STRJOIN(bits_to_keep,','),/remove)],imskh
      sxaddhist,['MAKE_CBCD_DIFFS.pro Ran '+systime(),'Difference between '+file_basename(cbcd_files[j]),' and '+file_basename(bcd_files[j])],cbcdh
      sxaddpar,cbcdh,'BITPIX',-64
      writefits,imask_2_files[j],imask,imskh
      writefits,diff_files[j],cbcd-bcd,cbcdh
    endfor
  endif
  
return
end
