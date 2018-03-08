pro run_pld_snapshots, planetname, chname, apradius

  planetinfo = create_planetinfo()
  if chname eq '2' then aorname= planetinfo[planetname, 'aorname_ch2'] else aorname = planetinfo[planetname, 'aorname_ch1'] 
  basedir = planetinfo[planetname, 'basedir']
  dirname = strcompress(basedir + planetname +'/')
  savefilename = strcompress(dirname + planetname +'_phot_ch'+chname+'_'+string(apradius)+'.sav',/remove_all) ;
  print, 'restoring ', savefilename
  restore, savefilename
  
  startaor = 0
  stopaor =  n_elements(aorname) - 1
  for aor = startaor,stopaor do begin
     ;put together one big  array to fit Jim's code input
     if aor eq 0 then begin
        flux = planethash[aorname(aor),'flux'] 
        sigma_flux = planethash[aorname(aor),'fluxerr'] 
        pi = planethash[aorname(aor),'pixvals'] 
        s = size(pi)
        pixgrid = findgen(7,7,n_elements(aorname)*s(3)) 
        pixgrid[aor*s(3)] = pi
     
     endif else begin
        flux = [flux, planethash[aorname(aor),'flux'] ]
        sigma_flux =[sigma_flux,  planethash[aorname(aor),'fluxerr'] ]
        pi = planethash[aorname(aor),'pixvals'] 
        pixgrid[aor*s(3)] = pi
     endelse

  endfor

;get rid of NAN's  (add_pixvals already checks for nans in pixgrid
  notnan = where(finite(flux) gt 0 , goodcount, ncomplement = badnan)
  notnans =where( finite(sigma_flux) gt 0, goodcounts, ncomplement = badnans)

  print, 'n badnan', badnan, badnans, goodcount, goodcounts
  if badnan gt 0 or badnans gt 0 then begin
;XXXthis doesn't work for pixgrid, 
;need to figure out which dimension to cut?
     flux = flux[notnan]
     sigma_flux = sigma_flux[notnan]
     pixgrid = pixgrid[notnan]
  endif

  help, flux
  help, sigma_flux
  help, pixgrid

  corrected_flux = PLD_CORRECT_WITHFIT(pixgrid,flux,sigma_flux)

;when this starts working....
; can I add this to the planethash as a new key
;but they are no longer split up by aor.


end
