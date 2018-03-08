pro phot_stability, chname
;use BD+67_1044 to see what the photometric stability is on year-long timescales.
;will need to only look at one channel, one frametime, one narrow position.


t1 = systime(1)
dirname =  '/Users/jkrick/irac_warm/calstars/pmap_star_'+chname+'/'
if chname eq 'ch2' then begin
   ra_ref= 269.72731
   dec_ref = 67.793348
   ftime = 0.1
   starname = 'bd67_1044'
endif else begin
   ra_ref = 269.84561
   dec_ref = 66.048782
   ftime = 0.4
   starname = 'KF09T1'
endelse

planethash = hash()

;need a list of all relevant AORs.
cd,dirname
command = 'ls -d */ > aorlist.txt'
spawn, command
readcol, 'aorlist.txt', aorname, format = 'A', /silent

for a = 0, n_elements(aorname) -1 do begin
   undefine, fitsname
   undefine, buncname
   print, 'working on aorname ', aorname(a) , a, n_elements(aorname)
   cd, aorname(a)
   command1 = 'ls ch1/bcd/*_bcd.fits > ch1bcdnames.txt'
   command2 = 'ls ch2/bcd/*_bcd.fits > ch2bcdnames.txt'
   command3 =  'ls ch1/bcd/*_bunc.fits > ch1buncnames.txt'
   command4 =  'ls ch2/bcd/*_bunc.fits > ch2buncnames.txt'
   spawn, command1
   spawn, command2
   spawn, command3
   spawn, command4

   readcol, strcompress(chname + 'bcdnames.txt',/remove_all), fitsname, format = 'A', /silent, count = bcdcount
   readcol, strcompress(chname + 'buncnames.txt',/remove_all), buncname, format = 'A', /silent
   print, 'n count', bcdcount
   if bcdcount gt 0 then begin  ; got some correct channel data
      header = headfits(fitsname(0))        ;
      naxis = sxpar(header, 'NAXIS')
      framtime = sxpar(header, 'FRAMTIME')

      if naxis gt 2 and framtime eq ftime then begin  ; got some subarray data at the right exposure time
         fluxarr = fltarr(61*n_elements(fitsname))
         fluxerrarr = fluxarr
         backarr = fluxarr
         backerrarr = fluxarr
         nparr = fluxarr
         bmjdarr = fluxarr
         xarr = fluxarr
         yarr = fluxarr

      for i = 0, n_elements(fitsname) - 1 do begin
;         print, 'working on fitsname', i
         fits_read, fitsname(i), im, h
         fits_read, buncname(i), unc, hunc
         bmjd_obs= sxpar(header, 'BMJD_OBS')
         aintbeg = sxpar(header, 'AINTBEG')
         atimeend = sxpar(header, 'ATIMEEND')

         ;now not strictly correct, but I am going to bin anyway
         deltatime = (atimeend - aintbeg) / 61.D ; real time for each of the 64 frames
         nt = dindgen(61)
         bmjdarr[i*61]= bmjd_obs + (deltatime*nt)/60./60./24.D   ; 0.5*(frametime/60./60./24.) + (frametime/60./60./24.)*nt


         get_centroids_for_calstar_jk,im, h, unc, ra_ref, dec_ref,  t, dt, hjd, xft, x3, y3, $
                                      x5, y5, x7, y7, xg, yg, xh, yh, f, b, x3s, y3s, x5s, y5s, $
                                      x7s, y7s, fs, bs, xp3, yp3, xp5, yp5, xp7, yp7, xp3s, yp3s, $
                                      xp5s, yp5s, xp7s, yp7s, fp, fps, np, flag, ns, sf, $
                                      xfwhm, yfwhm, /WARM,/silent

         fluxarr[i*61] = f[3:*,1]      ; chose 2 pixel aperture size
         fluxerrarr[i*61] = fs[3:*,1]
         backarr[i*61] = b[3:*,0]
         backerrarr[i*61] = bs[3:*,0]
         nparr[i*61] = np[3:*]
         xarr[i*61] = x3[3:*]
         yarr[i*61] = y3[3:*]
       
      endfor  ; for each fits image

   endif  ; for subarray 0.1s data

   endif ; for the right channel data

   keys =['xarr', 'yarr', 'fluxarr','fluxerrarr',  'bmjdarr',  'backarr', 'backerrarr','np']
   values=list(xarr, yarr, fluxarr, fluxerrarr, bmjdarr,  backarr, backerrarr, nparr)
   planethash[aorname(a)] = HASH(keys, values)


   cd, dirname  ;reset for the next AOR
endfor  ; for each AOR
savename = strcompress(dirname + starname + '_phot_'+chname+'.sav',/remove_all)
save, planethash, filename=savename
print, 'saving planethash', savename
print, 'time check', systime(1) - t1




end
