pro pixnoise

  
  ra =269.84614
  dec = 66.049296
  c = 0

  dirloc = ['/Users/jkrick/nutella/spitzer/planets_proposal/pixnoise']
  cd, dirloc
    
  command1  =  "find . -name 'SPITZER_*bcd.fits' > /Users/jkrick/nutella/spitzer/planets_proposal/pixnoise/bcdlist.txt"
  spawn, command1
  
  readcol,'/Users/jkrick/nutella/spitzer/planets_proposal/pixnoise/bcdlist.txt', fitsname, format = 'A', /silent
  
  numoobjects = n_elements(fitsname)      ; 23 image per AOR, 2 channels, 47 AORs in april calstar pixnoise directory
  calstar = replicate({calob, ra:0D,dec:0D,bcdxcen:fltarr(64),bcdycen:fltarr(64),bcdflux:fltarr(64),ch:0U,filename:' ', sclktime:0D, aor:' '},numoobjects)

  for i =0, n_elements(fitsname) - 1 do begin
     print, 'working on ', fitsname(i)
     ;testname =  '/Users/jkrick/nutella/spitzer/planets_proposal/pixnoise/r39167488/ch2/bcd/SPITZER_I2_39167488_0023_0000_2_bcd.fits'
     header = headfits(fitsname(i)) ;
     sclk = sxpar(header, 'SCLK_OBS')
     chnl = sxpar(header, 'CHNLNUM')
     aornum = strmid(fitsname(i), 3, 8)
     ;first check to make sure the star is on the frame
     adxy, header, ra,dec , xcenter, ycenter
     if xcenter gt 0 and xcenter lt 30 and ycenter gt 0 and ycenter lt 30 then begin
        get_centroids,fitsname(i), t, dt, abcdxcen, abcdycen, abcdflux, xs, ys, fs, b, /WARM, /APER, APRAD = [5], RA=ra, DEC=dec ,/silent

        calstar[c] ={calob, ra, dec,abcdxcen, abcdycen, abcdflux,chnl, fitsname(i),sclk, aornum}
        
        print, 'test progress',calstar[c].ch, calstar[c].aor
        c = c + 1
     endif

  endfor
calstar = calstar[0:c-1]
save, calstar, filename='/Users/jkrick/nutella/spitzer/planets_proposal/pixnoise/pixnoise.sav'
;help, calstar

;a = where (calstar.aor eq 3911344)
;plot, calstar[a].bcdxcen - fix(calstar[a].bcdxcen), calstar[a].bcdycen - fix(calstar[a].bcdycen), xrange = [-0.8,0.8], yrange = [-0.6,0.8], psym = 1



undefine, calstar
end
