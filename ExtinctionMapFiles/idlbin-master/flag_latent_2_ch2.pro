pro flag_latent_ch2

logdate=string(strmid(systime(0),8,2),format='(i2.2)')+strmid(systime(0),4,3)+strmid(systime(0),22,2)
print, logdate
logfile = '~/irac_warm/latent/flag_latent/idllog_ch2-'+logdate
if file_test(logfile) then print,'Log file found. no logging'
if not file_test(logfile) then print,'Beginning log file...' & journal,logfile

dirloc = '/Users/jkrick/irac_warm/latent/flag_latent/'


aorname = ['38763008_ch2','38929408_ch2','38933248_ch2','38929152_ch2']

;open file to contain master list of all the saturated stars
   mastername = dirloc  + '/master_satcorr_phot.tbl'
   openw,  masterlun, mastername,/get_lun
                                ;print a header for the log file
   printf, masterlun, '|  x    |  y    |    ra   |    dec  |replace|Pratio|Fratio| ferr | hsig | loc sky| skysig| radius|    raw  |   model  |        rectify | sclk | naxis_trig| duration'
   printf, masterlun, '|  -    |  -    |   deg   |    deg  | arcsec|  -   |  -   |  -   |  -   |  MJy/sr| MJy/sr| arcsec|    mJy  |   mJy    |        mJy     | s       | -  | s '
   ntot = 0

for a = 0,  n_elements(aorname) - 1 do begin  ;for each AOR


   CD,dirloc +aorname(a)   ; change directories to the correct AOR directory
   print, 'current dir ', dirloc + aorname(a)
   ;use unix to list all the filenames inside each AOR
   command1  =  "find . -name 'satcorr_phot*.tbl' > ./flag_latentlist.txt"
   spawn, command1
  ; *******original code used the bcds.
   command2  =  "find . -name 'cbcd*.fits' > ./cbcdlist.txt"
   spawn, command2

   ;read in all the individual satcorr_phot files, make them into one file with time stamps from sclk
   readcol,'./flag_latentlist.txt', latentlist, format = 'A', /silent
   readcol,'./cbcdlist.txt',cbcdlist, format = 'A', /silent
   
   print, 'n_elements(latentlist)', n_elements(latentlist)

   for i =0, n_elements(latentlist) - 1 do begin  ;read each file, output it to one common file

      ngood = file_lines(latentlist(i))   ;how many lines are in each file
      if ngood gt 0 then begin  ;if the file is not empty
         print, 'reading a saturated star'
         readcol,latentlist(i), x    ,  y    ,    ra   ,    dec  ,replace,Pratio,Fratio, ferr , hsig , loc_sky, skysig, radius,    raw  ,   $
                 model  ,  rectify , skipline = 3,/silent

         ;read in each cbcd header, also output the sclk time and the naxis so we
         ;know if it is subarray or full array
         header = headfits(cbcdlist(i))
         sclk = sxpar(header, 'SCLK_OBS')
         naxis_trig = sxpar(header, 'NAXIS')

         fits_read, cbcdlist(i), data, header
         
         for n = 0, ngood - 4 do begin  ;for each trigger in the individual file, output to master file
            print, 'ngood, ntot',  ngood, '  ', ntot, '  ', latentlist(i), x(n), y(n)
            senddata = data
            ;figure out how long to mask
            flux = quick_phot(x(n), y(n), senddata, 2,header)
                                ;known_flux_electron =
                                ;mjy_to_electron(rectify(lat(l)),
                                ;pixel_scale, gain, exptime,
                                ;flux_conv)
             print, 'flux', flux
             if flux le 1.96E7 then begin
                duration = 220
                print, 'flag for 220 seconds'
             endif
             if flux gt 1.96E7 and flux le 1.317E8 then begin
                duration = 320
                print, 'flag for 320 seconds'
             endif
             if flux gt 1.317E8 and flux le 3.74E8 then begin
                duration = 360
                print, 'flag for 360 seconds'
             endif
             if flux gt 3.74E8  then begin
                duration = 440
                print, 'flag for 440 seconds'
             endif
             
             print, 'duration', duration
             printf, masterlun, x(n)    ,  y(n)    ,    ra(n)   ,    dec(n)  ,replace(n),Pratio(n),Fratio(n), ferr (n), hsig(n) , $
                    loc_sky(n), skysig(n), radius(n),    raw(n) ,   model(n)  ,  rectify(n) , sclk, naxis_trig, duration,$
                    format = '(F10.2,F10.2,F10.5,F10.5,F10.2,F10.3,F10.3,F10.3,F10.3,F10.3,F10.3,F10.2,F10.2,F10.2,F10.2,D26.8, I10, I10)'
             
             
         endfor                 ;end for each trigger
         
      endif                     ; end if the file is not empty
      
   endfor                       ; number of images satcorr_phot lists

endfor                          ;end for each AOR in the PAO


close, masterlun
free_lun, masterlun

;-------------------------------
;now go through each image in each AOR and figure out which pixels need to be
;masked.

 ;read in the formatted list of triggers
fmt = '(F10.2,F10.2,F10.5,F10.5,F10.2,F10.3,F10.3,F10.3,F10.3,F10.3,F10.3,F10.2,F10.2,F10.2,F10.2,D26.8, I10, I10)'
readfmt, mastername, fmt, x    ,  y    ,    ra   ,    dec  ,replace,Pratio,Fratio, ferr , hsig , loc_sky, skysig, radius,    raw  ,   $
         model  ,  rectify , triggersclk, naxis_trig,duration, /silent

;replace is currently in units of arcseconds, change this to pixels
replace = replace / 1.2

for a = 0,  n_elements(aorname) - 1 do begin ;for each AOR
   CD,dirloc +aorname(a)                     ;change directories
   print, 'working on AOR', aorname(a)
   readcol, './cbcdlist.txt', cbcdlist, format = 'A'

   for j =0, n_elements(cbcdlist) - 1 do begin
      header = headfits(cbcdlist(j))
      currentsclk = sxpar(header, 'SCLK_OBS')
      naxis_trig = sxpar(header, 'NAXIS')
      
      ;in order to trigger a latent the current image has to be no more duration seconds  after the triggering event.
      lat = where(currentsclk gt triggersclk and currentsclk lt triggersclk + duration  , count)
;      print, cbcdlist(j), currentsclk, count
      for l = 0, count - 1 do begin ; start this loop if there is at least one trigger
                                ;        print, l, lat(l), cbcdlist(j), currentsclk, sclkarr(lat(l))
         ;print, 'flag ', cbcdlist(j), '  centered at x,y', x(lat(l)), y(lat(l))
         for x1 =fix(x(lat(l)) - replace(lat(l))) , fix(x(lat(l)) + replace(lat(l))) do begin
            for y1 = fix(y(lat(l)) - replace(lat(l))) , fix(y(lat(l)) + replace(lat(l))) do begin
               print, 'flag ', cbcdlist(j), x1, y1
            endfor
         endfor
         
      endfor  ; for each trigger
      
   endfor                       ; for each image in cbcdlist
   
endfor     ; for each AOR


journal

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function quick_phot, xcen, ycen, data, sq2, imheader
           ;get rid of NaN's
  f = where(finite(data) lt 1, fcount)
  print, 'n_finite',fcount

  if fcount gt 0 then data[where(finite(data) lt 1)] = 0
          ;convert the data from MJy/sr to electrons
  gain = sxpar(imheader, 'GAIN')
  gain = 3.7                    ;force it to the correct value since currently it is wrong in most headers
  exptime =  sxpar(imheader, 'EXPTIME') 
  fluxconv =  sxpar(imheader, 'FLUXCONV') 
  time = sxpar(imheader, 'SCLK_OBS')
  sbtoe = gain*exptime/fluxconv
  data = data*sbtoe
           
           ;do photometry in a 5x5 box around the object
  flux = total(data[xcen - sq2:xcen +sq2, ycen -sq2:ycen  + sq2])
  return, flux
end

