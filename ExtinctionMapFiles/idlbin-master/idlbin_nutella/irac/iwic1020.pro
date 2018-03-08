pro iwic1020

!P.multi = [0,2,2]
ps_open, filename='/Users/jkrick/IRAC/EGS/test/results.ps',/portrait,/square,/color

; to run this, download data into directory_100, etc.
;code will make a cdf directiry, put in your namelist, and run mopex

directory_name = '/Users/jkrick/IRAC/EGS/test'  ;parent directory of the bcd's
directory_100 = '/Users/jkrick/IRAC/EGS/test/test100'  ;parent directory of the bcd's
directory_200 = '/Users/jkrick/IRAC/EGS/test/test200'  ;parent directory of the bcd's
directory_400 = '/Users/jkrick/IRAC/EGS/test/test400'  ;parent directory of the bcd's
mopex_script_env = '~/bin/mopex/mopex-script-env.csh' ;where does mopex set environmental variables
;should check inside of the namelists that they are appropriate
;(pmasks, output directories, input file names)
namelist_ch1 = 'mosaic_ch1_namelist.nl'  ;favorite namelist ch1
namelist_ch2 = 'mosaic_ch2_namelist.nl'  ;favorite namelist ch2


;make cdf directory in each and copy namelists to that directory
dircommand1 = 'mkdir ' + directory_100 + '/cdf'
dircommand2 = 'mkdir ' + directory_200 + '/cdf'
dircommand3 = 'mkdir ' + directory_400 + '/cdf'
cpcommand1 = 'cp '+directory_name+'/'+namelist_ch1+'   ' +directory_100+'/cdf/'
cpcommand2 = 'cp '+directory_name+'/'+namelist_ch1+'   ' +directory_200+'/cdf/'
cpcommand3 = 'cp '+directory_name+'/'+namelist_ch1+'   ' +directory_400+'/cdf/'
cpcommand12 = 'cp '+directory_name+'/'+namelist_ch2+'   ' +directory_100+'/cdf/'
cpcommand22 = 'cp '+directory_name+'/'+namelist_ch2+'   ' +directory_200+'/cdf/'
cpcommand32 = 'cp '+directory_name+'/'+namelist_ch2+'   ' +directory_400+'/cdf/'

;spawn, dircommand1+'; '+dircommand2 + '; '+dircommand3 
;spawn, cpcommand1+'; '+cpcommand2 + '; '+cpcommand3 
;spawn, cpcommand12+'; '+cpcommand22 + '; '+cpcommand32 


;run mopex on the three datasets

all_dirs = [directory_100, directory_200, directory_400,'/Users/jkrick/IRAC/EGS/barmby', directory_100, directory_200, directory_400,'/Users/jkrick/IRAC/EGS/barmby']
for i = 0, n_elements(all_dirs) - 2 do begin
   junk = run_mopex(all_dirs(i), mopex_script_env, namelist_ch1)   ;ch1
   junk = run_mopex(all_dirs(i), mopex_script_env, namelist_ch2)   ;ch2
endfor


;re-name final mosaic names when reading them in
;--- ch1
fits_read, directory_100+'/mosaic_ch1/Combine/mosaic.fits',w_data100_ch1, w_header100_ch1
fits_read, directory_200+'/mosaic_ch1/Combine/mosaic.fits',w_data200_ch1, w_header200_ch1
fits_read, directory_400+'/mosaic_ch1/Combine/mosaic.fits',w_data400_ch1, w_header400_ch1
fits_read, directory_100+'/mosaic_ch1/Combine/mosaic_cov.fits',w_cov_data100_ch1, w_cov_header100_ch1
fits_read, directory_200+'/mosaic_ch1/Combine/mosaic_cov.fits',w_cov_data200_ch1, w_cov_header200_ch1
fits_read, directory_400+'/mosaic_ch1/Combine/mosaic_cov.fits',w_cov_data400_ch1, w_cov_header400_ch1
;--- ch2
fits_read, directory_100+'/mosaic_ch2/Combine/mosaic.fits',w_data100_ch2, w_header100_ch2
fits_read, directory_200+'/mosaic_ch2/Combine/mosaic.fits',w_data200_ch2, w_header200_ch2
fits_read, directory_400+'/mosaic_ch2/Combine/mosaic.fits',w_data400_ch2, w_header400_ch2
fits_read, directory_100+'/mosaic_ch2/Combine/mosaic_cov.fits',w_cov_data100_ch2, w_cov_header100_ch2
fits_read, directory_200+'/mosaic_ch2/Combine/mosaic_cov.fits',w_cov_data200_ch2, w_cov_header200_ch2
fits_read, directory_400+'/mosaic_ch2/Combine/mosaic_cov.fits',w_cov_data400_ch2, w_cov_header400_ch2

;also read in the Barmby data
fits_read, '/Users/jkrick/IRAC/EGS/barmby/egs_small_3.6.fits', barmbydata_ch1, barmbyheader_ch1
fits_read, '/Users/jkrick/IRAC/EGS/barmby/egs_small_cov_3.6.fits', barmbycovdata_ch1, barmbycovheader_ch1
;---ch2
fits_read, '/Users/jkrick/IRAC/EGS/barmby/egs_small_mos_4.5.fits', barmbydata_ch2, barmbyheader_ch2
fits_read, '/Users/jkrick/IRAC/EGS/barmby/egs_small_cov_4.5.fits', barmbycovdata_ch2, barmbycovheader_ch2

;change Nan's to zeros
;w_data100_ch1(where(~finite(w_data100_ch1))) = 0
;w_data200_ch1(where(~finite(w_data200_ch1))) = 0
;w_data400_ch1(where(~finite(w_data400_ch1))) = 0
;---
;w_data100_ch2(where(~finite(w_data100_ch2))) = 0
;w_data200_ch2(where(~finite(w_data200_ch2))) = 0
;w_data400_ch2(where(~finite(w_data400_ch2))) = 0

;make arrays of pointers to facilitate analyzing all data in same for
;loop.
;can't use structures here because the data arrays have
;different sizes and I don't know how to dynamically allocate
;data in structures.

;ch1data = (w_data100_ch1, w_data200_ch1, w_data400_ch1, barmbydata_ch1)
;ch1cov = (w_cov_data100_ch1, w_cov_data200_ch1, w_cov_data400_ch1, barbycovdata_ch1)
;ch1head = (w_header100_ch1, w_header200_ch1, w_header400_ch1, barmbyheader_ch1)

ch1dataptr = [Ptr_New(w_data100_ch1), Ptr_New(w_data200_ch1), Ptr_New(w_data400_ch1), Ptr_New(barmbydata_ch1)]
ch1covptr = [Ptr_New(w_cov_data100_ch1), Ptr_New(w_cov_data200_ch1), Ptr_New(w_cov_data400_ch1), Ptr_New(barmbycovdata_ch1)]
ch1headptr = [Ptr_New(w_header100_ch1), Ptr_New(w_header200_ch1), Ptr_New(w_header400_ch1), Ptr_New(barmbyheader_ch1)]
ch1covheadptr = [Ptr_New(w_cov_header100_ch1), Ptr_New(w_cov_header200_ch1), Ptr_New(w_cov_header400_ch1), Ptr_New(barmbycovheader_ch1)]

;---ch2
;ch2data = [w_data100_ch2, w_data200_ch2, w_data400_ch2, barmbyh_ch2]
;ch2cov = [w_cov_data100_ch2, w_cov_data200_ch2, w_cov_data400_ch2, barbycovdata_ch2]
;ch2head = [w_header100_ch2, w_header200_ch2, w_header400_ch2, barmbyheader_ch2]

ch2dataptr = [Ptr_New(w_data100_ch2), Ptr_New(w_data200_ch2), Ptr_New(w_data400_ch2), Ptr_New(barmbydata_ch2)]
ch2covptr = [Ptr_New(w_cov_data100_ch2), Ptr_New(w_cov_data200_ch2), Ptr_New(w_cov_data400_ch2), Ptr_New(barmbycovdata_ch2)]
ch2headptr = [Ptr_New(w_header100_ch2), Ptr_New(w_header200_ch2), Ptr_New(w_header400_ch2), Ptr_New(barmbyheader_ch2)]
ch2covheadptr = [Ptr_New(w_cov_header100_ch2), Ptr_New(w_cov_header200_ch2), Ptr_New(w_cov_header400_ch2), Ptr_New(barmbycovheader_ch2)]

;---
;alldata = [ch1data, ch2data]
;allcov = [ch1cov,ch2cov]
;allhead = [ch1head, ch2head]

alldataptr  = [Ptr_New(w_data100_ch1), Ptr_New(w_data200_ch1), Ptr_New(w_data400_ch1), Ptr_New(barmbydata_ch1), Ptr_New(w_data100_ch2), Ptr_New(w_data200_ch2), Ptr_New(w_data400_ch2), Ptr_New(barmbydata_ch2)]
allcovptr = [Ptr_New(w_cov_data100_ch1), Ptr_New(w_cov_data200_ch1), Ptr_New(w_cov_data400_ch1), Ptr_New(barmbycovdata_ch1), Ptr_New(w_cov_data100_ch2), Ptr_New(w_cov_data200_ch2), Ptr_New(w_cov_data400_ch2), Ptr_New(barmbycovdata_ch2)]
allheadptr = [Ptr_New(w_header100_ch1), Ptr_New(w_header200_ch1), Ptr_New(w_header400_ch1), Ptr_New(barmbyheader_ch1), Ptr_New(w_header100_ch2), Ptr_New(w_header200_ch2), Ptr_New(w_header400_ch2), Ptr_New(barmbyheader_ch2)]
allcovheadptr = [Ptr_New(w_cov_header100_ch1), Ptr_New(w_cov_header200_ch1), Ptr_New(w_cov_header400_ch1), Ptr_New(barmbycovheader_ch1),Ptr_New(w_cov_header100_ch2), Ptr_New(w_cov_header200_ch2), Ptr_New(w_cov_header400_ch2), Ptr_New(barmbycovheader_ch2)]


;alldataptr = [ch1dataptr, ch2dataptr]
;allcovptr = [ch1covptr, ch2covptr]
;allheadptr = [ch1headptr, ch2headptr]
;allcovheadptr = [ch1covheadptr, ch2covheadptr]

;--------------------------------------------------------------------------
;--------------------------------------------------------------------------
;--------------------------------------------------------------------------
; start analysis
;--------------------------------------------------------------------------
;
;median number of frames in coverage map
;XXX wait but the coverage header is different than the data header
medcov_ch1 = [fxpar(*ch1covheadptr(0), 'medcov'), fxpar(*ch1covheadptr(1), 'medcov'),fxpar(*ch1covheadptr(2), 'medcov'),17]
medcov_ch2 = [fxpar(*ch2covheadptr(0), 'medcov'), fxpar(*ch2covheadptr(1), 'medcov'),fxpar(*ch2covheadptr(2), 'medcov'),17]
all_medcov = [medcov_ch1, medcov_ch2]

;--------------------------------------------------------------------------
;--------------------------------------------------------------------------
;masure and compare background noise characteristics on 4 datasets per
;channel
;--------------------------------------------------------------------------
;start with working on just one channel
;XXX test this before switching to alldata , etc.
for i = 0, n_elements(alldataptr) - 1 do begin
   print, 'working on background noise ', i
;make an array of total flux values in each 36 pixel square region
   naxis1 = fxpar(*allheadptr[i], 'NAXIS1')
   naxis2 = fxpar(*allheadptr[i], 'NAXIS2')
   bkgobject = replicate({bkgob, xcenter:0D, ycenter:0D, total:0D, ra:0D, dec:0D, stddev:0D}, naxis1*naxis2)
   z = long(0)

   for x = 5, naxis1 - 5, 3 do begin
      for y = 5, naxis2 - 5, 3 do begin
;   for x = 500, 525, 3 do begin
;      for y = 500, 525, 3 do begin

         bool = 0

;get rid of blank areas, and bad combined pixels lt 0 counts anywhere
;in the aperture
         arr = (*alldataptr[i])[x-3:x+3, y-3:y+3]
 ;        if i eq 0 then print, 'arr', arr
         for anyval =0, n_elements(arr) - 1 do begin
            if arr(anyval) le 0 then bool = 1 ; no good
         endfor
;         print, 'bool after 0', bool
;only look at regions which are covered by more than the median number
;of frames.  take the info from the coverage map header
         if (*allcovptr[i])[x,y] lt all_medcov(i) then bool = 1 ;no good
;         print, 'bool after coverage', bool

         if bool lt 1 then begin   
            xyad, *allheadptr[i], x, y, ra, dec
            bkgobject[z].ra = ra
            bkgobject[z].dec = dec
            bkgobject[z].xcenter = x
            bkgobject[z].ycenter = y
            bkgobject[z].total = total((*alldataptr[i])[x-3:x+3,y-3:y+3])
            bkgobject[z].stddev = stddev((*alldataptr[i])[x-3:x+3,y-3:y+3])
            z = z + 1
         endif
      endfor
   endfor


   bkgobject = bkgobject[0:z-1]

;what is the limit where we are only looking at background areas?
   mmm, bkgobject.total, skymode, skysigma, skyskew

;diagnostic
;plot the distribution of total fluxes within all the apertures
;overplot mean +- sigma to see where the cutoff happens
   plothist, bkgobject.total, bin = 0.01, xrange=[-1,5], xtitle = 'Total counts in aperture', ytitle = 'number'
   oplot, fltarr(1000) + skymode, findgen(1000)
   oplot, fltarr(1000) + skymode+skysigma, findgen(1000)
   oplot, fltarr(1000) + skymode-skysigma, findgen(1000)
   bkglimit = skymode +skysigma

;diagnostic
;try plotting on ds9 the regions that meet this criterion
   openw, outlunred, all_dirs(i) + '/bkg.reg', /get_lun
   printf, outlunred, 'fk5'
   for c=long(0), z -2 do begin
      if bkgobject[c].total lt bkglimit then printf, outlunred, 'circle( ', bkgobject[c].ra, bkgobject[c].dec, ' 2")'
   endfor
   close, outlunred
   free_lun, outlunred

;what is the distribution of standard deviations inside of the
;background apertures
   good = where(bkgobject.total lt bkglimit)
   plothist, bkgobject[good].stddev, bin = 0.001, xtitle = 'Standard deviations within background apertures', ytitle = 'Number'
   
   print, 'mean stddev of background regions', mean(bkgobject[good].stddev)

endfor

;--------------------------------------------------------------------------------------
;--------------------------------------------------------------------------------------
;--------------------------------------------------------------------------------------
;now look at the objects, not just the background
;--------------------------------------------------------------------------------------

;XXXneed to think about how to put this together
;want to have in the end a plot of 100vsbarmby, 200vsbarmby,
;400vsbarmby for both channels.

;change this if necessary
gain = [3.3, 3.7]; ch1, ch2-4=, 3.7, 3.8, 3.8]

;only need to do the Barmby photometry once for each channel
find, barmbydata_ch1, xcen, ycen, flux, sharp, round, .01, 1.8, [0.0,1.0],[-1.0,1.0],/silent
xyad, barmbyheader_ch1, xcen, ycen, ch1_ra, ch1_dec
N= barmbycovdata_ch1[xcen, ycen]
phpadu = N*gain
aper,  barmbydata_ch1, xcen, ycen, ch1_mags, ch1_errap, sky, skyerr, phpadu, [3.3], [15,25], [-100,1000], /nan,/exact,/flux,/silent
good = where(barmbycovdata_ch1(xcen, ycen) gt 17)
ch1_ra = ch1_ra(good)
ch1_dec = ch1_dec(good)
ch1_mags = ch1_mags(good)
ch1_errap = ch1_errap(good)
;----
find, barmbydata_ch2, xcen, ycen, flux, sharp, round, .01, 1.8, [0.0,1.0],[-1.0,1.0],/silent
xyad, barmbyheader_ch2, xcen, ycen, ch2_ra, ch2_dec
N= barmbycovdata_ch2[xcen, ycen]
phpadu = N*gain
aper,  barmbydata_ch2, xcen, ycen, ch2_mags, ch2_errap, sky, skyerr, phpadu, [3.3], [15,25], [-100,1000], /nan,/exact,/flux,/silent
good = where(barmbycovdata_ch2(xcen, ycen) gt 17)
ch2_ra = ch2_ra(good)
ch2_dec = ch2_dec(good)
ch2_mags = ch2_mags(good)
ch2_errap = ch2_errap(good)

;----------------------------------------------------------------------------;
for i = 0, n_elements(alldataptr) - 1 do begin

   find, *alldataptr(i), warmxcen, warmycen, warmflux, warmsharp, warmround, .02, 1.8, [0.0,1.0],[-1.0,1.0],/silent
   xyad, *allheadptr(i), warmxcen, warmycen, warmra, warmdec


   N= (*allcovptr[i])[warmxcen, warmycen]
   warmphpadu = N*gain
   ;NaN = alog10(-1)
   ;a = where(*ch1dataptr[i] eq 0)
  ; (*ch1dataptr[i])[a] = NaN

   aper,  *alldataptr[i], warmxcen, warmycen, warmmags, warmerrap, warmsky, warmskyerr, warmphpadu, [3.3], [15,25], [-100,1000], /nan,/exact,/flux,/silent

   good = where((*allcovptr[i])[warmxcen, warmycen] gt 17)
   warmra = warmra(good)
   warmdec = warmdec(good)
   warmmags = warmmags(good)
   warmerrap = warmerrap(good)

;----------------------------------------------------------------------------;
;need to match the catalogs

   m=n_elements(warmra)
   mmatch=fltarr(m)
   mmatch[*]=-999
   dist=mmatch
   dist[*]=0

   matchedra = fltarr(n_elements(warmra))
   matcheddec = fltarr(n_elements(warmra))
   matchedwarmflux = fltarr(n_elements(warmra))
   matchedbarmbyflux = fltarr(n_elements(warmra))
   matchedwarmerrap = fltarr(n_elements(warmra))
   matchedbarmbyerrap = fltarr(n_elements(warmra))
   total = 0

;get the correct comparison catalog
   if i lt 3 then begin
      ra = ch1_ra
      dec = ch1_dec
      mags = ch1_mags
      errap = ch1_errap
   endif else begin
      ra = ch2_ra
      dec = ch2_dec
      mags = ch2_mags
      errap = ch2_errap
   endelse

   for q=0,m-1 do begin
      dist=sphdist( warmra(q), warmdec(q),ra,dec,/degrees)
      sep=min(dist,ind)
      if (sep LT 0.0008 and finite(warmmags(q)) gt 0 and finite(mags(ind)) gt 0)  then begin
         matchedra[total] = ra(ind)
         matcheddec[total]= dec(ind)
         matchedwarmflux[total] = warmmags(q)
         matchedbarmbyflux[total] = mags(ind)
         matchedbarmbyerrap[total] = errap(ind)
         matchedwarmerrap[total]= warmerrap(q)
         total = total + 1
      endif 
   endfor

   matchedra = matchedra[0:total - 2]
   matcheddec = matcheddec[0:total - 2]
   matchedwarmflux = matchedwarmflux[0:total - 2]
   matchedwarmerrap = matchedwarmerrap[0:total - 2]
   matchedbarmbyflux = matchedbarmbyflux[0:total - 2]
   matchedbarmbyerrap = matchedbarmbyerrap[0:total - 2]
   
   ;----------------------------------------------------------------------------;
;plot the photometry against each other
   plot, matchedwarmflux, matchedbarmbyflux, xtitle = 'warm photometry', ytitle = 'Barmby photometry', psym = 3, /xlog, /ylog
;oplot, findgen(10000) /10., findgen(10000)/10.
 ;---------------------------------------------------------------------------;
;measure the scatter in this relation

;try robust_lineft and fitexy both.  should be doing similar things,
;maybe robust_linefit has outlier rejections which is good but not
;errors on photometry.  fitexy knows about errors in both variables.
;go with that one, keep rest for reference.

;sort fluxes?
   sortbarmby = matchedbarmbyflux[sort(matchedwarmflux)]
   sortbarmbyerr = matchedbarmbyerrap[sort(matchedwarmflux)]
   sortwarmerr = matchedwarmerrap[sort(matchedwarmflux)]
   sortwarm = matchedwarmflux[sort(matchedwarmflux)]

;coeff1 = ROBUST_LINEFIT(sortwarm, sortbarmby, yfit, sig, coeff_sig)
;print, 'sig from robust', coeff1, sig


   fitexy, sortwarm, sortbarmby, intercept, slope, X_sigma =sortwarmerr, Y_sigma =sortbarmbyerr, sigma_A_B, chisq, q
   print, 'from fitexy: ', 'intercept ' , intercept, '    slope', slope, '  sigma(intercept)', sigma_A_B(0), '  sigma(slope)', sigma_A_B(1),'  chisq',  chisq, '  q',q

   oplot, sortwarm, slope*(sortwarm) + intercept

;six linear regression techniques from idlastro.  again they
;don't know about errors on the variables.
;sixlin, sortwarm, sortbarmby, a, siga, b, sigb
;print,'a', a ,'siga', siga, 'b',b,'sigb', sigb

endfor

;--------------------------------------------------------------------------
;--------------------------------------------------------------------------
;get a real measure of what the exposure time is for each frame time
;given that CR will be rejected and thereby decrease the exptime
;perpixel over a larger number of seconds in the longer exposure
;possible negating the advantage of more exptime per set time.
;--------------------------------------------------------------------------

;what is the effective exposure time of each frame time set?

;can measure the total exposure time over all pixels, and an average
;per pixel
;frame_time = [100,200,400,200]
frame_time = [100,200,400,200,100,200,400,200]

totalcov_100_ch1 = total(w_cov_data100_ch1) * 100.  ;in seconds
totalcov_200_ch1 = total(w_cov_data200_ch1) * 200.  ;in seconds
totalcov_400_ch1 = total(w_cov_data400_ch1) * 400.  ;in seconds
totalcov_barmby_ch1 = total(barmbycovdata_ch1)*200.
;---ch2
totalcov_100_ch2 = total(w_cov_data100_ch2) * 100.  ;in seconds
totalcov_200_ch2 = total(w_cov_data200_ch2) * 200. ;in seconds
totalcov_400_ch2 = total(w_cov_data400_ch2) * 400.  ;in seconds
totalcov_barmby_ch2 = total(barmbycovdata_ch2)*200.

print, 'total effective exptime, ch1; 100, 200, 400s, barmby(200) frame times ', $
       totalcov_100_ch1,totalcov_200_ch1,totalcov_400_ch1, totalcov_barmby_ch1
print, 'total effective exptime, ch2; 100, 200, 400s, barmby(200) frame times ', $
       totalcov_100_ch2,totalcov_200_ch2,totalcov_400_ch2, totalcov_barmby_ch2

;can also plot a cumulative distribution of exptimes/pixel.;
;don't know how easy it will be to read this, can get an
;average easily visually
;check just ch1 first.

for i = 0, n_elements(alldataptr) - 1 do begin
;need to first sort the coverage values
   sortcov = (*allcovptr[i])[sort(*allcovptr[i])]
   sortcov = sortcov(where(sortcov gt 5.)) ;don't really care about the mosaics with only a few images
   sortcov = sortcov*frame_time(i)
    N1 = n_elements(sortcov)
   f1 = (findgen(N1) + 1.)/ N1
   if i eq 0 then begin
      plot, sortcov, f1, xtitle = 'effective exposure time per pixel', ytitle='Cumulative Fraction', linestyle = i, yrange=[0,1], xrange=[100,30000]
   endif else begin
      oplot, sortcov, f1, linestyle = i
   endelse

   
endfor
ps_close, /noprint,/noid

end
