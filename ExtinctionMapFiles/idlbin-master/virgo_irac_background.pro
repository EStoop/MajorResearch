pro virgo_irac_background

;want to determine the true zero background level.
for ch = 1, 1 do begin 

if ch eq 0 then begin
   fits_read, '/Users/jkrick/virgo/irac/same_fif/ch1_Combine-mosaic/mosaic.fits', data_raw, header
   fits_read, '/Users/jkrick/virgo/irac/same_fif/ch1_Combine-mosaic/mosaic_cov.fits', covdata, covheader


   a = where (finite(data_raw) eq 0)  ;data that are nan's
   data_nonan = data_raw
   data_nonan(a) = 0
   fits_write,  '/Users/jkrick/virgo/irac/same_fif/ch1_Combine-mosaic/mosaic_nonan.fits', data_nonan, header


   cd, '/Users/jkrick/virgo/irac/same_fif/ch1_Combine-mosaic'
   command = 'sex mosaic_nonan.fits -c ../../default.sex'
   spawn, command

   fits_read, '/Users/jkrick/Virgo/IRAC/same_fif/ch1_Combine-mosaic/segmentation.fits', segdata, segheader
endif

if ch eq 1 then begin
   fits_read, '/Users/jkrick/virgo/irac/same_fif/ch2/ch2_Combine-mosaic/mosaic.fits', data_raw, header
   fits_read, '/Users/jkrick/virgo/irac/same_fif/ch2/ch2_Combine-mosaic/mosaic_cov.fits', covdata, covheader
 
   a = where (finite(data_raw) eq 0)  ;data that are nan's
   data_nonan = data_raw

   data_nonan(a) = 0
   fits_write,  '/Users/jkrick/virgo/irac/same_fif/ch2/ch2_Combine-mosaic/mosaic_nonan.fits', data_nonan, header

   cd, '/Users/jkrick/virgo/irac/same_fif/ch2/ch2_Combine-mosaic'
   command = 'sex mosaic_nonan.fits -c ../../../default.sex'
   spawn, command

   fits_read, '/Users/jkrick/Virgo/IRAC/same_fif/ch2/ch2_Combine-mosaic/segmentation.fits', segdata, segheader
endif

;first use the segmentation image as a mask
data = data_raw
nan = alog10(-1)
segmask = where(segdata gt 0)
data(segmask) = nan

; a good big background region
;chosen based on optical image, ir image, and coverage map.  No
;saturated stars.

if ch eq 0 then begin
;raroi = [186.68976,186.64672,186.70705,186.75103]
;decroi =[13.572017,13.444642,13.428045,13.550216]
   raroi=[186.69731,186.65426,186.7146,186.75857]
   decroi  = [13.569359,13.441985,13.425386,13.547556]
endif
if ch eq 1 then begin
   raroi=[186.65561,186.63648,186.69681,186.78523,186.80576]
   decroi  = [13.467328,13.4053,13.388705,13.358774,13.415448]
endif


adxy, header, raroi, decroi, px2, py2

;roi photometry
naxis1 = fxpar(header, 'naxis1')
naxis2 = fxpar(header, 'naxis2')

photroi = Obj_New('IDLanROI', px2, py2) 
      
;create a mask out of the roi
mask = photroi->ComputeMask(Dimensions = [naxis1, naxis2], mask_rule = 1)
      
;apply the mask to make all data outside the roi = 0
onlyplume = data*(mask GT 0)
      
     
n = where(onlyplume eq 0)
onlyplume_nan = onlyplume
onlyplume_nan(n) = nan
      
c = where(covdata lt 25)
onlyplume_nan(c) = nan

a = where(onlyplume_nan  eq nan, counta)
goodarea = n_elements(where(onlyplume_nan gt 0 ))
goodfinite = n_elements(where(finite(onlyplume_nan) gt 0))
print, 'goodarea', goodarea, goodfinite, '    nan', counta
 

plumesb = total(onlyplume_nan,/nan)
plumesigma = stddev(onlyplume_nan, /nan)
print, 'background sb, background sigma   ', plumesb/goodarea, plumesigma


;ok, subtract
data_raw = data_raw - plumesb/goodarea
if ch eq 0 then fits_write, '/Users/jkrick/virgo/irac/same_fif/ch1_Combine-mosaic/mosaic_bkgd.fits', data_raw, header
if ch eq 1 then fits_write, '/Users/jkrick/virgo/irac/same_fif/ch2/ch2_Combine-mosaic/mosaic_bkgd.fits', data_raw, header

endfor

end
