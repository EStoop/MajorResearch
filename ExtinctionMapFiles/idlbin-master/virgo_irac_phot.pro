pro virgo_irac_phot
;measure accurate photometry and background errors for IRAC data of
;ICL plume in the virgo cluster
;ps_open, filename='/Users/jkrick/Virgo/IRAC/background.ps',/portrait,/square,/color

;***Actually I do want to do this on at least a somewhat masked image.  Go for segmentation image as a first test.

!P.multi = [0,2,2]
greencolor = FSC_COLOR("Green", !D.Table_Size-4)
nan = alog10(-1)


araroi = [187.32094,187.27515,187.23817,187.31383,187.4071]
adecroi = [13.130524,13.084224,13.089371,12.87497,12.899801]

braroi = [187.13428,187.10259,187.27072,187.29272]
bdecroi = [13.079912,13.075608,12.843246,12.863825]

e1raroi = [187.40936,187.29532,187.31513,187.41925]
e1decroi = [12.880349,12.838504,12.761167,12.83201]

d1raroi = [186.88879,186.86213,186.85729,186.86701,186.88317,186.89125]
d1decroi =[13.22866,13.229458,13.190105,13.14446,13.136594,13.163355]

d2raroi=[186.88879,186.79582,186.80071,186.84598,186.85569,186.88317,186.86701,186.85729,186.86213]
d2decroi=[13.22866,13.254614,13.193228,13.168064,13.146031,13.136594,13.14446,13.190105,13.229458]

gasraroi = [186.54169,186.49866,186.56285,186.61082]
gasdecroi = [12.786601,12.682722,12.698795,12.733999]

for ch = 0, 1 do begin 
 offsb = fltarr(100)

;run SExtractor
   CD, strcompress('/Users/jkrick/Virgo/irac/ch' + string(ch+ 1) + '/ch' + string(ch+1) + '_Combine-mosaic',/remove_all)
   command = 'sex mosaic_bkgd.fits -c ../../default.sex'
   ;spawn, command

  if ch eq 0 then begin
      fits_read, '/Users/jkrick/Virgo/IRAC/ch1/ch1_Combine-mosaic/mosaic_bkgd.fits', data, header
      fits_read, '/Users/jkrick/Virgo//IRAC/ch1/ch1_Combine-mosaic/segmentation.fits', segdata, segheader
      print, '****ch1*****'
   endif

   if ch eq 1 then begin
      fits_read, '/Users/jkrick/Virgo/IRAC/ch2/ch2_Combine-mosaic/mosaic_bkgd.fits', data, header
      ;fits_read, '/Users/jkrick/Virgo/IRAC/ch2/ch2_Combine-mosaic/mosaic_cov.fits', covdata, covheader
      fits_read, '/Users/jkrick/Virgo//IRAC/ch2/ch2_Combine-mosaic/segmentation.fits', segdata, segheader
      print, '****ch2*****'
   endif


;first use the segmentation image as a mask
;   segmask = where(segdata gt 0)
;   data(segmask) = nan
;   fits_write, strcompress('/Users/jkrick/Virgo/IRAC/s18p14/ch' + string(ch + 1) + '/ch' + string(ch + 1) + '_mosaic_seg.fits',/remove_all), data, header

;try increasing the mask size of all objects uniformly
   readcol, 'test.cat', NUMBER                 , X_WORLD                , Y_WORLD                , X                , Y                , MAG_AUTO               , FLUX_BEST              , FLUXERR_BEST           , MAG_BEST               , MAGERR_BEST            , BACKGROUND             , FLUX_MAX               , ISOAREA_IMAGE          , ALPHA_J2000            , DELTA_J2000            , A_IMAGE                , B_IMAGE                , THETA_IMAGE            , MU_THRESHOLD           , MU_MAX                 , FLAGS                  , FWHM_IMAGE             , CLASS_STAR             , ELLIPTICITY            

   print, ellipticity(24)

   sma = sqrt(isoarea_image/(!PI*(1-ellipticity)))
   smb = isoarea_image/(!PI*sma)
   sma = 2.*sma                ; 3.0,2.3 ,1.6                     3.2;3.0*sma;2.6*sma;2*sma
   smb = 2.*smb                ; 3.0,2.3,1.6                     3.2;3.0*smb;2.6*smb;2*smb
   
;convert theta from degres into radians
   theta = theta_image*(!PI/180)
   increment =  (2*!PI)/500
  
   xmax = 6261                  ;1601;2049;1900;2866;1900;2866;2048;2000;2800;1000
   ymax = 6155                  ;1801;3148;2100;3486;2100;3486;3148;2000;2700;1000
   
   ;for each detection
   for d = long(0), n_elements(number) - 1 do begin

;calculate which pixels need to be masked
;to do this change to a new coordinate system (ang)
; 	then work arond the ellipse, each time
;	incrementing angle by some small amount
;	smaller increments for dimmer (smaller) objs.
      FOR ang = 0.0, (2*!PI),increment DO BEGIN
         xedge = sma(d)*cos(theta(d))*cos(ang) - smb(d)*sin(theta(d))*sin(ang)
         yedge = sma(d)*sin(theta(d))*cos(ang) + smb(d)*cos(theta(d))*sin(ang)
;   print, "x, y,xedge, yedge", x,y,xedge,yedge
         
                                ; first get the center pixel
         IF (x(d) GT 0 AND x(d) LT xmax) THEN BEGIN
            IF (y(d) GT 0 AND y(d) LT ymax) THEN BEGIN
               data[x(d),y(d)] = alog10(-1)
                                ;make sure the other pixels are not off the edge
                                ;then mask the little suckers
               IF (xedge+x(d) GT 0 AND xedge+x(d) LT xmax) THEN BEGIN
                  IF (yedge+y(d) GT 0 AND yedge+y(d) LT ymax) THEN BEGIN
                                ; have to do four cases to get everything 
                                ;between x and xegde and y and yedge
                     IF ((xedge+x(d)) LT x(d)) THEN BEGIN
                        If (yedge+y(d) LT y(d)) THEN BEGIN
                           data[xedge+x(d):x(d),yedge+y(d):y(d)] = alog10(-1)
                        ENDIF ELSE BEGIN
                           data[xedge+x(d):x(d),y(d):yedge+y(d)] = alog10(-1)
                        ENDELSE
                     ENDIF ELSE BEGIN
                        If (yedge+y(d) LT y(d)) THEN BEGIN
                           data[x(d):xedge+x(d),yedge+y(d):y(d)] = alog10(-1)
                        ENDIF ELSE BEGIN
                           data[x(d):xedge+x(d),y(d):yedge+y(d)] =alog10(-1)
                        ENDELSE
                     ENDELSE
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         
      ENDFOR                    ; end for each angle
   endfor                       ;end for each detection
   
   ;fits_write, strcompress('/Users/jkrick/Virgo/IRAC/s18p14/ch' + string(ch + 1) + '/ch' + string(ch + 1) + '_mosaic_myseg_2.fits',/remove_all), data, header
   
;roi photometry
   naxis1 = fxpar(header, 'naxis1')
   naxis2 = fxpar(header, 'naxis2')

   photlist = [ 'a','b', 'd1', 'd2', 'e1','gas']

   for r = 0, 5 do begin; n_elements(photlist) - 1 do begin
      case r of 
         0: begin
            raroi = araroi
            decroi = adecroi
         end
          1: begin
            raroi =braroi
            decroi = bdecroi
         end
         2: begin
            raroi = d1raroi
            decroi = d1decroi
         end
         3: begin
            raroi = d2raroi
            decroi = d2decroi
         end
         4: begin
            raroi = e1raroi
            decroi = e1decroi
         end
        5: begin
            raroi = gasraroi
            decroi = gasdecroi
         end
       endcase

      adxy, header, raroi, decroi, px2, py2
      
      photroi = Obj_New('IDLanROI', px2, py2) 
;   draw_roi,photroi             ;,/line_fill,color='444444'x 
      
;create a mask out of the roi
      mask = photroi->ComputeMask(Dimensions = [naxis1, naxis2], mask_rule = 1)
      
;apply the mask to make all data outside the roi = 0

      onlyplume = data*(mask GT 0)
      ;fits_write, '/Users/jkrick/Virgo/IRAC/s18p14/ch1/test_1.fits', onlyplume, header
                                ;try turning 0's into  Nan's.
      n = where(onlyplume eq 0)
      ;print, 'n  ', n_elements(n)
      onlyplume(n) = nan
      fits_write, strcompress('/Users/jkrick/Virgo/IRAC/ch'+ string(ch + 1)+'/onlyplume_2.fits',/remove_all), onlyplume, header
    
      a = where(finite(onlyplume) gt 0 , goodarea)
      a = where(finite(onlyplume)  lt 1, counta)
      print, 'goodarea', goodarea;, '       nan', counta
 
      plumesb = total(onlyplume,/nan)
      plumesigma = stddev(onlyplume, /nan)
      print,  photlist(r),  '  plumesb/goodarea, binnning length arcsec ',plumesb/goodarea, sqrt(goodarea) * 0.6
   endfor  ; end for r = 0, n_elements(_)


endfor ;end for each channel;
;ps_close, /noprint,/noid

end
;off_ch1 = [0.274932, 0.275366, 0.275689, 0.275391, 0.275030, 0.274892,0.274892,0.274942,0.275026,0.275089,0.274320]
;off_ch2 = [0.153731,0.153296, 0.154138, 0.154593, 0.154481, 0.154044,0.153741,0.153360,0.153192,0.153585,0.153432]


;-------------------

;aoff1raroi = [187.43426,187.42032,187.32534,187.41149,187.4775]
;aoff1decroi = [13.139112,13.171745,13.13138,12.900655,12.928049]

;aoff2raroi = [187.232,187.12153,187.0685,187.29536,187.30943]
;aoff2decroi = [13.091025,13.22886,13.185331,12.868053,12.873195]

;aoff3raroi = [187.14547,187.0603,187.12463,187.16636,187.26971]
;aoff3decroi = [13.012886,12.951521,12.867829,12.784827,12.840673]

;myaoff1raroi =[186.91182,186.79724,186.80218,186.91581]
;myaoff1decroi = [13.420772,13.417875,13.260217,13.222331]

;myaoff2raroi =[187.00689,186.91567,186.92141,186.99772,187.08518]
;myaoff2decroi = [13.379199,13.3917,13.210558,13.232393,13.304051]

;myaoff3raroi =[186.79638,186.70323,186.7084,186.82937,186.79579,186.79752]
;myaoff3decroi = [13.381629,13.382352,13.143142,13.144273,13.189525,13.260209]


;myaoff4raroi = [186.63038,186.63056,186.7916,186.77092]
;myaoff4decroi =[13.517481,13.387489,13.38565,13.53296]

;myaoff5raroi = [186.79286,186.77367,186.79705,186.917,186.917,186.99033,186.9828]
;myaoff5decroi =[13.538309,13.533629,13.422985,13.425036,13.395037,13.385041,13.489704]


;-----------------

;doff1raroi =[186.93606,186.88838,186.89289,186.9364]
;doff1decroi = [13.276319,13.276494,13.198924,13.198981]

;doff2raroi =[186.88581,186.75671,186.75698,186.88808]
;doff2decroi = [13.27661,13.279725,13.269017,13.230924]

;doff3raroi = [186.83216,186.7984,186.79413,186.7567,186.76451,186.82084]
;doff3decroi = [13.173966,13.19247,13.255522,13.266637,13.147909,13.146938]



;myoff1raroi = [186.94171,186.94022,186.97255,186.97434]
;myoff1decroi = [13.378422,13.22689,13.225427,13.37807]

;myoff2raroi = [186.77617,186.7777,186.93795,186.93972]
;myoff2decroi = [13.373961,13.343949,13.349141,13.378801]

;myoff3raroi = [186.77845,186.77697,186.93645,186.93747]
;myoff3decroi = [13.339557,13.313935,13.316199,13.344395]
