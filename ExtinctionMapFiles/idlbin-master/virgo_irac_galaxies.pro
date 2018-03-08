pro virgo_irac_galaxies
;measure IRAC colors of
;ICL plume in the virgo cluster
;ps_open, filename='/Users/jkrick/Virgo/IRAC/background.ps',/portrait,/square,/color


!P.multi = [0,1,1]
redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
greencolor = FSC_COLOR("Green", !D.Table_Size-4)
yellowcolor = FSC_COLOR("Yellow", !D.Table_Size-5)
cyancolor = FSC_COLOR("cyan", !D.Table_Size-6)
orangecolor = FSC_COLOR("orange", !D.Table_Size-7)
purplecolor = FSC_COLOR("purple", !D.Table_Size-8)
colorname = [yellowcolor, redcolor, greencolor, cyancolor, bluecolor, yellowcolor, orangecolor, purplecolor]
nan = alog10(-1)


araroi = [187.32094,187.27515,187.23817,187.31383,187.4071]
adecroi = [13.130524,13.084224,13.089371,12.87497,12.899801]

myaoff1raroi =[186.91182,186.79724,186.80218,186.91581]
myaoff1decroi = [13.420772,13.417875,13.260217,13.222331]

myaoff2raroi =[187.00689,186.91567,186.92141,186.99772,187.08518]
myaoff2decroi = [13.379199,13.3917,13.210558,13.232393,13.304051]

myaoff3raroi =[186.79638,186.70323,186.7084,186.82937,186.79579,186.79752]
myaoff3decroi = [13.381629,13.382352,13.143142,13.144273,13.189525,13.260209]


myaoff4raroi = [186.63038,186.63056,186.7916,186.77092]
myaoff4decroi =[13.517481,13.387489,13.38565,13.53296]

myaoff5raroi = [186.79286,186.77367,186.79705,186.917,186.917,186.99033,186.9828]
myaoff5decroi =[13.538309,13.533629,13.422985,13.425036,13.395037,13.385041,13.489704]

;and for 2 that are directly next to A
;myaoff6raroi = [187.19046,187.13334,187.23529,187.3074,187.23668]
;myaoff6decroi = [13.044475,13.029887,12.845322,12.868652,13.081086]

myaoff6raroi=[187.14547,187.0603,187.12463,187.16636,187.26971]
myaoff6decroi = [13.012886,12.951521,12.867829,12.784827,12.840673]

myaoff7raroi = [187.38011,187.33077,187.41559,187.50334]
myaoff7decroi = [13.211107,13.13061,12.900681,13.178655]

braroi = [187.13428,187.10259,187.27072,187.29272]
bdecroi = [13.079912,13.075608,12.843246,12.863825]

e1raroi = [187.40936,187.29532,187.31513,187.41925]
e1decroi = [12.880349,12.838504,12.761167,12.83201]

d1raroi = [186.88879,186.86213,186.85729,186.86701,186.88317,186.89125]
d1decroi =[13.22866,13.229458,13.190105,13.14446,13.136594,13.163355]

d2raroi=[186.88879,186.79582,186.80071,186.84598,186.85569,186.88317,186.86701,186.85729,186.86213]
d2decroi=[13.22866,13.254614,13.193228,13.168064,13.146031,13.136594,13.14446,13.190105,13.229458]


readcol, '/Users/jkrick/Virgo/irac/ch1/ch1_Combine-mosaic/test.cat', NUMBER_ch1                 , X_WORLD_ch1                , Y_WORLD_ch1                , X_ch1                , Y_ch1                , MAG_AUTO_ch1               , FLUX_BEST_ch1              , FLUXERR_BEST_ch1           , MAG_BEST_ch1               , MAGERR_BEST_ch1            , BACKGROUND_ch1             , FLUX_MAX_ch1               , ISOAREA_IMAGE_ch1          , ALPHA_J2000_ch1            , DELTA_J2000_ch1            , A_IMAGE_ch1                , B_IMAGE_ch1                , THETA_IMAGE_ch1            , MU_THRESHOLD_ch1           , MU_MAX_ch1                 , FLAGS_ch1                  , FWHM_IMAGE_ch1             , CLASS_STAR_ch1             , ELLIPTICITY_ch1  , format = "I10, D10.11, D10.11, F10.5, F10.5,      F10.5, F10.5,    F10.5, F10.5,    F10.5, F10.5,      F10.5, F10.5,    F10.5, F10.5, F10.5, F10.5,      F10.5, F10.5,    F10.5, I10, F10.5, F10.5, F10.5"

readcol, '/Users/jkrick/Virgo/irac/ch2/ch2_Combine-mosaic/test.cat', NUMBER_ch2                 , X_WORLD_ch2                , Y_WORLD_ch2                , X_ch2                , Y_ch2                , MAG_AUTO_ch2               , FLUX_BEST_ch2              , FLUXERR_BEST_ch2           , MAG_BEST_ch2               , MAGERR_BEST_ch2            , BACKGROUND_ch2             , FLUX_MAX_ch2               , ISOAREA_IMAGE_ch2          , ALPHA_J2000_ch2            , DELTA_J2000_ch2            , A_IMAGE_ch2                , B_IMAGE_ch2                , THETA_IMAGE_ch2            , MU_THRESHOLD_ch2           , MU_MAX_ch2                 , FLAGS_ch2                  , FWHM_IMAGE_ch2             , CLASS_STAR_ch2             , ELLIPTICITY_ch2       , format = "I10, D10.11, D10.11, F10.5, F10.5,      F10.5, F10.5,    F10.5, F10.5,    F10.5, F10.5,      F10.5, F10.5,    F10.5, F10.5, F10.5, F10.5,      F10.5, F10.5,    F10.5, I10, F10.5, F10.5, F10.5"      


print, X_world_ch1(35)

for r = 0, 6 do begin           ; n_elements(photlist) - 1 do begin
   case r of 
      0: begin
            raroi = myaoff3raroi
            decroi = myaoff3decroi
         end
      1: begin
         raroi =myaoff1raroi
         decroi = myaoff1decroi
      end
      2: begin
         raroi = myaoff2raroi
         decroi = myaoff2decroi
      end
      3: begin
         raroi = myaoff3raroi
         decroi = myaoff3decroi
      end
      4: begin
         raroi = myaoff4raroi
         decroi = myaoff4decroi
      end
     5: begin
         raroi = myaoff5raroi
         decroi = myaoff5decroi
      end
     6: begin
         raroi = myaoff6raroi
         decroi = myaoff6decroi
      end
     7: begin
         raroi = myaoff7raroi
         decroi = myaoff7decroi
      end
   endcase
   
   
   
   for ch = 0, 1 do begin 
      offsb = fltarr(100)
      
      CD, strcompress('/Users/jkrick/Virgo/irac/ch' + string(ch+ 1) + '/ch' + string(ch+1) + '_Combine-mosaic',/remove_all)
      
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
      
;roi photometry
      naxis1 = fxpar(header, 'naxis1')
      naxis2 = fxpar(header, 'naxis2')
      
      photlist = [ 'a','b', 'd1', 'd2', 'e1']
      
      
      adxy, header, raroi, decroi, px2, py2
      
      photroi = Obj_New('IDLanROI', px2, py2) 
      
      if ch eq 0 then ptTest = photroi -> ContainsPoints(x_ch1,y_ch1) 
      if ch eq 1 then ptTest = photroi -> ContainsPoints(x_ch2,y_ch2) 
      
      containResults = [ $ 
                       'Point lies outside ROI', $ 
                       'Point lies inside ROI', $ 
                       'Point lies on the edge of the ROI', $ 
                       'Point lies on vertex of the ROI'] 
      
      PRINT, 'Result =',ptTest[46],':   ', containResults[ptTest[46]] 
      
      if ch eq 0 then inside_ch1 = where(ptTest eq 1, incount)
      if ch eq 1 then inside_ch2 = where(ptTest eq 1, incount)
      
      
      
   endfor                       ;end for each channel;
   
   ra_ch1 = x_world_ch1(inside_ch1)
   dec_ch1 = y_world_ch1(inside_ch1)
   ra_ch2 = x_world_ch1(inside_ch2)
   dec_ch2 = y_world_ch1(inside_ch2)
   
   
; now match ch1 and ch2 catalog of sources which are inside the plume:
   
; match ch1 to ch2
   i36=n_elements(x_ch1(inside_ch1))
   i36match=fltarr(i36)
   i36match[*]=-999
   print,'Matching '
   print,"Starting at "+systime()
   dist=i36match
   dist[*]=0
   colorarr = fltarr(i36)
   p = 0
   print, 'n_insidech1', n_elements(inside_ch1)
   
   for q=0,i36-1 do begin
      
      dist=sphdist(ra_ch1[q],dec_ch1[q],ra_ch2,dec_ch2, /degrees)
      sep=min(dist,ind)
      if (sep LT 0.0008) then begin
         i36match[q]=ind
                                ;iracobject[ind].iracmatch = q
                                ;iracobject[ind].iracmatchdist = sep
;         print, 'match', ra_ch1[q], dec_ch1[q], ra_ch2[ind], dec_ch2[ind]
         colorarr[p] = mag_best_ch1[inside_ch1[q]] - mag_best_ch2[inside_ch2[ind]]
         p = p + 1
      endif
      if (sep ge 0.0008 and sep lt .001) then begin
         print, 'not a match', ra_ch1[q], dec_ch1[q], ra_ch2[ind], dec_ch2[ind]
      endif

   endfor
   
   print,"Finished at "+systime()
   matched=where(i36match GE 0)
   nonmatched = where(i36match lt 0)
   print, n_elements(matched),"matched"
   print, n_elements(nonmatched),"nonmatched"
   
   
;make a histogram of irac colors inside each plume, and then also in
;the background regions.
   colorarr = colorarr[0:p-1]
   
   if r eq 0 then begin
      plothist, colorarr , bin = 0.5, peak = 1, thick = 2
      color_a = colorarr
   endif else begin
      plothist, colorarr , bin = 0.5,/overplot, peak = 1, color = colorname[r]
      color_aoff = colorarr

;are these populations drawn from the same distribution?
      kstwo, color_a, color_aoff, ksdist, prob
      print, '-------------------------------------------------------------------------------'
      print, 'ksdist for a with background, prob', ksdist, prob
      
   endelse
   

endfor



;ps_close, /noprint,/noid
end
