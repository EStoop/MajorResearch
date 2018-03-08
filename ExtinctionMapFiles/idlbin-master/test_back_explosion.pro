pro test_back_explosion
;  ra_ref = 133.14754
;  dec_ref = 28.330195
  fwhm = 1.415

  ;ch2
  aorname = ['pcrs_planets/hd7924/r44605184','pcrs_planets/55cnc/r43981056', 'pcrs_planets/55cnc/r43981312', 'pcrs_planets/55cnc/r39524608', 'pcrs_planets/55cnc/r43981568', 'pcrs_planets/55cnc/r43981824', 'staring/r34769408', 'staring/r34869504', 'staring/r42000384', 'staring/r42790912', 'staring/r44273920']
  binsize = [64, 64, 100*64L,64L, 64, 64, 7, 7]

  ;ch1
;  aorname = ['pcrs_planets/wasp33/r39445760', 'pcrs_planets/wasp33/r45383680', 'pcrs_planets/wasp33/r45383936', 'pcrs_planets/wasp33/r45384192','staring/r34776832', 'staring/r38110464', 'staring/r34769152', 'staring/r38110720']
;  binsize = [7,7,7,7,7,7,64, 64]

  for a = 0,  0 do begin; n_elements(aorname) - 1 do begin
     print, 'working on aorname ', aorname(a)
     dir = '/Users/jkrick/irac_warm/' + aorname(a) + '/ch2/bcd/'
     CD, dir                    ; change directories to the corrct AOR directory
     command  =  "find . -name '*_bcd.fits' > /Users/jkrick/irac_warm/staring/bcdlist.txt"
     spawn, command
     readcol,'/Users/jkrick/irac_warm/staring/bcdlist.txt',fitsname, format = 'A', /silent

     xarr = dblarr(64*n_elements(fitsname))
     yarr = dblarr(64*n_elements(fitsname))
     skyarr = dblarr(64*n_elements(fitsname))
     fluxarr = dblarr(64*n_elements(fitsname))
     backarr = dblarr(64*n_elements(fitsname))
     backsigarr = dblarr(64*n_elements(fitsname))
     sclkarr = dblarr(64*n_elements(fitsname))

     VRSTUCCarr=dblarr(64*n_elements(fitsname))
     AVRSTBEGarr=dblarr(64*n_elements(fitsname))
     AVDETCarr= dblarr(64*n_elements(fitsname))
     AVDETBEGarr=dblarr(64*n_elements(fitsname))
     AVGG1Carr= dblarr(64*n_elements(fitsname))
     AVGG1BEGarr=dblarr(64*n_elements(fitsname))
     AVDDUCCarr= dblarr(64*n_elements(fitsname))
     AVDDUBEGarr=dblarr(64*n_elements(fitsname))
     AVGGCLCarr= dblarr(64*n_elements(fitsname))
     AVGGCBEGarr=dblarr(64*n_elements(fitsname))
     AHTRIBEGarr=dblarr(64*n_elements(fitsname))
     AHTRVBEGarr=dblarr(64*n_elements(fitsname))
     AFPAT2Barr=dblarr(64*n_elements(fitsname))
     AFPAT2BTarr=dblarr(64*n_elements(fitsname))
     AFPAT2Earr=dblarr(64*n_elements(fitsname))
     AFPAT2ETarr=dblarr(64*n_elements(fitsname))
     ACTENDTarr=dblarr(64*n_elements(fitsname))
     AFPECTEarr=dblarr(64*n_elements(fitsname))
     AFPEATEarr=dblarr(64*n_elements(fitsname))
     ASHTEMPEarr=dblarr(64*n_elements(fitsname))
     ATCTEMPEarr=dblarr(64*n_elements(fitsname))
     ACETEMPEarr=dblarr(64*n_elements(fitsname))
     APDTEMPEarr=dblarr(64*n_elements(fitsname))
     ACATMP1Earr=dblarr(64*n_elements(fitsname))
     ACATMP2Earr=dblarr(64*n_elements(fitsname))
     ACATMP3Earr=dblarr(64*n_elements(fitsname))
     ACATMP4Earr=dblarr(64*n_elements(fitsname))
     ACATMP5Earr=dblarr(64*n_elements(fitsname))
     ACATMP6Earr=dblarr(64*n_elements(fitsname))
     ACATMP7Earr=dblarr(64*n_elements(fitsname))
     ACATMP8Earr=dblarr(64*n_elements(fitsname))


     i = 0L
     print, 'n fits', n_elements(fitsname)
     ap3 = [3.]
     back3 = [11.,15.]           ; [11., 15.5]
     
     for f =0.D,  n_elements(fitsname) - 1 do begin ;read each cbcd file, find centroid, keep track
        fits_read, fitsname(f), data, header  
        gain = sxpar(header, 'GAIN')
        fluxconv = sxpar(header, 'FLUXCONV')
        exptime = sxpar(header, 'EXPTIME')
        ronoise = sxpar(header, 'RONOISE')
        aorlabel = sxpar(header, 'OBJECT')
        sclk = sxpar(header, 'SCLK_OBS')

        ;pull a bunch of other telemetry
        VRSTUCC= sxpar(header, 'VRSTUCC')
        AVRSTBEG= sxpar(header, 'AVRSTBEG')
        AVDETC = sxpar(header, 'AVDETC')
        AVDETBEG= sxpar(header, 'AVDETBEG')
        AVGG1C = sxpar(header, 'AVGG1C')
        AVGG1BEG= sxpar(header, 'AVGG1BEG')
        AVDDUCC = sxpar(header, 'AVDDUCC')
        AVDDUBEG= sxpar(header, 'AVDDUBEG')
        AVGGCLC = sxpar(header, 'AVGGCLC')
        AVGGCBEG= sxpar(header, 'AVGGCBEG')
        AHTRIBEG= sxpar(header, 'AHTRIBEG')
        AHTRVBEG= sxpar(header, 'AHTRVBEG')
        AFPAT2B=  sxpar(header, 'AFPAT2B')
        AFPAT2BT= sxpar(header, 'AFPAT2BT')
        AFPAT2E=  sxpar(header, 'AFPAT2E')
        AFPAT2ET= sxpar(header, 'AFPAT2ET')
        ACTENDT=  sxpar(header, 'ACTENDT')
        AFPECTE = sxpar(header, 'AFPECTE')
        AFPEATE = sxpar(header, 'AFPEATE')
        ASHTEMPE= sxpar(header, 'ASHTEMPE')
        ATCTEMPE= sxpar(header, 'ATCTEMPE')
        ACETEMPE= sxpar(header, 'ACETEMPE')
        APDTEMPE= sxpar(header, 'APDTEMPE')
        ACATMP1E= sxpar(header, 'ACATMP1E')
        ACATMP2E= sxpar(header, 'ACATMP2E')
        ACATMP3E= sxpar(header, 'ACATMP3E')
        ACATMP4E= sxpar(header, 'ACATMP4E')
        ACATMP5E= sxpar(header, 'ACATMP5E')
        ACATMP6E= sxpar(header, 'ACATMP6E')
        ACATMP7E= sxpar(header, 'ACATMP7E')
        ACATMP8E= sxpar(header, 'ACATMP8E')


        convfac = gain*exptime/fluxconv
;        adxy, header, ra_ref, dec_ref, xinit, yinit
        xinit = 15.5
        yinit = 15.5
        for j = 0, 63 do begin
           indim = data[*,*,j]
           indim = indim*convfac
           irac_box_centroider, indim, xinit, yinit, 3, 6, 3, ronoise,xcen, ycen, box_f, box_sky, box_c, box_cb, box_np,/MMM
           xarr[i] = xcen
           yarr[i] = ycen
           skyarr[i] = box_sky
           
;now run aper at the centroids found, and use that background estimate.
           aper, indim,xcen, ycen, xf, xfs, xb, xbs, 1.0, ap3, back3, $
                 /FLUX, /EXACT, /SILENT, /NAN, READNOISE=ronoise
           fluxarr[i] = xf
           backarr[i] = xb
           backsigarr[i] =xbs
           sclkarr[i] = sclk

           VRSTUCCarr[i]=VRSTUCC
           AVRSTBEGarr[i]=AVRSTBEG
           AVDETCarr[i] =AVDETC
           AVDETBEGarr[i]=AVDETBEG
           AVGG1Carr[i] =AVGG1C
           AVGG1BEGarr[i]=AVGG1BEG
           AVDDUCCarr[i] =AVDDUCC
           AVDDUBEGarr[i]=AVDDUBEG
           AVGGCLCarr[i] =AVGGCLC
           AVGGCBEGarr[i]=AVGGCBEG
           AHTRIBEGarr[i]=AHTRIBEG
           AHTRVBEGarr[i]=AHTRVBEG
           AFPAT2Barr[i]=AFPAT2B
           AFPAT2BTarr[i]=AFPAT2BT
           AFPAT2Earr[i]=AFPAT2E
           AFPAT2ETarr[i]=AFPAT2ET
           ACTENDTarr[i]=ACTENDT
           AFPECTEarr[i]=AFPECTE
           AFPEATEarr[i]=AFPEATE
           ASHTEMPEarr[i]=ASHTEMPE
           ATCTEMPEarr[i]=ATCTEMPE
           ACETEMPEarr[i]=ACETEMPE
           APDTEMPEarr[i]=APDTEMPE
           ACATMP1Earr[i]=ACATMP1E
           ACATMP2Earr[i]=ACATMP2E
           ACATMP3Earr[i]=ACATMP3E
           ACATMP4Earr[i]=ACATMP4E
           ACATMP5Earr[i]=ACATMP5E
           ACATMP6Earr[i]=ACATMP6E
           ACATMP7Earr[i]=ACATMP7E
           ACATMP8Earr[i]=ACATMP8E


           i = i + 1

        endfor
     endfor
     print, 'i', i
     xarr = xarr[0:i-1]
     yarr = yarr[0:i-1]
     skyarr = skyarr[0:i-1]
     fluxarr = fluxarr[0:i-1]
     backarr = backarr[0:i-1]
     backsigarr = backsigarr[0:i-1]
     sclkarr = sclkarr[0:i-1]
     numberarr = findgen(n_elements(xarr))
     
     VRSTUCCarr = VRSTUCCarr[0:i-1]
     AVRSTBEGarr = AVRSTBEGarr[0:i-1]
     AVDETCarr = AVDETCarr[0:i-1]
     AVDETBEGarr = AVDETBEGarr[0:i-1]
     AVGG1Carr = AVGG1Carr[0:i-1]
     AVGG1BEGarr = AVGG1BEGarr[0:i-1]
     AVDDUCCarr = AVDDUCCarr[0:i-1]
     AVDDUBEGarr = AVDDUBEGarr[0:i-1]
     AVGGCLCarr = AVGGCLCarr[0:i-1]
     AVGGCBEGarr = AVGGCBEGarr[0:i-1]
     AHTRIBEGarr = AHTRIBEGarr[0:i-1]
     AHTRVBEGarr = AHTRVBEGarr[0:i-1]
     AFPAT2Barr = AFPAT2Barr[0:i-1]
     AFPAT2BTarr = AFPAT2BTarr[0:i-1]
     AFPAT2Earr = AFPAT2Earr[0:i-1]
     AFPAT2ETarr = AFPAT2ETarr[0:i-1]
     ACTENDTarr = ACTENDTarr[0:i-1]
     AFPECTEarr = AFPECTEarr[0:i-1]
     AFPEATEarr = AFPEATEarr[0:i-1]
     ASHTEMPEarr = ASHTEMPEarr[0:i-1]
     ATCTEMPEarr = ATCTEMPEarr[0:i-1]
     ACETEMPEarr = ACETEMPEarr[0:i-1]
     APDTEMPEarr = APDTEMPEarr[0:i-1]
     ACATMP1Earr = ACATMP1Earr[0:i-1]
     ACATMP2Earr = ACATMP2Earr[0:i-1]
     ACATMP3Earr = ACATMP3Earr[0:i-1]
     ACATMP4Earr = ACATMP4Earr[0:i-1]
     ACATMP5Earr = ACATMP5Earr[0:i-1]
     ACATMP6Earr = ACATMP6Earr[0:i-1]
     ACATMP7Earr = ACATMP7Earr[0:i-1]
     ACATMP8Earr = ACATMP8Earr[0:i-1]

     
;  a = plot(numberarr, xarr,'1.', xtitle = 'number', ytitle = 'x cntrd', yrange = [14.4, 14.9])
;  a = plot(numberarr, yarr,'1.', xtitle = 'number', ytitle = 'y cntrd', yrange = [14.9, 15.1])
;     a = plot(numberarr, skyarr, '1.', xtitle = 'number', ytitle = 'box_sky')
;  a = plot(numberarr, fluxarr,  '1.', xtitle = 'number', ytitle = 'aper_flux')
;     a = plot(numberarr, backarr,  '1.', xtitle = 'number', ytitle = 'aper_bkgd')
     
;----------------------------
;try binning
;histogram numberarr by 64 as a technique to bin
     h = histogram(numberarr, OMIN=om, binsize = binsize(a), reverse_indices = ri)
                                ;print, 'h', h
     print, 'omin', om, 'nh', n_elements(h)
     
;mean together the flux values in each phase bin
     binxcen = dblarr(n_elements(h))
     binnum = dblarr(n_elements(h))
     binycen= dblarr(n_elements(h))
     binbkgd =  dblarr(n_elements(h))
     binflux=  dblarr(n_elements(h))
     binback=  dblarr(n_elements(h))
     bintime = dblarr(n_elements(h))
     binbacksig = dblarr(n_elements(h))

     binVRSTUCCarr = dblarr(n_elements(h))
     binAVRSTBEGarr = dblarr(n_elements(h))
     binAVDETCarr = dblarr(n_elements(h))
     binAVDETBEGarr = dblarr(n_elements(h))
     binAVGG1Carr = dblarr(n_elements(h))
     binAVGG1BEGarr = dblarr(n_elements(h))
     binAVDDUCCarr = dblarr(n_elements(h))
     binAVDDUBEGarr = dblarr(n_elements(h))
     binAVGGCLCarr = dblarr(n_elements(h))
     binAVGGCBEGarr = dblarr(n_elements(h))
     binAHTRIBEGarr = dblarr(n_elements(h))
     binAHTRVBEGarr = dblarr(n_elements(h))
     binAFPAT2Barr = dblarr(n_elements(h))
     binAFPAT2BTarr = dblarr(n_elements(h))
     binAFPAT2Earr = dblarr(n_elements(h))
     binAFPAT2ETarr = dblarr(n_elements(h))
     binACTENDTarr = dblarr(n_elements(h))
     binAFPECTEarr = dblarr(n_elements(h))
     binAFPEATEarr = dblarr(n_elements(h))
     binASHTEMPEarr = dblarr(n_elements(h))
     binATCTEMPEarr = dblarr(n_elements(h))
     binACETEMPEarr = dblarr(n_elements(h))
     binAPDTEMPEarr = dblarr(n_elements(h))
     binACATMP1Earr = dblarr(n_elements(h))
     binACATMP2Earr = dblarr(n_elements(h))
     binACATMP3Earr = dblarr(n_elements(h))
     binACATMP4Earr = dblarr(n_elements(h))
     binACATMP5Earr = dblarr(n_elements(h))
     binACATMP6Earr = dblarr(n_elements(h))
     binACATMP7Earr = dblarr(n_elements(h))
     binACATMP8Earr = dblarr(n_elements(h))


     c = 0L
     for j = 0L, n_elements(h) - 2 do begin
        
;get rid of the bins with no values and low numbers, meaning low overlap
        if (ri[j+1] gt ri[j])  then begin
                                ;print, 'binning together', n_elements(numberarr[ri[ri[j]:ri[j+1]-1]])
                                ;print, 'binning', numberarr[ri[ri[j]:ri[j+1]-1]]
           
;           meanclip, xarr[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
;           binxcen[c] = meanx   ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
           
;           meanclip, yarr[ri[ri[j]:ri[j+1]-1]], meany, sigmay
;           binycen[c] = meany   ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
           
;           meanclip, skyarr[ri[ri[j]:ri[j+1]-1]], meansky, sigmasky
;           binbkgd[c] = meansky ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
           
           binnum[c] = mean(numberarr[ri[ri[j]:ri[j+1]-1]])
           
           meanclip, fluxarr[ri[ri[j]:ri[j+1]-1]], meanflux, sigmaflux
           binflux[c] = meanflux ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
           
           meanclip, backarr[ri[ri[j]:ri[j+1]-1]], meanback, sigmaback
           mmm, backarr[ri[ri[j]:ri[j+1]-1]],skymod, skysigma, skew, /silent
           ;print, 'meanclip, mmm', sigmaback, skysigma
           binback[c] = skymod; meanback ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])

           meanclip, backsigarr[ri[ri[j]:ri[j+1]-1]], meanbacksig, sigmabacksig

           binbacksig[c]= skysigma; sigmaback

           bintime[c] =mean(sclkarr[ri[ri[j]:ri[j+1]-1]])

           ;plot some histograms
          ; plothist,backarr[ri[ri[j]:ri[j+1]-1]], xhist, yhist, bin = 0.4, /noplot
          ; if c eq 0 then begin
          ;    sp = plot(xhist, yhist, xtitle = 'background (electrons)') 
          ; endif else begin
          ;    sp = plot(xhist, yhist - j*50., /overplot)
         ;  endelse

           binVRSTUCCarr[c] = mean(VRSTUCCarr[ri[ri[j]:ri[j+1]-1]])
           binAVRSTBEGarr[c] = mean(AVRSTBEGarr[ri[ri[j]:ri[j+1]-1]])
           binAVDETCarr[c] = mean(AVDETCarr[ri[ri[j]:ri[j+1]-1]])
           binAVDETBEGarr[c] = mean(AVDETBEGarr[ri[ri[j]:ri[j+1]-1]])
           binAVGG1Carr[c] = mean(AVGG1Carr[ri[ri[j]:ri[j+1]-1]])
           binAVGG1BEGarr[c] = mean(AVGG1BEGarr[ri[ri[j]:ri[j+1]-1]])
           binAVDDUCCarr[c] = mean(AVDDUCCarr[ri[ri[j]:ri[j+1]-1]])
           binAVDDUBEGarr[c] = mean(AVDDUBEGarr[ri[ri[j]:ri[j+1]-1]])
           binAVGGCLCarr[c] = mean(AVGGCLCarr[ri[ri[j]:ri[j+1]-1]])
           binAVGGCBEGarr[c] = mean(AVGGCBEGarr[ri[ri[j]:ri[j+1]-1]])
           binAHTRIBEGarr[c] = mean(AHTRIBEGarr[ri[ri[j]:ri[j+1]-1]])
           binAHTRVBEGarr[c] = mean(AHTRVBEGarr[ri[ri[j]:ri[j+1]-1]])
           binAFPAT2Barr[c] = mean(AFPAT2Barr[ri[ri[j]:ri[j+1]-1]])
           binAFPAT2BTarr[c] = mean(AFPAT2BTarr[ri[ri[j]:ri[j+1]-1]])
           binAFPAT2Earr[c] = mean(AFPAT2Earr[ri[ri[j]:ri[j+1]-1]])
           binAFPAT2ETarr[c] = mean(AFPAT2ETarr[ri[ri[j]:ri[j+1]-1]])
           binACTENDTarr[c] = mean(ACTENDTarr[ri[ri[j]:ri[j+1]-1]])
           binAFPECTEarr[c] = mean(AFPECTEarr[ri[ri[j]:ri[j+1]-1]])
           binAFPEATEarr[c] = mean(AFPEATEarr[ri[ri[j]:ri[j+1]-1]])
           binASHTEMPEarr[c] = mean(ASHTEMPEarr[ri[ri[j]:ri[j+1]-1]])
           binATCTEMPEarr[c] = mean(ATCTEMPEarr[ri[ri[j]:ri[j+1]-1]])
           binACETEMPEarr[c] = mean(ACETEMPEarr[ri[ri[j]:ri[j+1]-1]])
           binAPDTEMPEarr[c] = mean(APDTEMPEarr[ri[ri[j]:ri[j+1]-1]])
           binACATMP1Earr[c] = mean(ACATMP1Earr[ri[ri[j]:ri[j+1]-1]])
           binACATMP2Earr[c] = mean(ACATMP2Earr[ri[ri[j]:ri[j+1]-1]])
           binACATMP3Earr[c] = mean(ACATMP3Earr[ri[ri[j]:ri[j+1]-1]])
           binACATMP4Earr[c] = mean(ACATMP4Earr[ri[ri[j]:ri[j+1]-1]])
           binACATMP5Earr[c] = mean(ACATMP5Earr[ri[ri[j]:ri[j+1]-1]])
           binACATMP6Earr[c] = mean(ACATMP6Earr[ri[ri[j]:ri[j+1]-1]])
           binACATMP7Earr[c] = mean(ACATMP7Earr[ri[ri[j]:ri[j+1]-1]])
           binACATMP8Earr[c] = mean(ACATMP8Earr[ri[ri[j]:ri[j+1]-1]])




           c = c + 1
                                ;print, 'testing', j, phasearr[ri[ri[j]:ri[j+1]-1]]
        endif
     endfor
     
     binxcen = binxcen[0:c-1]
     binycen = binycen[0:c-1]
     binnum = binnum[0:c-1]
     binbkgd = binbkgd[0:c-1]
     bintime = bintime[0:c-1]
     binflux = binflux[0:c-1]
     binback = binback[0:c-1]
     binbacksig = binbacksig[0:c-1]

     binVRSTUCCarr = binVRSTUCCarr[0:c-1]
     binAVRSTBEGarr = binAVRSTBEGarr[0:c-1]
     binAVDETCarr = binAVDETCarr[0:c-1]
     binAVDETBEGarr = binAVDETBEGarr[0:c-1]
     binAVGG1Carr = binAVGG1Carr[0:c-1]
     binAVGG1BEGarr = binAVGG1BEGarr[0:c-1]
     binAVDDUCCarr = binAVDDUCCarr[0:c-1]
     binAVDDUBEGarr = binAVDDUBEGarr[0:c-1]
     binAVGGCLCarr = binAVGGCLCarr[0:c-1]
     binAVGGCBEGarr = binAVGGCBEGarr[0:c-1]
     binAHTRIBEGarr = binAHTRIBEGarr[0:c-1]
     binAHTRVBEGarr = binAHTRVBEGarr[0:c-1]
     binAFPAT2Barr = binAFPAT2Barr[0:c-1]
     binAFPAT2BTarr = binAFPAT2BTarr[0:c-1]
     binAFPAT2Earr = binAFPAT2Earr[0:c-1]
     binAFPAT2ETarr = binAFPAT2ETarr[0:c-1]
     binACTENDTarr = binACTENDTarr[0:c-1]
     binAFPECTEarr = binAFPECTEarr[0:c-1]
     binAFPEATEarr = binAFPEATEarr[0:c-1]
     binASHTEMPEarr = binASHTEMPEarr[0:c-1]
     binATCTEMPEarr = binATCTEMPEarr[0:c-1]
     binACETEMPEarr = binACETEMPEarr[0:c-1]
     binAPDTEMPEarr = binAPDTEMPEarr[0:c-1]
     binACATMP1Earr = binACATMP1Earr[0:c-1]
     binACATMP2Earr = binACATMP2Earr[0:c-1]
     binACATMP3Earr = binACATMP3Earr[0:c-1]
     binACATMP4Earr = binACATMP4Earr[0:c-1]
     binACATMP5Earr = binACATMP5Earr[0:c-1]
     binACATMP6Earr = binACATMP6Earr[0:c-1]
     binACATMP7Earr = binACATMP7Earr[0:c-1]  
     

     save, /all, filename = '/Users/jkrick/irac_warm/pcrs_planets/55cnc/back_save_bin.sav'
     
     print, 'c', c
                                ;and plot
;  a = plot(binnum, binxcen, '1.', xtitle = 'number', ytitle = 'binned x cntrd', yrange = [14.4, 14.9])
;  a = plot(binnum, binycen, '1.', xtitle = 'number', ytitle = 'binned y cntrd')
;     a = plot(binnum, binbkgd, '1.', xtitle = 'number', ytitle = 'binned box_bkgdd')
;  a = plot(binnum, binflux, '1.', xtitle = 'number', ytitle = 'binned aper_flux')
     titlename = string(aorname(a)) +' '+ string(aorlabel) + string(exptime)
     b = plot((bintime - bintime(0))/60./60., binback, '1.', xtitle = 'Time (hrs)', ytitle = 'binned background (electrons)', title = titlename, yrange = [-0.5, 3])

     


;now for the correlation plots:
p = plot((bintime - bintime(0)) / 60./60., binbacksig, '1s', sym_size = 0.1, ytitle = 'background sigma',xtitle = 'time(hrs)')

p = plot(binback, binbacksig, '1s',sym_size = 0.1, xtitle = 'binned background (electrons)', ytitle = 'background sigma')

;p = plot(binbacksig,binVRSTUCCarr, '1s', xtitle = 'background sigma', ytitle = 'VRSTUCC')
;p = plot(binbacksig,binAVRSTBEGarr, '1s', xtitle = 'background sigma', ytitle = 'AVRSTBEG')
;p = plot(binbacksig,binAVDETCarr, '1s', xtitle = 'background sigma', ytitle = 'AVDETC')
;p = plot(binbacksig,binAVDETBEGarr, '1s', xtitle = 'background sigma', ytitle = 'AVDETBEG')
;p = plot(binbacksig,binAVGG1Carr, '1s', xtitle = 'background sigma', ytitle = 'AVGG1C')
;p = plot(binbacksig,binAVGG1BEGarr, '1s', xtitle = 'background sigma', ytitle = 'AVGG1BEG')
;p = plot(binbacksig,binAVDDUCCarr, '1s', xtitle = 'background sigma', ytitle = 'AVDDUCC')
;p = plot(binbacksig, binAVDDUBEGarr,'1s', xtitle = 'background sigma', ytitle = 'AVDDUBEG')
;p = plot(binbacksig, binAVGGCLCarr,'1s', xtitle = 'background sigma', ytitle = 'AVGGCLC')
;p = plot(binbacksig, binAVGGCBEGarr,'1s', xtitle = 'background sigma', ytitle = 'AVGGCBEG')
;p = plot(binbacksig, binAHTRIBEGarr,'1s', xtitle = 'background sigma', ytitle = 'AHTRIBEG')
;p = plot(binbacksig, binAHTRVBEGarr,'1s', xtitle = 'background sigma', ytitle = 'AHTRVBEG')
;p = plot(binbacksig, binAFPAT2Barr,'1s', xtitle = 'background sigma', ytitle = 'AFPAT2B')
;p = plot(binbacksig, binAFPAT2BTarr,'1s', xtitle = 'background sigma', ytitle = 'AFPAT2BT')
;p = plot(binbacksig, binAFPAT2Earr,'1s', xtitle = 'background sigma', ytitle = 'AFPAT2E')
;p = plot(binbacksig, binAFPAT2ETarr,'1s', xtitle = 'background sigma', ytitle = 'AFPAT2ET')
;p = plot(binbacksig, binACTENDTarr,'1s', xtitle = 'background sigma', ytitle = 'ACTENDT')
;p = plot(binbacksig, binAFPECTEarr,'1s', xtitle = 'background sigma', ytitle = 'AFPECTE')
;p = plot(binbacksig, binAFPEATEarr,'1s', xtitle = 'background sigma', ytitle = 'AFPEATE')
;p = plot(binbacksig, binASHTEMPEarr,'1s', xtitle = 'background sigma', ytitle = 'ASHTEMPE')
;p = plot(binbacksig, binATCTEMPEarr,'1s', xtitle = 'background sigma', ytitle = 'ATCTEMPE')
;p = plot(binbacksig, binACETEMPEarr,'1s', xtitle = 'background sigma', ytitle = 'ACETEMPE')
;p = plot(binbacksig, binAPDTEMPEarr,'1s', xtitle = 'background sigma', ytitle = 'APDTEMPE')
;p = plot(binbacksig, binACATMP1Earr,'1s', xtitle = 'background sigma', ytitle = 'ACATMP1E')
;p = plot(binbacksig, binACATMP2Earr,'1s', xtitle = 'background sigma', ytitle = 'ACATMP2E')
;p = plot(binbacksig, binACATMP3Earr,'1s', xtitle = 'background sigma', ytitle = 'ACATMP3E')
;p = plot(binbacksig, binACATMP4Earr,'1s', xtitle = 'background sigma', ytitle = 'ACATMP4E')
;p = plot(binbacksig, binACATMP5Earr,'1s', xtitle = 'background sigma', ytitle = 'ACATMP5E')
;p = plot(binbacksig, binACATMP6Earr,'1s', xtitle = 'background sigma', ytitle = 'ACATMP6E')
;p = plot(binbacksig, binACATMP7Earr,'1s', xtitle = 'background sigma', ytitle = 'ACATMP7E')
;

p = plot((bintime - bintime(0)) / 60./60.,binVRSTUCCarr, '1s', xtitle = 'time(hrs)', ytitle = 'VRSTUCC')
p = plot((bintime - bintime(0)) / 60./60.,binAVRSTBEGarr, '1s', xtitle = 'time(hrs)', ytitle = 'AVRSTBEG')
p = plot((bintime - bintime(0)) / 60./60.,binAVDETCarr, '1s', xtitle = 'time(hrs)', ytitle = 'AVDETC')
p = plot((bintime - bintime(0)) / 60./60.,binAVDETBEGarr, '1s', xtitle = 'time(hrs)', ytitle = 'AVDETBEG')
p = plot((bintime - bintime(0)) / 60./60.,binAVGG1Carr, '1s', xtitle = 'time(hrs)', ytitle = 'AVGG1C')
p = plot((bintime - bintime(0)) / 60./60.,binAVGG1BEGarr, '1s', xtitle = 'time(hrs)', ytitle = 'AVGG1BEG')
p = plot((bintime - bintime(0)) / 60./60.,binAVDDUCCarr, '1s', xtitle = 'time(hrs)', ytitle = 'AVDDUCC')
p = plot((bintime - bintime(0)) / 60./60., binAVDDUBEGarr,'1s', xtitle = 'time(hrs)', ytitle = 'AVDDUBEG')
p = plot((bintime - bintime(0)) / 60./60., binAVGGCLCarr,'1s', xtitle = 'time(hrs)', ytitle = 'AVGGCLC')
p = plot((bintime - bintime(0)) / 60./60., binAVGGCBEGarr,'1s', xtitle = 'time(hrs)', ytitle = 'AVGGCBEG')
p = plot((bintime - bintime(0)) / 60./60., binAHTRIBEGarr,'1s', xtitle = 'time(hrs)', ytitle = 'AHTRIBEG')
p = plot((bintime - bintime(0)) / 60./60., binAHTRVBEGarr,'1s', xtitle = 'time(hrs)', ytitle = 'AHTRVBEG')
p = plot((bintime - bintime(0)) / 60./60., binAFPAT2Barr,'1s', xtitle = 'time(hrs)', ytitle = 'AFPAT2B')
p = plot((bintime - bintime(0)) / 60./60., binAFPAT2BTarr,'1s', xtitle = 'time(hrs)', ytitle = 'AFPAT2BT')
p = plot((bintime - bintime(0)) / 60./60., binAFPAT2Earr,'1s', xtitle = 'time(hrs)', ytitle = 'AFPAT2E')
p = plot((bintime - bintime(0)) / 60./60., binAFPAT2ETarr,'1s', xtitle = 'time(hrs)', ytitle = 'AFPAT2ET')
p = plot((bintime - bintime(0)) / 60./60., binACTENDTarr,'1s', xtitle = 'time(hrs)', ytitle = 'ACTENDT')
p = plot((bintime - bintime(0)) / 60./60., binAFPECTEarr,'1s', xtitle = 'time(hrs)', ytitle = 'AFPECTE')
p = plot((bintime - bintime(0)) / 60./60., binAFPEATEarr,'1s', xtitle = 'time(hrs)', ytitle = 'AFPEATE')
p = plot((bintime - bintime(0)) / 60./60., binASHTEMPEarr,'1s', xtitle = 'time(hrs)', ytitle = 'ASHTEMPE')
p = plot((bintime - bintime(0)) / 60./60., binATCTEMPEarr,'1s', xtitle = 'time(hrs)', ytitle = 'ATCTEMPE')
p = plot((bintime - bintime(0)) / 60./60., binACETEMPEarr,'1s', xtitle = 'time(hrs)', ytitle = 'ACETEMPE')
p = plot((bintime - bintime(0)) / 60./60., binAPDTEMPEarr,'1s', xtitle = 'time(hrs)', ytitle = 'APDTEMPE')
p = plot((bintime - bintime(0)) / 60./60., binACATMP1Earr,'1s', xtitle = 'time(hrs)', ytitle = 'ACATMP1E')
p = plot((bintime - bintime(0)) / 60./60., binACATMP2Earr,'1s', xtitle = 'time(hrs)', ytitle = 'ACATMP2E')
p = plot((bintime - bintime(0)) / 60./60., binACATMP3Earr,'1s', xtitle = 'time(hrs)', ytitle = 'ACATMP3E')
p = plot((bintime - bintime(0)) / 60./60., binACATMP4Earr,'1s', xtitle = 'time(hrs)', ytitle = 'ACATMP4E')
p = plot((bintime - bintime(0)) / 60./60., binACATMP5Earr,'1s', xtitle = 'time(hrs)', ytitle = 'ACATMP5E')
p = plot((bintime - bintime(0)) / 60./60., binACATMP6Earr,'1s', xtitle = 'time(hrs)', ytitle = 'ACATMP6E')
p = plot((bintime - bintime(0)) / 60./60., binACATMP7Earr,'1s', xtitle = 'time(hrs)', ytitle = 'ACATMP7E')


endfor
end
 
