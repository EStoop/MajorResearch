pro plot_back_explosion

  restore,  '/Users/jkrick/irac_warm/pcrs_planets/55cnc/back_save.sav' ;r43981568
  ;  restore,  '/Users/jkrick/irac_warm/pcrs_planets/55cnc/back_save.sav' ;r43981568

;color code those things with high background and those with low background.

high = where(binback gt 1.4,nhigh)
low = where(binback le 1.4,nlow)
lowlow = where(binback lt 0.6, nlowlow)
highhigh = where(binback gt 2.1, nhighhigh)
print, 'nlow, nhigh', nlow, nhigh

;a= plot((bintime(high) - bintime(0))/60./60., binback(high), '1sb', xtitle = 'Time (hrs)', ytitle = 'binned background (electrons)', title ='r43981568' )

;a = plot((bintime(low) - bintime(0))/60./60., binback(low), '1sr',/overplot)

;b = plot((bintime(high) - bintime(0))/60./60, binbacksig(high), '1sb', xtitle = 'Time (hrs)', ytitle = 'background sigma (electrons)', title ='r43981568' )
;b = plot((bintime(low) - bintime(0))/60./60., binbacksig(low), '1sr',/overplot)

;c = plot(binback(high), binbacksig(high),  '1sb',  ytitle = 'background sigma (electrons)', xtitle = 'binned background (electrons)', title ='r43981568' )
;c = plot(binback(low), binbacksig(low),  '1sr',  /overplot)

a= plot((bintime(high) - bintime(0))/60./60., binback(high), xtitle = 'Time (hrs)', ytitle = 'binned background (electrons)', title ='r43981568', xrange = [0.5, 3.0] )
a= plot((bintime(high) - bintime(0))/60./60., binback(high), '1*', /overplot)
;--------------------------
;now try correlations with header telemetry keywords


;sp = plot(binbacksig,binVRSTUCCarr, '1s', xtitle = 'background sigma', ytitle = 'VRSTUCC')
;sp = plot((bintime - bintime(0)) / 60./60.,binVRSTUCCarr, '1s', xtitle = 'time(hrs)', ytitle = 'VRSTUCC')
;p = plot(binbacksig,binAVRSTBEGarr, '1s', xtitle = 'background sigma', ytitle = 'AVRSTBEG')
;p = plot((bintime - bintime(0)) / 60./60.,binAVRSTBEGarr, '1s', xtitle = 'time(hrs)', ytitle = 'AVRSTBEG')
;sp = plot(binbacksig,binAVDETCarr, '1s', xtitle = 'background sigma', ytitle = 'AVDETC')
;sp = plot((bintime - bintime(0)) / 60./60.,binAVDETCarr, '1s', xtitle = 'time(hrs)', ytitle = 'AVDETC')
;p = plot(binbacksig,binAVDETBEGarr, '1s', xtitle = 'background sigma', ytitle = 'AVDETBEG')
;p = plot((bintime - bintime(0)) / 60./60.,binAVDETBEGarr, '1s', xtitle = 'time(hrs)', ytitle = 'AVDETBEG')
;sp = plot(binbacksig,binAVGG1Carr, '1s', xtitle = 'background sigma', ytitle = 'AVGG1C')
;sp = plot((bintime - bintime(0)) / 60./60.,binAVGG1Carr, '1s', xtitle = 'time(hrs)', ytitle = 'AVGG1C')
;p = plot(binbacksig,binAVGG1BEGarr, '1s', xtitle = 'background sigma', ytitle = 'AVGG1BEG')
;p = plot((bintime - bintime(0)) / 60./60.,binAVGG1BEGarr, '1s', xtitle = 'time(hrs)', ytitle = 'AVGG1BEG')
;sp = plot(binbacksig,binAVDDUCCarr, '1s', xtitle = 'background sigma', ytitle = 'AVDDUCC')
;sp = plot((bintime - bintime(0)) / 60./60.,binAVDDUCCarr, '1s', xtitle = 'time(hrs)', ytitle = 'AVDDUCC')
;p = plot(binbacksig, binAVDDUBEGarr,'1s', xtitle = 'background sigma', ytitle = 'AVDDUBEG')
;p = plot((bintime - bintime(0)) / 60./60., binAVDDUBEGarr,'1s', xtitle = 'time(hrs)', ytitle = 'AVDDUBEG')
;sp = plot(binbacksig, binAVGGCLCarr,'1s', xtitle = 'background sigma', ytitle = 'AVGGCLC')
;sp = plot((bintime - bintime(0)) / 60./60., binAVGGCLCarr,'1s', xtitle = 'time(hrs)', ytitle = 'AVGGCLC')
;sp = plot(binbacksig, binAVGGCBEGarr,'1s', xtitle = 'background sigma', ytitle = 'AVGGCBEG')
;sp = plot((bintime - bintime(0)) / 60./60., binAVGGCBEGarr,'1s', xtitle = 'time(hrs)', ytitle = 'AVGGCBEG')
;p = plot(binbacksig, binAHTRIBEGarr,'1s', xtitle = 'background sigma', ytitle = 'AHTRIBEG')
;p = plot((bintime - bintime(0)) / 60./60., binAHTRIBEGarr,'1s', xtitle = 'time(hrs)', ytitle = 'AHTRIBEG')
;p = plot(binbacksig, binAHTRVBEGarr,'1s', xtitle = 'background sigma', ytitle = 'AHTRVBEG')
;p = plot((bintime - bintime(0)) / 60./60., binAHTRVBEGarr,'1s', xtitle = 'time(hrs)', ytitle = 'AHTRVBEG')
;p = plot(binbacksig, binAFPAT2Barr,'1s', xtitle = 'background sigma', ytitle = 'AFPAT2B')
;p = plot((bintime - bintime(0)) / 60./60., binAFPAT2Barr,'1s', xtitle = 'time(hrs)', ytitle = 'AFPAT2B')
;sp = plot(binbacksig, binAFPAT2BTarr,'1s', xtitle = 'background sigma', ytitle = 'AFPAT2BT')
;sp = plot((bintime - bintime(0)) / 60./60., binAFPAT2BTarr,'1s', xtitle = 'time(hrs)', ytitle = 'AFPAT2BT')
;p = plot(binbacksig, binAFPAT2Earr,'1s', xtitle = 'background sigma', ytitle = 'AFPAT2E')
;p = plot((bintime - bintime(0)) / 60./60., binAFPAT2Earr,'1s', xtitle = 'time(hrs)', ytitle = 'AFPAT2E')
;sp = plot(binbacksig, binAFPAT2ETarr,'1s', xtitle = 'background sigma', ytitle = 'AFPAT2ET')
;sp = plot((bintime - bintime(0)) / 60./60., binAFPAT2ETarr,'1s', xtitle = 'time(hrs)', ytitle = 'AFPAT2ET')
;p = plot(binbacksig, binACTENDTarr,'1s', xtitle = 'background sigma', ytitle = 'ACTENDT')
;p = plot((bintime - bintime(0)) / 60./60., binACTENDTarr,'1s', xtitle = 'time(hrs)', ytitle = 'ACTENDT')
;p = plot(binbacksig, binAFPECTEarr,'1s', xtitle = 'background sigma', ytitle = 'AFPECTE')
;p = plot((bintime - bintime(0)) / 60./60., binAFPECTEarr,'1s', xtitle = 'time(hrs)', ytitle = 'AFPECTE')
;p = plot(binbacksig, binAFPEATEarr,'1s', xtitle = 'background sigma', ytitle = 'AFPEATE')
;p = plot((bintime - bintime(0)) / 60./60., binAFPEATEarr,'1s', xtitle = 'time(hrs)', ytitle = 'AFPEATE')
;p = plot(binbacksig, binASHTEMPEarr,'1s', xtitle = 'background sigma', ytitle = 'ASHTEMPE')
;p = plot((bintime - bintime(0)) / 60./60., binASHTEMPEarr,'1s', xtitle = 'time(hrs)', ytitle = 'ASHTEMPE')
;p = plot(binbacksig, binATCTEMPEarr,'1s', xtitle = 'background sigma', ytitle = 'ATCTEMPE')
;p = plot((bintime - bintime(0)) / 60./60., binATCTEMPEarr,'1s', xtitle = 'time(hrs)', ytitle = 'ATCTEMPE')
;p = plot(binbacksig, binACETEMPEarr,'1s', xtitle = 'background sigma', ytitle = 'ACETEMPE')
;p = plot((bintime - bintime(0)) / 60./60., binACETEMPEarr,'1s', xtitle = 'time(hrs)', ytitle = 'ACETEMPE')
;p = plot(binbacksig, binAPDTEMPEarr,'1s', xtitle = 'background sigma', ytitle = 'APDTEMPE')
;p = plot((bintime - bintime(0)) / 60./60., binAPDTEMPEarr,'1s', xtitle = 'time(hrs)', ytitle = 'APDTEMPE')
;p = plot(binbacksig, binACATMP1Earr,'1s', xtitle = 'background sigma', ytitle = 'ACATMP1E')
;p = plot((bintime - bintime(0)) / 60./60., binACATMP1Earr,'1s', xtitle = 'time(hrs)', ytitle = 'ACATMP1E')
;p = plot(binbacksig, binACATMP2Earr,'1s', xtitle = 'background sigma', ytitle = 'ACATMP2E')
;p = plot((bintime - bintime(0)) / 60./60., binACATMP2Earr,'1s', xtitle = 'time(hrs)', ytitle = 'ACATMP2E')
;p = plot(binbacksig, binACATMP3Earr,'1s', xtitle = 'background sigma', ytitle = 'ACATMP3E')
;p = plot((bintime - bintime(0)) / 60./60., binACATMP3Earr,'1s', xtitle = 'time(hrs)', ytitle = 'ACATMP3E')
;p = plot(binbacksig, binACATMP4Earr,'1s', xtitle = 'background sigma', ytitle = 'ACATMP4E')
;p = plot((bintime - bintime(0)) / 60./60., binACATMP4Earr,'1s', xtitle = 'time(hrs)', ytitle = 'ACATMP4E')
;p = plot(binbacksig, binACATMP5Earr,'1s', xtitle = 'background sigma', ytitle = 'ACATMP5E')
;p = plot((bintime - bintime(0)) / 60./60., binACATMP5Earr,'1s', xtitle = 'time(hrs)', ytitle = 'ACATMP5E')
;p = plot(binbacksig, binACATMP6Earr,'1s', xtitle = 'background sigma', ytitle = 'ACATMP6E')
;p = plot((bintime - bintime(0)) / 60./60., binACATMP6Earr,'1s', xtitle = 'time(hrs)', ytitle = 'ACATMP6E')
;p = plot(binbacksig, binACATMP7Earr,'1s', xtitle = 'background sigma', ytitle = 'ACATMP7E')
;p = plot((bintime - bintime(0)) / 60./60., binACATMP7Earr,'1s', xtitle = 'time(hrs)', ytitle = 'ACATMP7E')

end
