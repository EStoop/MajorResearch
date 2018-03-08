pro read_pmap

  restore, '/Users/jkrick/irac_warm/pmap/pmap_ch2_0p1s_x4_rmulti_s3_7_photometry.sav'

  ;;clean the data a bit; only keep sweet spot data
  ;;would prefer these were sigma clipped

  meanclip, ch2dat.xpos, meanxpos, sigmaxpos
  print,'meanx, sigmax',  meanxpos - 3*sigmaxpos, meanxpos + 3*sigmaxpos
  meanclip, ch2dat.ypos, meanypos, sigmaypos
  print,'meany, sigmay',  meanypos - 3*sigmaypos, meanypos + 3*sigmaypos
  meanclip, ch2dat.flux, meanflux, sigmaflux
  print,'meanflux, sigma',  meanflux - 1*sigmaflux, meanflux + 1*sigmaflux
  
  meanclip, ch2dat.dn_peak, meandn_peak, sigmadn_peak
  print,'meandn_peak, sigma',  meandn_peak - 3*sigmadn_peak, meandn_peak + 3*sigmadn_peak


  meanclip, ch2dat.sigma_flux, meansigma_flux, sigmasigma_flux
  print,'meansigma_flux, sigma',  meansigma_flux - 3*sigmasigma_flux, meansigma_flux + 3*sigmasigma_flux


    meanclip, ch2dat.bg_flux, meanbg_flux, sigmabg_flux
  print,'meanbg_flux, sigma',  meanbg_flux - 3*sigmabg_flux, meanbg_flux + 3*sigmabg_flux


    meanclip, ch2dat.sigma_bg_flux, meansigma_bg_flux, sigmasigma_bg_flux
  print,'meansigma_bg_flux, sigma',  meansigma_bg_flux - 3*sigmasigma_bg_flux, meansigma_bg_flux + 3*sigmasigma_bg_flux

 ;; ss = where(ch2dat.xpos ge 14.5 and ch2dat.xpos le 15.5 and $
 ;;            ch2dat.ypos ge 14.5 and ch2dat.ypos le 15.5 and $
 ;;            ch2dat.flux[*,5] ge 0.39 and ch2dat.flux[*,5] le 0.41 $
 ;;            and ch2dat.dn_peak lt 3500 $
 ;;            and ch2dat.sigma_flux lt 0.005 $
 ;;            and ch2dat.bg_flux le 0.0004 $
 ;;            and ch2dat.sigma_bg_flux le 0.0004)
  ss = where(ch2dat.xpos gt 14.5 and ch2dat.xpos lt 15.5 and $
             ch2dat.ypos gt 14.5 and ch2dat.ypos lt 15.5 and $
             ch2dat.flux[*,5] ge (meanflux - 1*sigmaflux) and ch2dat.flux[*,5] le 0.41 $
             and ch2dat.dn_peak lt 3500 $
             and ch2dat.sigma_flux gt (meansigma_flux - 3*sigmasigma_flux) and ch2dat.sigma_flux lt (meansigma_flux + 3*sigmasigma_flux)$
            ;; and ch2dat.bg_flux gt (meanbg_flux - 3*sigmabg_flux) and ch2dat.bg_flux lt (meanbg_flux + 3*sigmabg_flux) $
             and ch2dat.sigma_bg_flux gt (meansigma_bg_flux - 3*sigmasigma_bg_flux) and ch2dat.sigma_bg_flux lt (meansigma_bg_flux + 3*sigmasigma_bg_flux))
  
  xpos = ch2dat.xpos(ss)
  xerr = ch2dat.xerr(ss)
  ypos = ch2dat.ypos(ss)
  yerr = ch2dat.yerr(ss)
  xycov = ch2dat.xycov(ss)
  flux = ch2dat.flux[*,5]
  flux = flux(ss)
  sigma_flux = ch2dat.sigma_flux[*,5]
  sigma_flux = sigma_flux(ss)
  noise_pix = ch2dat.noise_pix(ss)
  xfwhm = ch2dat.xfwhm(ss)
  yfwhm = ch2dat.yfwhm(ss)
  dn_peak = ch2dat.dn_peak(ss)
  bmjd = ch2dat.bmjd(ss)
  t_cernox = ch2dat.t_cernox(ss)
  bg_flux = ch2dat.bg_flux(ss)
  sigma_bg_flux = ch2dat.sigma_bg_flux(ss)
  pix2_2 = reform(ch2dat.pix7x7[2,2,*])
  pix2_2 = pix2_2(ss)
  pix2_3 = reform(ch2dat.pix7x7[2,3,*])
  pix2_3 = pix2_3(ss)
  pix2_4 = reform(ch2dat.pix7x7[2,4,*])
  pix2_4 = pix2_4(ss)

  pix3_2 = reform(ch2dat.pix7x7[3,2,*])
  pix3_2 = pix3_2(ss)
  pix3_3 = reform(ch2dat.pix7x7[3,3,*])
  pix3_3 = pix3_3(ss)
  pix3_4 = reform(ch2dat.pix7x7[3,4,*])
  pix3_4 = pix3_4(ss)

  pix4_2 = reform(ch2dat.pix7x7[4,2,*])
  pix4_2 = pix4_2(ss)
  pix4_3 = reform(ch2dat.pix7x7[4,3,*])
  pix4_3 = pix4_3(ss)
  pix4_4 = reform(ch2dat.pix7x7[4,4,*])
  pix4_4 = pix4_4(ss)


 ; plothist, sigma_bg_flux, xhist, yhist, bin = 0.00001, /noprint,/noplot,/NAN
 ; a = barplot(xhist, yhist, title = 'sigma_bg_flux')

  

  
  ;;make a 2D array to hold variables instead of struct

  all_arr = [[xpos],[xerr], [ypos],[yerr],[xycov], [flux], [sigma_flux], [noise_pix], [xfwhm], [yfwhm], [dn_peak], [bmjd], [t_cernox], [bg_flux], [sigma_bg_flux], [pix2_2] , [pix2_3], [pix2_4], [pix3_2] , [pix3_3], [pix3_4], [pix4_2] , [pix4_3], [pix4_4]]
  help, all_arr

  trans_arr = transpose(all_arr)
  help, trans_arr
  header = ['xpos','xerr','ypos','yerr','xycov','flux', 'fluxerr','np', 'xfwhm','yfwhm','dn_peak', 'bmjd', 't_cernox', 'bg_flux', 'sigma_bg_flux','pix1','pix2', 'pix3','pix4','pix5','pix6','pix7','pix8','pix9']


  
  

  
  ;;use write_csv
  
  write_csv, '//Users/jkrick/irac_warm/pmap/pmap_ch2_0p1s_x4_rmulti_s3_7.csv', trans_arr, header = header
  
  print, n_elements(xpos), max(flux)
end
; [pix2_2] , [pix2_3], [pix2_4], [pix3_2] , [pix3_3], [pix3_4], [pix4_2] , [pix4_3], [pix4_4]
