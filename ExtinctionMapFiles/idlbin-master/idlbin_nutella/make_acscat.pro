pro make_acscat

close,/all
openw, outlun, '/Users/jkrick/hst/acscat_merged.txt', /get_lun
filelist = ['/Users/jkrick/hst/a_catalog.fits', '/Users/jkrick/hst/b_catalog.fits', '/Users/jkrick/hst/c_catalog.fits', '/Users/jkrick/hst/d_catalog.fits', '/Users/jkrick/hst/e_catalog.fits', '/Users/jkrick/hst/f_catalog.fits', '/Users/jkrick/hst/g_catalog.fits']
imagelist = ['/Users/jkrick/hst/a_drz.fits', '/Users/jkrick/hst/b_drz.fits', '/Users/jkrick/hst/c_drz.fits', '/Users/jkrick/hst/d_drz.fits', '/Users/jkrick/hst/e_drz.fits', '/Users/jkrick/hst/f_drz.fits', '/Users/jkrick/hst/g_drz.fits']
count = 0l
for i = 0, n_elements(filelist) - 1 do begin
   print, "working on ", filelist(i)
   ftab_ext,filelist(i), 'Concentration,Gini,mag_auto,SurfaceBrightness,x_image,y_image,class_star,a_image,b_image', c,g,mag,mu,x,y,s,a,b
   ftab_ext, filelist(i), 'flags,x_world,y_world,isoarea_image,theta_world,signaltonoiseratio,fluxcounts,sky,skyrms',flag,ra,dec,isoarea,theta,snr,flux, sky,skyrms
   ftab_ext, filelist(i), 'centralsurfacebrightness,rotationalasymmetry,halflightradiusinpixels,m20,class_star,magerr_auto',cmu,asym,hlr,m20,class,magerr

   print, n_elements(c)
   fits_read, imagelist(i), data, header
   ;also need the header for the extent of the image
   naxis1=sxpar(header,'NAXIS1')-2
   naxis2=sxpar(header,'NAXIS2')-2

   for j = 0l, n_elements(c) -1 do begin

      if mag(j) lt 20. and x(j) gt 100. and y(j) gt 100. and  x(j) lt naxis1 - 100. and y(j) lt naxis2 -100. then begin

         if data[x(j)-100.,y(j)] ne 0. and data[x(j)+100.,y(j)] ne 0. and data[x(j),y(j)-100.] ne 0. and data[x(j),y(j)+100.] ne 0. then begin
            printf, outlun,  format='(F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.5,F10.5,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)',c(j),g(j),mag(j),mu(j),x(j),y(j),s(j),a(j),b(j),flag(j),ra(j),dec(j),isoarea(j),theta(j),snr(j),flux(j),sky(j),skyrms(j),cmu(j),asym(j),hlr(j),m20(j),class(j),magerr(j)
            count = count + 1
         endif

      endif
      if mag(j) lt 23. and mag(j) ge 20. and x(j) gt 100. and y(j) gt 100. and  x(j) lt naxis1 - 100. and y(j) lt naxis2 -100. then begin

         if data[x(j)-70.,y(j)] ne 0. and data[x(j)+70.,y(j)] ne 0. and data[x(j),y(j)-70.] ne 0. and data[x(j),y(j)+70.] ne 0. then begin
            printf, outlun, format='(F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.5,F10.5,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)',c(j),g(j),mag(j),mu(j),x(j),y(j),s(j),a(j),b(j),flag(j),ra(j),dec(j),isoarea(j),theta(j),snr(j),flux(j),sky(j),skyrms(j),cmu(j),asym(j),hlr(j),m20(j),class(j),magerr(j)
            count = count + 1

         endif

      endif
      if mag(j) ge 23. and x(j) gt 100. and y(j) gt 100. and  x(j) lt naxis1 - 100. and y(j) lt naxis2 -100. then begin

         if data[x(j)-50.,y(j)] ne 0. and data[x(j)+50.,y(j)] ne 0. and data[x(j),y(j)-50.] ne 0. and data[x(j),y(j)+50.] ne 0. then begin
            printf, outlun, format='(F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.5,F10.5,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2)',c(j),g(j),mag(j),mu(j),x(j),y(j),s(j),a(j),b(j),flag(j),ra(j),dec(j),isoarea(j),theta(j),snr(j),flux(j),sky(j),skyrms(j),cmu(j),asym(j),hlr(j),m20(j),class(j),magerr(j)
            count = count + 1

         endif

      endif

   endfor


endfor
print, "count", count
close, outlun
free_lun, outlun
end
;ftab_ext, filelist(i), 'NUMBER,MAG_AUTO,MAGERR_AUTO,BACKGROUND,THRESHOLD,MU_THRESHOLD,FLUX_MAX,MU_MAX AREA_IMAGE,ISOAREA_WORLD,XMIN_IMAGE,YMIN_IMAGE,XMAX_IMAGE,YMAX_IMAGE,X_IMAGE,Y_IMAGE,X_WORLD,Y_WORLD,PEAK_IMAGE,PEAK_IMAGE,XPEAK_WORLD,YPEAK_WORLD,ALPHA_J2000,DELTA_J2000,X2_IMAGE, Y2_IMAGE,XY_IMAGE,X2_WORLD,Y2_WORLD, XY_WORLD,A_IMAGE,B_IMAGE,A_WORLD,B_WORLD,THETA_IMAGE,THETA_WORLD,FLAGS,CLASS_STAR,UPhi,ZeroPoint,Area,SegmentNumber,SignalToNoiseRatio,IsophotalMagnitude,Sky,MaxFluxCounts,CR58,SurfaceBrightness,CR16,CR17,R18,MinimumFluxY,MinimumFluxX,ApertureMagnitude,CR18,CentroidX0,CentralSurfaceBrightness,R16,AxisRatio,Phi,Gini,MinimumAsymmetryY0,ErrorFlags,I1,UA,BadValueIndicator,I5,RawAsymmetry,I7,I6,NX,NY,R58,PixelScale,CR56,CR68,R78,X0,CR57,R56,R57,SkyRMS,CR67,I8,M20,RotationalAsymmetry,MaximumFluxY,MaximumFluxX,CorrectionForSkyAsymmetry,HalfLightRadiusInPixels,R15,CR78,R17,CentroidY0,CentralPixelCounts,MinimumAsymmetryX0,CR15,MinorAxis,Concentration,FluxCounts,R67,MajorAxis,ApertureRadiusInPixels,Y0,QuasiPetrosianFraction,R68,MinFluxCounts,UseCentroid,UB,KEY', v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14, v15,v16,v17,v18,v19,v20,v21,v22,v23,v24,v25,v26,v27,v28,v29,v30,v31,v32,v33,v34,v35,v36,v37,v38,v39,v40,v41,v42,v43,v44,v45,v46,v47,v48,v49,v50,v51,v52,v53,v54,v55,v56,v57,v58,v59,v60,v61,v62,v63,v64,v65,v66,v67,v68,v69,v70,v71,v72,v73,v74,v75,v76,v77,v78,v79,v80,v81,v82,v83,v84,v85,v86,v87,v88,v89,v90,v91,v92,v93,v94,v95,v96,v97,v98,v99,v100,v101,v102,v103,v104,v105,v106,v107,v108,v109,v110
