pro plot_hd158460_snap
   colorarr = [  'SILVER'   , 'SKY_BLUE'  , 'SLATE_BLUE' , 'SLATE_GREY' , ' red',  'PLUM' ,  'SPRING_GREEN'  ,  'STEEL_BLUE'   ,  'green', 'TEAL' ,  'THISTLE' ,  'TOMATO'  ,  'TURQUOISE' ,  'VIOLET'  , 'YELLOW'  ,  'YELLOW_GREEN', 'PAPAYA_WHIP'  ,' PEACH_PUFF','   PERU',   ' POWDER_BLUE' ]
 
;!P.MULTI = [0, 1, 2]

;January 2012 snapshot take two:
aorname = ['0045184256','0045184512','0045184768','0045185024','0045185280','0045185536','0045185792','0045186048','0045186304','0045186560','0045186816']
restore, '/Users/jkrick/irac_warm/pmap/hd158460_v2.sav'


;  restore, '/Users/jkrick/irac_warm/snapshots/allsnap_corr.sav'

; aorname = ['r44497408','r44497152','r44497664','r44497920','r44498176','r44498432','r44498688','r44498944','r44499200','r44499456','r44499712', 'r44499968','r44500224','r44500480','r44500736', 'r44500992', 'r44501248', 'r44501504' ]
  ;remove those with no pmap coverage
;  aorname = ['r44497152','r44497664','r44497920','r44498176','r44498432','r44498688','r44498944','r44499200','r44499456','r44500224','r44500480','r44500736', 'r44500992', 'r44501248', 'r44501504' ]

 ;read in the coverage map for the pixel phase effect
  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_500x500_0043_111129/occu_ch2_500x500_0043_111129.fits', occdata, occheader

  ;do some normalization
 ; for a = 0, n_elements(aorname) - 1 do begin
 ;    snapshots[a].flux = snapshots[a].flux / snapshots[a].flux(0)
 ;    snapshots[a].fluxerr = snapshots[a].fluxerr / snapshots[a].flux(0)
 ;    snapshots[a].corrflux = snapshots[a].corrflux / snapshots[a].corrflux(0)
 ;    snapshots[a].corrfluxerr = snapshots[a].corrfluxerr / snapshots[a].corrflux(0)
 ; endfor

z = pp_multiplot(multi_layout=[n_elements(aorname), 1], global_xtitle='X (pixels)',global_ytitle='Y (pixels)')

  ;have to plot each row one at a time
  for a = 0,  n_elements(aorname) - 1 do begin
     xy=z.plot(snapshots[a].xcen, snapshots[a].ycen, '6r1.',color = 'black', thick = 2, xminor = 0, xtickinterval = 0.5)
     xy.yrange = [14.5, 15.5]
     xy.xrange = [14.5, 15.5]
  endfor
 
 o = pp_multiplot(multi_layout=[n_elements(aorname),1],global_xtitle='Time (hours)',global_ytitle='Flux (mJy)')

 for a = 0,  n_elements(aorname) - 1 do begin
    st = o.errorplot(snapshots[a].timearr, snapshots[a].flux,snapshots[a].fluxerr,  color = colorarr[a], xminor = 0, xtickinterval = 0.5, multi_index = a) ;
    st.yrange = [0.97,1.1]
    st.xrange = [0,0.5]
    sto = o.errorplot(snapshots[a].timearr, snapshots[a].corrflux- 0.04, snapshots[a].corrfluxerr, multi_index = a)
    sto.yrange = [0.97,1.1]
    sto.xrange = [0,0.5]
                                ; ;  st.xminor = 0
 endfor
  
;  plotname.xrange = [14.5, 15.5]
;  print, 'xranges', m.xranges

;  m.sync_axes, 9
;  print, 'xranges after sync', m.xranges

;another plot
  for a = 0, n_elements(aorname) - 1 do begin
     if a eq 0 then begin
        an = plot(snapshots[a].xcen, snapshots[a].ycen, '6r1.', thick = 2, xrange = [14.5, 15.5], yrange = [14.5, 15.5], xtitle = 'X (pixel)', ytitle = 'Y (pixel)', color = colorarr[a], aspect_ratio = 1)
     endif

     if a gt 0 then begin
        an = plot(snapshots[a].xcen, snapshots[a].ycen, '6r1.', thick = 2, color = colorarr[a],/overplot)
     endif
  endfor

  xsweet = 15.120
  ysweet = 15.085  
  box_x = [xsweet-0.1, xsweet-0.1, xsweet + 0.1, xsweet + 0.1, xsweet -0.1]
  box_y = [ysweet-0.1, ysweet +0.1, ysweet +0.1, ysweet - 0.1,ysweet -0.1]
  line4 = polyline(box_x, box_y, thick = 2, color = !color.black,/data)
 
end
