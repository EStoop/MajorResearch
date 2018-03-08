pro bc03_color

!P.multi = [0,1,1]
;readlargecol, '/Users/jkrick/bin/bc03/models/Padova1994/salpeter/lr_m62_salp_ssp.color_F117_F124', z, ltt, tz, dm, c_rf, c_ne, c_ev, c_ab_rf, c_ab_ne, c_ab_ev, ekcor, kcor,format="A"

;readlargecol, '/Users/jkrick/bin/bc03/models/Padova1994/salpeter/lr_m62_salp_ssp.color_F123_F124', z1, ltt1, tz1, dm1, c_rf1, c_ne1, c_ev1, c_ab_rf1, c_ab_ne1, c_ab_ev1, ekcor1, kcor1,format="A"

;readlargecol, '/Users/jkrick/bin/bc03/models/Padova1994/salpeter/lr_m62_salp_ssp.color_F116_F123', z2, ltt2, tz2, dm2, c_rf2, c_ne2, c_ev2, c_ab_rf2, c_ab_ne2, c_ab_ev2, ekcor2, kcor2,format="A"

;readlargecol, '/Users/jkrick/bin/bc03/models/Padova1994/salpeter/lr_m62_salp_ssp.color_F116_F117', z3, ltt3, tz3, dm3, c_rf3, c_ne3, c_ev3, c_ab_rf3, c_ab_ne3, c_ab_ev3, ekcor3, kcor3,format="A"

ps_open, filename='/Users/jkrick/virgo/colorpredict.ps',/portrait,/square,/color

redcolor = FSC_COLOR("Red", !D.Table_Size-2);
orangecolor = FSC_COLOR("Orange", !D.Table_Size-3);
greencolor = FSC_COLOR("Green", !D.Table_Size-4);
bluecolor = FSC_COLOR("Blue", !D.Table_Size-5)

readlargecol, '/Users/jkrick/bc03/models/Padova1994/salpeter/bc2003_lr_m62_salp_ssp.color_F124_F129', z, ltt, tz, dm, c_rf, c_ne, c_ev, c_ab_rf, c_ab_ne, c_ab_ev, ekcor, kcor,format="A"

plot, tz, c_ab_ne, xtitle = 'age of galaxy', ytitle = 'ch1 - ch2 color', charthick = 3, xthick = 3, ythick = 3, thick = 3

readlargecol, '/Users/jkrick/bc03/models/Padova1994/salpeter/bc2003_lr_m52_salp_ssp.color_F124_F129', z, ltt, tz, dm, c_rf, c_ne, c_ev, c_ab_rf, c_ab_ne, c_ab_ev, ekcor, kcor,format="A"

oplot, tz, c_ab_ne, color = redcolor

readlargecol, '/Users/jkrick/bc03/models/Padova1994/salpeter/bc2003_lr_m42_salp_ssp.color_F124_F129', z, ltt, tz, dm, c_rf, c_ne, c_ev, c_ab_rf, c_ab_ne, c_ab_ev, ekcor, kcor,format="A"

oplot, tz, c_ab_ne, color = orangecolor


readlargecol, '/Users/jkrick/bc03/models/Padova1994/salpeter/bc2003_lr_m32_salp_ssp.color_F124_F129', z, ltt, tz, dm, c_rf, c_ne, c_ev, c_ab_rf, c_ab_ne, c_ab_ev, ekcor, kcor,format="A"

oplot, tz, c_ab_ne, color = greencolor


readlargecol, '/Users/jkrick/bc03/models/Padova1994/salpeter/bc2003_lr_m22_salp_ssp.color_F124_F129', z, ltt, tz, dm, c_rf, c_ne, c_ev, c_ab_rf, c_ab_ne, c_ab_ev, ekcor, kcor,format="A"

oplot, tz, c_ab_ne, color = bluecolor



readlargecol, '/Users/jkrick/bc03/models/Padova1994/salpeter/bc2003_hr_m52_salp_ssp.color_F124_F129', z, ltt, tz, dm, c_rf, c_ne, c_ev, c_ab_rf, c_ab_ne, c_ab_ev, ekcor, kcor,format="A"

oplot, tz, c_ab_ne, linestyle = 2, color = redcolor


readlargecol, '/Users/jkrick/bc03/models/Padova1994/salpeter/bc2003_hr_m62_salp_ssp.color_F124_F129', z, ltt, tz, dm, c_rf, c_ne, c_ev, c_ab_rf, c_ab_ne, c_ab_ev, ekcor, kcor,format="A"

oplot, tz, c_ab_ne, linestyle = 2




;oplot, z, c_ab_ne1, color = redcolor, thick = 3
;oplot, z , c_ab_ne3, color = bluecolor, thick = 3
;xyouts, 0.1, 7.2, 'ACS F814 - ch1', color  = redcolor, charthick = 3;xyouts, 0.1, 7.5, "palomar r' - ch1", charthick = 3
;xyouts, 0.1, 7.0, "palomar g' - r' ", charthick = 3, color = bluecolor
ps_close, /noprint,/noid

end
