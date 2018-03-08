pro polletta
close, /all
redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
greencolor = FSC_COLOR("Green", !D.Table_Size-4)
yellowcolor = FSC_COLOR("Yellow", !D.Table_Size-5)
cyancolor = FSC_COLOR("cyan", !D.Table_Size-6)
orangecolor = FSC_COLOR("orange", !D.Table_Size-7)
purplecolor = FSC_COLOR("purple", !D.Table_Size-8)

readlargecol,'/Users/jkrick/ZPHOT/templates/XMM_LSS/Sa_flux_z.dat',z, Ul,U,B,V,R,I,up,gp,rp,ip,zp,J,Ks,irac1,irac2,irac3,irac4,mips24,Ulmag,Umag,Bmag,Vmag,Rmag,Imag,upmag,gpmag,rpmag,ipmag,zpmag,Jmag,Ksmag,irac1mag,irac2mag,irac3mag,irac4mag,mips24mag,format="A"

count = [16,23,27,29,31]

plot, alog10(irac3(count)/irac1(count)), alog10(irac4(count)/irac2(count)), psym = 4, xrange=[-1,1], yrange=[-1.5,1]
oplot, alog10(irac3(count)/irac1(count)), alog10(irac4(count)/irac2(count)), color = redcolor
xyouts,  alog10(irac3(16)/irac1(16)), alog10(irac4(16)/irac2(16)), ".", color = redcolor

readlargecol,'/Users/jkrick/ZPHOT/templates/XMM_LSS/Sd_flux_z.dat',z, Ul,U,B,V,R,I,up,gp,rp,ip,zp,J,Ks,irac1,irac2,irac3,irac4,mips24,Ulmag,Umag,Bmag,Vmag,Rmag,Imag,upmag,gpmag,rpmag,ipmag,zpmag,Jmag,Ksmag,irac1mag,irac2mag,irac3mag,irac4mag,mips24mag,format="A"

count = [16,23,27,29,31]

oplot, alog10(irac3(count)/irac1(count)), alog10(irac4(count)/irac2(count)), psym = 4
oplot, alog10(irac3(count)/irac1(count)), alog10(irac4(count)/irac2(count)), color = bluecolor
xyouts,  alog10(irac3(16)/irac1(16)), alog10(irac4(16)/irac2(16)), ".", color = bluecolor

readlargecol,'/Users/jkrick/ZPHOT/templates/XMM_LSS/Sb_flux_z.dat',z, Ul,U,B,V,R,I,up,gp,rp,ip,zp,J,Ks,irac1,irac2,irac3,irac4,mips24,Ulmag,Umag,Bmag,Vmag,Rmag,Imag,upmag,gpmag,rpmag,ipmag,zpmag,Jmag,Ksmag,irac1mag,irac2mag,irac3mag,irac4mag,mips24mag,format="A"

count = [16,23,27,29,31]

oplot, alog10(irac3(count)/irac1(count)), alog10(irac4(count)/irac2(count)), psym = 4
oplot, alog10(irac3(count)/irac1(count)), alog10(irac4(count)/irac2(count)), color = greencolor
xyouts,  alog10(irac3(16)/irac1(16)), alog10(irac4(16)/irac2(16)), ".", color = greencolor



end
