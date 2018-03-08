 FUNCTION galpars, galnum,gal,s, t, filearr, j
 
IF (t EQ 5) THEN BEGIN
     gal[galnum].xcenter2 = s[0]
     gal[galnum].ycenter2 = s[1]
     gal[galnum].chi2= s[2]
     gal[galnum].mag2 = s[3]
     gal[galnum].re2 = s[4]
     j = j+1
     s = strsplit(filearr[j], /extract)
     
     gal[galnum].axisratio2 = s[0]
     gal[galnum].pa2 = s[1]
     gal[galnum].bkgd2 = s[2]
     gal[galnum].dsx2= s[2]
     gal[galnum].dsy2 = s[3]
 ENDIF

 IF (t EQ 9) THEN BEGIN
     gal[galnum].xcenter2 = s[0]
     gal[galnum].ycenter2 = s[1]
     gal[galnum].chi2= s[2]
     gal[galnum].mag2 = s[3]
     gal[galnum].re2 = s[4]
     gal[galnum].axisratio2 = s[5]
     gal[galnum].pa2 = s[6]
     gal[galnum].bkgd2 = s[7]
     gal[galnum].dsx2 = s[8]
     j = j+1
     s = strsplit(filearr[j], /extract)
     gal[galnum].dsy2 = s[0]
 ENDIF

 IF (t EQ 10) THEN BEGIN
     gal[galnum].xcenter2 = s[0]
     gal[galnum].ycenter2 = s[1]
     gal[galnum].chi2= s[2]
     gal[galnum].mag2 = s[3]
     gal[galnum].re2 = s[4]
     gal[galnum].axisratio2 = s[5]
     gal[galnum].pa2 = s[6]
     gal[galnum].bkgd2 = s[7]
     gal[galnum].dsx2 = s[8]
     gal[galnum].dsy2= s[9]
 ENDIF

return, 0

END

