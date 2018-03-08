pro zcat

close, /all
;device, true=24
;device, decomposed=0
;colors = GetColor(/load, Start=1)

;!p.multi = [0, 0, 1]
;ps_open, file = "/Users/jkrick/ZPHOT/zhist.swire.ps", /portrait, xsize = 4, ysize = 4,/color


restore, '/Users/jkrick/idlbin/objectnew.sav'

readcol,'/Users/jkrick/ZPHOT/hyperz#.z_phot',idz, zphota, chia, proba, specta, nagea, agea, $
        ava, ba, zinf99a,zsup99a,zinf90a,zsup90a,zinf68a,zsup68a,zwma,probwma,Mabsa,$
        zphot2a,prob2a,format="A"

;plot, chia, proba, psym = 3, xtitle = "chia", ytitle = "proba", xrange = [0,10]
plothist, zphota, xrange=[0,6], bin=0.1

readlargecol, '/Users/jkrick/ZPHOT/hyperz#.obs_sed',idz, u,g,r,i,acs,z,j,wj,tj,wh,th,wk,tk,one,two,three, four,$
         uerr, gerr, rerr, ierr,zerr,acserr, jerr, wjerr,tjerr,wherr,therr,wkerr,tkerr,oneerr, twoerr, threeerr, fourerr,format="A",/debug

print, n_elements(zphota), n_elements(u), n_elements(objectnew.rmag)
for count = 0l, n_elements(zphota) - 1 do begin
;   print, count
   objectnew[idz(count)].zphot = zphota(count)
   objectnew[idz(count)].chisq = chia(count)
   objectnew[idz(count)].prob = proba(count)
   objectnew[idz(count)].spt = specta(count)
   objectnew[idz(count)].age = nagea(count)
   objectnew[idz(count)].av = ava(count)
   objectnew[idz(count)].Mabs = Mabsa(count) 
;   objectnew[idz(count)].zphot = zphota(count)
;   objectnew[idz(count)].chisq = chia(count)
;   objectnew[idz(count)].prob = proba(count)
;   objectnew[idz(count)].spt = specta(count)
;   objectnew[idz(count)].age = nagea(count)
;   objectnew[idz(count)].av = ava(count)
;   objectnew[idz(count)].Mabs = Mabsa(count) 
endfor


;-----------------------------------------------------------------------

;ps_close, /noprint, /noid


;now it has the photo-z info in it.
save, objectnew, filename='/Users/jkrick/idlbin/objectnew.sav'

END

