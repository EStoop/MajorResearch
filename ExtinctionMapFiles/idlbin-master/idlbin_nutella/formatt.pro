pro formatt
  restore, '/Users/jkrick/nutella/idlbin/objectnew.sav'

  readcol, '/Users/jkrick/selection4matt.txt', n, c, ra, dec, mean_t, std_t, magauto,magerrauto,class_star, format="A"
  mmatch = findobject(ra, dec)

;get rid of the three bad ones
  a = where(mmatch gt 0)
  mmatch_good = mmatch(a)
  n = n(a)
  c = c(a)
  ra = ra(a)
  dec = dec(a)
  mean_t = mean_t(a)
  std_t = std_t(a)
  magauto = magauto(a)
  magerrauto = magerrauto(a)
  class_star = class_star(a)

;run hyperz including the new J band data
;  plothyperz, mmatch_good, '/Users/jkrick/nutella/nep/formatt.ps'

;read in hyperz fit info
  readcol,'/Users/jkrick/nutella/ZPHOT/hyperz_swire_target.z_phot',idz, zphota, chia, proba, specta, nagea, agea, $
          ava, ba, zinf99a,zsup99a,zinf90a,zsup90a,zinf68a,zsup68a,zwma,probwma,Mabsa,$
          zphot2a,prob2a,format="A"

zphotarr = fltarr(n_elements(zphota))
zphot2arr = fltarr(n_elements(zphota))
;print out all the information for matt
  for i = 0, n_elements(n) - 1 do begin
     zphot = objectnew[mmatch_good(i)].zphot
     prob = objectnew[mmatch_good(i)].prob
     specz = objectnew[mmatch_good(i)].specz
     zphot2 = zphota(i)
;     print,  zphot, prob, zphot2, proba(i)

     zphotarr(i) = zphot
     zphot2arr(i) = zphot2
     print, format='(I10, I10, I10, A, F10.5, F10.5, A, F10.5, A, F10.5, F10.5, F10.5, F10.5, F10.2, F10.2, F10.2)',n(i), c(i), mmatch_good(i), "  ", ra(i), dec(i), "  ", mean_t(i), "  ", std_t(i), magauto(i),magerrauto(i),class_star(i),zphot2, proba(i), specz
  endfor
!p.multi = [0, 0, 1]

plothist, zphotarr/zphot2arr, bin = 0.01, xrange = [0.5, 1.5]
end

