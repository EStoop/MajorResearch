pro membersum

close,/all

n = 0
rflux = 0.
vflux = 0.

openr, lun, "/n/Godiva2/jkrick/A141/cmd/member.txt", /get_lun
WHILE (NOT eof(lun)) DO BEGIN
    READF,lun, x, y, rmag, vmag
    IF (sqrt((x - 1652)^2 + (y-1601)^2) LT 860) AND (vmag LT 26.9) $
     AND (rmag LT 27.4) THEN BEGIN
        n = n + 1
        rflux = rflux + 10^((24.6 - rmag)/2.5)
        vflux = vflux + 10^((24.3 - Vmag)/2.5)
    ENDIF
ENDWHILE

print, "rflux, Vflux ", rflux, 100^((4.74 - (24.6 - 2.5*alog10(rflux) - 41))/5), vflux, 100^((4.74 - (24.3 - 2.5*alog10(vflux) - 41))/5)
close, lun
free_lun, lun

END
