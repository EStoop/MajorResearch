pro fukugita_sample
  
;read in the fukugita 2007 et al. sample.  Split out the Ellipticals.
;then put the output in a file for cross-matching with GALEX.


readcol, '/Users/jkrick/Virgo/Galex/UVX/fukugita.csv',   ID      , DR3     , RAdeg   , DEdeg   , EDR     , ISam    , ITar    , T       ,e_T       , Iph     , rmag    , Spec    , z       , Con_z  , Com     , PGC,     format = 'I4    ,A19   ,D10.7  ,F9.5  ,A17   , I1   ,A2    ,F4.1  ,F3.1  ,A2    ,F5.2  ,A13   ,F6.4  ,F5.3  ,A15   ,I5    '

a = where(T eq 0 and z lt 0.1)

print, n_elements(a)

print, RAdeg, DEdeg

openw, outlun, '/Users/jkrick/Virgo/Galex/fukugita_0.txt', /get_lun
for i = 0, n_elements(a) - 1 do begin
   printf, outlun,ID(a(i)),',',  RAdeg(a(i)),',', DEdeg(a(i))
endfor

close, outlun
free_lun, outlun

end

;    1-  4 I4     ---     ID      Running identification number
;   6- 24 A19    ---     DR3     The DR3 photometric identification (1)
;  26- 34 F9.5   deg     RAdeg   Right Ascension in decimal degrees (J2000)
 ; 36- 44 F9.5   deg     DEdeg   Declination in decimal degrees (J2000)
;  46- 62 A17    ---     EDR     The EDR photometric identification (1)
;      64 I1     ---     ISam    Flag for catalog inclusion (2)
;  66- 67 A2     ---     ITar    Flag for spectroscopic target selection (3)
;  69- 72 F4.1   ---     T       Morphology index (4)
;  74- 76 F3.1   ---   e_T       The 1{sigma} uncertainty in T
;  78- 79 A2     ---     Iph     Photometric quality (5)
;  81- 85 F5.2   mag     rmag    Petrosian r magnitude (6)
;  87- 99 A13    ---     Spec    Spectroscopic identification (7)
; 101-106 F6.4   ---     z       Heliocentric redshift
; 108-112 F5.3   ---     Con(z)  Confidence level of redshift measurement
 ;114-128 A15    ---     Com     Additional comments (8)
; 130-134 I5     ---     PGC     ? The PGC catalog number
