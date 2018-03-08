;notes from jim's emails


Sure.  It's a structure array named 'PCRS' with the following form:

PCRS[i] = {RELEASED:'',JD_TIME:0d0,CH:0,ROMODE:'',FOV:'',PID:0L,REQKEY:0L,PTITLE:'',REQTITLE:'',ESTDURATION:0.0,$                   TGT_RA:'',TGT_DEC:'',TGT_PMRA:0.0,TGT_PMDEC:0.0,TGT_EPOCH:0.0,TGT_DROW:!values.f_nan,TGT_DCOL:!values.f_nan,$
FRAMETIME:0.0,NFRAMES:0,PU_VMAG:!values.f_nan,PU_RA:'',PU_DEC:'',PU_PMRA:0.0,PU_PMDEC:0.0,$
PU_DRA:!values.f_nan,PU_DDEC:!values.f_nan,PU_EPOCH:0.0,SCLK_START:0D0,SCET_START:'',CAMPAIGN:'',POSANGLE_START:0.0}

You invoke the various elements (eg., peakup V magnitude) using PCRS[*].PU_VMAG.


By the way, you can determine if a measurement is a self-PU or guide star PU where

(pcrs[*].pu_dra EQ 0 AND pcrs[*].pu_ddec EQ 0 )

equals 1.
