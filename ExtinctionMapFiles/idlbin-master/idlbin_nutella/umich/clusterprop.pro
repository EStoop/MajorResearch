PRO clusterprop

close,/all
z = fltarr(21)
bm = fltarr(21)
rich = fltarr(21)
frac = fltarr(21)


z = [0.150 ,0.023 ,0.143 ,0.029 ,0.084 ,0.183 ,0.077 ,0.077 ,0.228 ,0.415 ,0.004 ,0.005 ,0.192 ,0.166 ,0.165 ,0.171 ,0.103 ,0.066 ,0.004 ,0 ,0.003,0.004]
						
bm = [1.5  ,2  ,1 ,1 ,1.5 ,2.5 ,1.5 ,1  ,1 ,3  ,3 ,1 ,2.5 ,3  ,3 ,2  ,2.5 ,2.5 ,3  ,0 ,0,3 ]	

rich = [2 , 0  , 3  , -1 , 0  , 4  , 0   , 0   , 1  , 0  , 0  , 0  , 2  , 2  , 2  , 2  , 2  , 1  , 0,0,0,0  ]  

frac = [10  ,50 ,13 ,21 ,5  ,0  ,30 ,10  ,2  ,20  ,15  ,40 ,3  ,5  ,5  ,7  ,20  ,20 ,10 ,1,,1.,10]

mydevice = !D.NAME
!p.multi = [0, 0, 2]
SET_PLOT, 'ps'

device, filename = '/n/Godiva7/jkrick/paper/clusterprop.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

plot, z, frac, psym = 2, thick = 3, xtitle= "redshift", ytitle = "fraction ICL", charthick = 3, xthick = 3, ythick = 3
;plot, rich,frac,psym= 2
plot, bm,frac, psym = 2, thick = 3, xtitle = "bautz-morgan type", ytitle = "fraction ICL", charthick = 3, xthick = 3, ythick = 3

device, /close
set_plot, mydevice

end
