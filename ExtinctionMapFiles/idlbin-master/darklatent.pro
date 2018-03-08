pro darklatent
;are the latents in the darks becoming more prevalent?

;0 indicates no latent
;1 indicates latents

l1=[0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,0,0,0,0,-1,-1,0,1,0,0,0,0,0,0,0,1,0,0]
l2=[0,0,1,0,1,0,0,1,0,0,0,1,1,0,1,0,0,0,1,1,0,0,0,0,1,1,0,1,1,1,0,0,0,0,1,0,1,1,-1,-1,1,0,0,0,-1,-1,1,1,1,0,0,1,1,1,1,1,1,1]
;got through including pc20
t = findgen(n_elements(l2)) ; an array representing the 30 dark sets since pc5
!P.multi = [0,1,1]
plot, t, l2, psym = 2, charthick= 1

end
