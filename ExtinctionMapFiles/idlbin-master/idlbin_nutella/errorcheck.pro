pro errorcheck

u = 0.2 /3.
sigu = 12.6
v=56.8 /3.
sigv=13.9

x = u/v

sigx = sqrt(x^2 * ((sigu^2)/(u^2) + (sigv^2)/(v^2)))

print, x, sigx


;alternatively

dxu = ((u+sigu) / v )- x
dxv = abs((u / (v+ sigv)) - x)

sigx = sqrt( (dxu)^2 + (dxv)^2 )

print, x, sigx
end
