
PRO diff

colors = GetColor(/Load, Start=1)

max = 10.
df = fltarr(max-1)
dfsum = fltarr(max-1)
x = findgen(max)
x = x+1
y = 10^(-x^(1/4.))

print, x
print, total(x(2:5))
FOR i = 0,max - 2,1 DO BEGIN

    df(i) = y(i) - y(i+1)
    dfsum(i) =  y(i) - (1./(max-1.-i))*(total(y(i:max-1.))) 
ENDFOR
print, dfsum
plot, y, yrange = [-0.1,0.1]
oplot, df, color = colors.blue
oplot, deriv(y), color = colors.yellow
oplot, dfsum , color = colors.orange
end
