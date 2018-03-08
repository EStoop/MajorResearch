pro test_lnp_test

X = [ 1.0, 2.0, 5.0, 7.0, 8.0, 9.0, $
   10.0, 11.0, 12.0, 13.0, 14.0, 15.0, $
   16.0, 17.0, 18.0, 19.0, 20.0, 22.0, $
   23.0, 24.0, 25.0, 26.0, 27.0, 28.0]
Y = [ 0.69502, -0.70425, 0.20632, 0.77206, -0.08339, 0.97806, $
   1.77324, 2.34086, 0.91354, 2.04189, 0.53560, -0.05348, $
   -0.76308, -0.84501, -0.06507, -0.12260, 1.83075, 1.41403, $
   -0.26438, -0.48142, -0.50929, 0.01942, -1.29268, 0.29697]

x = findgen(50)
y = sin(x)

; Test the hypothesis that X and Y represent a significant periodic
; signal against the hypothesis that they represent random noise:
result = LNP_TEST(X, Y, WK1 = wk1, WK2 = wk2, JMAX = jmax)
PRINT, result
print, 'jmax', jmax
a = plot(x, y)

b = plot(1/wk1, wk2)
end
