function rms, x


;x = [15.5,15.48,15.32, 15.55, 15.51, 15.44, 15.48, 15.49, 15.5, 15.51]
N = n_elements(x)
mu = (total(x)) / (N)

deviation  = x - mu
deviation2 = deviation^2
sigma2 = (total(deviation2)) /( N - 1)
sigma = sqrt(sigma2)

return, sigma

end
