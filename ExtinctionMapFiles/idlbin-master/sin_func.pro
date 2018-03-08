FUNCTION sin_func, X, P

return, P(0)*sin(X/P(2) + P(3)) + P(1)

;P(0) = amplitude,
;P(1) = y-intercept
;P(2) = x amplitude = period
;p(3) = phase
end
