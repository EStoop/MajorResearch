FUNCTION gauss, X, P
return, P(2)/sqrt(2.*!Pi) * exp(-0.5*((X - P(0))/P(1))^2.)
;return, P(0)/sqrt(2.*!Pi) * exp(-0.5*((X - 0)/0.0009)^2.)

;P(1) = sigma
;P(0) = mean
;P(2) = area
end
