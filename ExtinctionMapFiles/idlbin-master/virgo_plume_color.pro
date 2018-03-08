pro virgo_plume_color
;Vega
;Plume A

V_vega = 14.4761
V_vega_noise = 0.804
V_ab = V_vega -0.044
V_ab_noise =V_vega_noise


B_vega = 15.38
B_vega_noise = 1.03
B_ab = B_vega - 0.163
B_ab_noise = B_vega_noise

B_V_vega = B_vega - V_vega
B_V_ab = B_ab - V_ab

print, 'B_V_vega, B_V_ab', B_V_vega, B_V_ab
;all vega
ch1 = 11.58
ch2 = 11.07
ch1_ch2 = ch1 - ch2
ch1_noise = 0.47
ch2_noise = 0.64
ch1_ch2_noise = sqrt(ch1_noise^2 +ch2_noise^2)

;want g- 3.6 

;sdss g' in AB
sdss_g_ab = V_ab + 0.60*(B_V_ab) - 0.12
sdss_g_vega = sdss_g_ab +0.093

g_vega = sdss_g_vega  ; for now all we have is this as a good guess
g_noise = V_vega_noise


g_ch1_vega = g_vega - ch1
g_ch1_vega_noise = sqrt( g_noise^2 + ch1_noise^2)
print, ch1_ch2, '+-', ch1_ch2_noise, g_ch1_vega, '+-', g_ch1_vega_noise

end
