pro bias_v2
ps_start, filename= '/Users/jkrick/iwic/noise_darkcal28k.ps'

!p.multi = [0,3,3]
redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
blackcolor = FSC_COLOR("black", !D.Table_Size-4)

exptime = [0.1,0.4,2,6,12,100,200,0.1,0.4,2,6,12,100,200,0.1,0.4,2,6,12,100,200,0.1,0.4,2,6,12,100,200,0.1,0.4,2,6,12,100,200,0.1,0.4,2,6,12,100,200,0.1,0.4,2,6,12,100,200,0.1,0.4,2,6,12,100,200,0.4,12,30,0.1,0.4,2,6,12,100,200,0.1,0.4,2,6,12,100,200,0.1,0.4,2,6,12,100,200,0.1,0.4,2,6,12,100,200,0.4,12,30,0.1,0.4,2,6,12,100,200,0.1,0.4,2,6,12,100,200,0.1,0.4,2,6,12,100,200,0.1,0.4,2,6,12,100,200,0.1,0.4,2,6,12,100,200,0.1,0.4,2,6,12,100,200]

channel = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]

temp = [27.8,27.8,27.8,27.8,27.8,27.8,27.8,27.8,27.8,27.8,27.8,27.8,27.8,27.8,27.8,27.8,27.8,27.8,27.8,27.8,27.8,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28,28.1,28.1,28.1,28.1,28.1,28.1,28.1,28.1,28.1,28.1,28.2,28.2,28.2,28.2,28.2,28.2,28.2,28.2,28.2,28.2,28.2,28.2,28.2,28.2,28.2,28.2,28.2,28.2,28.2,28.2,28.2,28.7,28.7,28.7,28.7,28.7,28.7,28.7,28.7,28.7,28.7,28.7,28.7,28.7,28.7,28.7,28.7,28.7,28.7,28.7,28.7,28.7,28.7,28.7,28.7,28.9,28.9,28.9,28.9,28.9,28.9,28.9,28.9,28.9,28.9,28.9,28.9,28.9,28.9,28.9,28.9,28.9,28.9,28.9,28.9,28.9]

vreset = [-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55]

bias = [450,450,450,450,450,450,450,500,500,500,500,500,500,500,550,550,550,550,550,550,550,500,500,500,500,500,500,500,450,450,450,450,450,450,450,450,450,450,450,450,450,450,550,550,550,550,550,550,550,550,550,550,550,550,550,550,750,750,750,450,450,450,450,450,450,450,450,450,450,450,450,450,450,550,550,550,550,550,550,550,550,550,550,550,550,550,550,500,500,500,600,600,600,600,600,600,600,450,450,450,450,450,450,450,750,750,750,750,750,750,750,450,450,450,450,450,450,450,500,500,500,500,500,500,500,550,550,550,550,550,550,550]

vdduc = [-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.6,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55,-3.55]

vgg1 = [-3.68,-3.68,-3.68,-3.68,-3.68,-3.68,-3.68,-3.68,-3.68,-3.68,-3.68,-3.68,-3.68,-3.68,-3.68,-3.68,-3.68,-3.68,-3.68,-3.68,-3.68,-3.7,-3.7,-3.7,-3.7,-3.7,-3.7,-3.7,-3.7,-3.7,-3.7,-3.7,-3.7,-3.7,-3.7,-3.61,-3.61,-3.61,-3.61,-3.61,-3.61,-3.61,-3.61,-3.61,-3.61,-3.61,-3.61,-3.61,-3.61,-3.61,-3.61,-3.61,-3.61,-3.61,-3.61,-3.61,-3.65,-3.65,-3.65,-3.75,-3.75,-3.75,-3.75,-3.75,-3.75,-3.75,-3.68,-3.68,-3.68,-3.68,-3.68,-3.68,-3.68,-3.68,-3.68,-3.68,-3.68,-3.68,-3.68,-3.68,-3.68,-3.68,-3.68,-3.68,-3.68,-3.68,-3.68,-3.65,-3.65,-3.65,-3.75,-3.75,-3.75,-3.75,-3.75,-3.75,-3.75,-3.75,-3.75,-3.75,-3.75,-3.75,-3.75,-3.75,-3.75,-3.75,-3.75,-3.75,-3.75,-3.75,-3.75,-3.61,-3.61,-3.61,-3.61,-3.61,-3.61,-3.61,-3.61,-3.61,-3.61,-3.61,-3.61,-3.61,-3.61,-3.61,-3.61,-3.61,-3.61,-3.61,-3.61,-3.61]

meannoise = [4.35585,1.52307,0.205295,0.07455,0.04359,0.007355,0.004249,5.3597,1.80651,0.242274,0.0846618,0.049815,0.0083257,0.004834,4.5824,1.48073,0.19828,0.068979,0.0403427,0.0070693,0.0041934,4.9541,1.70037,0.201764,0.066758,0.0415608,0.010937,0.007581,4.8714,1.64912,0.20097,0.067186,0.041589,0.01056,0.0071524,4.81386,1.64091,0.195416,0.06464,0.0401553,0.0104281,0.007138,5.02209,1.7353,0.203555,0.06623,0.041566,0.011458,0.0079478,5.04936,1.73134,0.204145,0.067846,0.042381,0.0114204,0.007872,1.74345,0.04834,0.0223477,4.8391,1.61096,0.22475,0.0831506,0.0491105,0.0077657,0.0043701,4.71015,1.54279,0.214919,0.077314,0.045276,0.007602,0.0043331,4.66335,1.50398,0.2027,0.070401,0.04155,0.0072077,0.004294,4.59043,1.48114,0.20109,0.070666,0.041262,0.0070809,0.004195,1.8212305,0.045053,0.022786,5.09934,1.65615,0.222884,0.078505,0.045598,0.007854,0.0044909,5.146429,1.727408,0.237279,0.086269,0.0498516,0.00771637,0.004323,5.06685,1.60866,0.21524,0.076037,0.046622,0.010209,0.0074856,4.81494,1.64163,0.197985,0.065134,0.0402216,0.010311,0.007084,5.23264,1.784307,0.213445,0.0704467,0.043546,0.011264,0.007855,5.11414,1.74424,0.207972,0.068868,0.042861,0.011521,0.0079034]

sigmanoise = [0.8775,0.2897,0.0345208,0.01189,0.0070414,0.001204,0.0007668,1.2176,0.39405,0.050138,0.0173077,0.010114,0.001665,0.001073,0.91264,0.27843,0.0343,0.0113273,0.006581,0.0011633,0.0007645,1.0247,0.34494,0.039024,0.0121446,0.0069103,0.001753,0.0014104,0.991873,0.34496,0.037105,0.012062,0.006939,0.0017111,0.00131559,1.00127,0.35537,0.036834,0.011006,0.0066697,0.0016973,0.001304,1.064,0.382992,0.0415776,0.0117,0.006827,0.001882,0.0014519,1.0615,0.362469,0.039697,0.012141,0.0071472,0.001863,0.001355,0.3507,0.008659,0.004057,0.90885,0.28227,0.0382,0.013619,0.008271,0.00127666,0.00078825,0.91824,0.279204,0.037066,0.0127899,0.00746,0.001294,0.00077715,0.93506,0.27354,0.035595,0.011014,0.006988,0.001194,0.000783,0.91865,0.27228,0.034882,0.01184,0.0067584,0.001161,0.000764,0.39694,0.00793,0.00382,0.95387,0.30262,0.03757,0.013244,0.007355,0.00147,0.0008345,0.92909,0.294481,0.038185,0.014262,0.008025,0.0012886,0.0008057,0.98344,0.298706,0.038668,0.013085,0.008431,0.00214,0.001946,1.01005,0.344904,0.039149,0.011141,0.006645,0.0017032,0.001238,1.23011,0.42154,0.048799,0.014495,0.008946,0.002112,0.001721,1.08334,0.369018,0.041533,0.0119434,0.007259,0.001811,0.0014072]

noise_darkcal = [5.02648,1.6878,0.227503,0.082612,0.048306,0.00815078,0.004709,5.00823,1.68805,0.226396,0.07911,0.046548,0.0077798,0.004517,5.0954,1.64597,0.2204077,0.076676,0.044845,0.007858,0.0046614,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.8455,1.65513,0.197111,0.065199,0.040502,0.01051867,0.0072302,5.198566,1.796303,0.2107119,0.06856,0.04303,0.011861,0.0082273,5.18175,1.77674,0.209051,0.0696247,0.04349,0.01172,0.008078,0,0,0,0,0,0,0,0,0,0,4.9180578,1.6109,0.22440645,0.08073,0.04727,0.0079381,0.004524,4.9092,1.5833,0.21339,0.07411,0.04374,0.007587,0.00452,5.10188,1.64614,0.223488,0.07854,0.045859,0.00787,0.00466,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4.8781,1.66318,0.200583,0.06596,0.04075,0.010446,0.007177,4.8793,1.66377,0.19903,0.06569,0.040604,0.010503,0.0073247,5.1263,1.74843,0.20847,0.069033,0.042964,0.011549,0.007922]

sigma_darkcal = [0.97244,0.32109,0.038255,0.013177,0.007803,0.001334,0.0008498,1.13778,0.36821,0.043838,0.01617,0.00945,0.0015558,0.001003,1.0162,0.3095,0.03813,0.01259,0.007315,0.00129,0.0008499,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.02004,0.358408,0.037153,0.011101,0.006758,0.001712,0.001315,1.10121,0.39645,0.04304,0.012111,0.00707,0.001948,0.001503,1.0893,0.37197,0.04074,0.0121595,0.007336,0.0019119,0.00139,0,0,0,0,0,0,0,0,0,0,0.958794,0.29153,0.038703,0.01335,0.0077889,0.0013513,0.0008115,0.9844,0.28796,0.03747,0.011594,0.007356,0.001257,0.000824,1.0212,0.30262,0.038764,0.013155,0.00751,0.00129,0.0009496,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1.02329,0.34943,0.039663,0.01129,0.006732,0.001726,0.001254,1.146999,0.39306,0.045503,0.013516,0.00834,0.001969,0.001605,1.0858,0.36988,0.04163,0.01197,0.007276,0.001815,0.001411]



noisy4x=[0,0,0,0,0,0.0010884,0,0,0,0,0,0,0.00201055,0,0,0,0,0,0,0.00291072,0,0,0,0,0,0,0.00128086,0,0,0,0,0,0,0.000823243,0,0,0,0,0,0,0.000704851,0,0,0,0,0,0,0.00176742,0,0,0,0,0,0,0.00178627,0,0,0,0,0,0,0,0,0,0.00122142,0,0,0,0,0,0,0.00123833,0,0,0,0,0,0,0.00373557,0,0,0,0,0,0,0.0037872,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.00103318,0,0,0,0,0,0,0.0020081,0,0,0,0,0,0,0.00251804,0]

latent_strength=[0.00226049,0.00226049,0.00226049,0.00226049,0.00226049,0.00226049,0.00226049,0.00235513,0.00235513,0.00235513,0.00235513,0.00235513,0.00235513,0.00235513,0.00225598,0.00225598,0.00225598,0.00225598,0.00225598,0.00225598,0.00225598,0.00145,0.00145,0.00145,0.00145,0.00145,0.00145,0.00145,0.00149077,0.00149077,0.00149077,0.00149077,0.00149077,0.00149077,0.00149077,0.0065591,0.0065591,0.0065591,0.0065591,0.0065591,0.0065591,0.0065591,0.00195328,0.00195328,0.00195328,0.00195328,0.00195328,0.00195328,0.00195328,0.00195688,0.00195688,0.00195688,0.00195688,0.00195688,0.00195688,0.00195688,0,0,0,0.00149077,0.00149077,0.00149077,0.00149077,0.00149077,0.00149077,0.00149077,0.0065591,0.0065591,0.0065591,0.0065591,0.0065591,0.0065591,0.0065591,0.00195328,0.00195328,0.00195328,0.00195328,0.00195328,0.00195328,0.00195328,0.00195688,0.00195688,0.00195688,0.00195688,0.00195688,0.00195688,0.00195688,0,0,0,0.00221165,0.00221165,0.00221165,0.00221165,0.00221165,0.00221165,0.00221165,0.00204855,0.00204855,0.00204855,0.00204855,0.00204855,0.00204855,0.00204855,0.00170674,0.00170674,0.00170674,0.00170674,0.00170674,0.00170674,0.00170674,0.00226049,0.00226049,0.00226049,0.00226049,0.00226049,0.00226049,0.00226049,0.00235513,0.00235513,0.00235513,0.00235513,0.00235513,0.00235513,0.00235513,0.00225598,0.00225598,0.00225598,0.00225598,0.00225598,0.00225598,0.00225598]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;looking at bias
;ch1

ch1t28r355g368d355_0p4s = where(temp ge 28.0 and temp le 28.2 and vreset eq -3.55 and vdduc eq -3.55 and vgg1 eq -3.68 and exptime eq 0.4 and channel eq 1)
ch1t28r355g368d355_12s = where(temp ge 28.0 and temp le 28.2 and vreset eq -3.55 and vdduc eq -3.55 and vgg1 eq -3.68 and exptime eq 12 and channel eq 1)
ch1t28r355g368d355_100s = where(temp ge 28.0 and temp le 28.2 and vreset eq -3.55 and vdduc eq -3.55 and vgg1 eq -3.68 and exptime eq 100 and channel eq 1)
;ch1t29r355g368d355_12s = where(temp ge 28.7 and temp le 28.9 and vreset eq -3.55 and vdduc eq -3.55 and vgg1 eq -3.68 and exptime eq 12 and channel eq 1)
;ch1t29r355g368d355_100s = where(temp ge 28.7 and temp le 28.9 and vreset eq -3.55 and vdduc eq -3.55 and vgg1 eq -3.68 and exptime eq 100 and channel eq 1)


;-----
;ch2
ch2t28r355g368d355_0p4s = where(temp ge 28.0 and temp le 28.2 and vreset eq -3.55 and vdduc eq -3.55 and vgg1 eq -3.61 and exptime eq 0.4 and channel eq 2)
ch2t28r355g368d355_12s = where(temp ge 28.0 and temp le 28.2 and vreset eq -3.55 and vdduc eq -3.55 and vgg1 eq -3.61 and exptime eq 12 and channel eq 2)
ch2t28r355g368d355_100s = where(temp ge 28.0 and temp le 28.2 and vreset eq -3.55 and vdduc eq -3.55 and vgg1 eq -3.61 and exptime eq 100 and channel eq 2)

ch2t29r355g368d355_0p4s = where(temp ge 28.7 and temp le 28.9 and vreset eq -3.55 and vdduc eq -3.55 and vgg1 eq -3.61 and exptime eq 0.4 and channel eq 2)
ch2t29r355g368d355_12s = where(temp ge 28.7 and temp le 28.9 and vreset eq -3.55 and vdduc eq -3.55 and vgg1 eq -3.61 and exptime eq 12 and channel eq 2)
ch2t29r355g368d355_100s = where(temp ge 28.7 and temp le 28.9 and vreset eq -3.55 and vdduc eq -3.55 and vgg1 eq -3.61 and exptime eq 100 and channel eq 2)

;---------
;plots
;0.4s
plot, bias(ch1t28r355g368d355_0p4s), meannoise(ch1t28r355g368d355_0p4s), psym = 4, xtitle = 'bias', ytitle = 'mean noise(MJy/sr)', title = '0.4s r355g368_361d355',yrange=[1.4,2.0], charsize=1.5
oplot, bias(ch1t28r355g368d355_0p4s), meannoise(ch1t28r355g368d355_0p4s)

oplot, bias(ch2t28r355g368d355_0p4s ), meannoise(ch2t28r355g368d355_0p4s ), psym = 4, color = bluecolor
oplot, bias(ch2t28r355g368d355_0p4s ), meannoise(ch2t28r355g368d355_0p4s ), color = bluecolor
oplot, bias(ch2t29r355g368d355_0p4s ), meannoise(ch2t29r355g368d355_0p4s ), psym = 2, color =bluecolor
oplot, bias(ch2t29r355g368d355_0p4s ), meannoise(ch2t29r355g368d355_0p4s ), color = bluecolor
;----
plot, bias(ch1t28r355g368d355_0p4s), noise_darkcal(ch1t28r355g368d355_0p4s), psym = 4, xtitle = 'bias', ytitle = 'mean noise(MJy/sr)', title = '0.4s r355g368_361d355',yrange=[1.4,2.0], charsize=1.5
oplot, bias(ch1t28r355g368d355_0p4s), noise_darkcal(ch1t28r355g368d355_0p4s)

oplot, bias(ch2t28r355g368d355_0p4s ), noise_darkcal(ch2t28r355g368d355_0p4s ), psym = 4, color = bluecolor
oplot, bias(ch2t28r355g368d355_0p4s ), noise_darkcal(ch2t28r355g368d355_0p4s ), color = bluecolor
oplot, bias(ch2t29r355g368d355_0p4s ), noise_darkcal(ch2t29r355g368d355_0p4s ), psym = 2, color =bluecolor
oplot, bias(ch2t29r355g368d355_0p4s ), noise_darkcal(ch2t29r355g368d355_0p4s ), color = bluecolor


plot, bias(ch1t28r355g368d355_0p4s), noisy4x(ch1t28r355g368d355_0p4s), psym = 4, xtitle = 'bias', ytitle = 'fraction of noisy pix 4x zodi', title = '0.4s r355g368_361d355' , charsize=1.5, /ylog, yrange=[.001, .01]
oplot, bias(ch1t28r355g368d355_0p4s), noisy4x(ch1t28r355g368d355_0p4s)

oplot, bias(ch2t28r355g368d355_0p4s ), noisy4x(ch2t28r355g368d355_0p4s ), psym = 4, color = bluecolor
oplot, bias(ch2t28r355g368d355_0p4s ), noisy4x(ch2t28r355g368d355_0p4s ), color = bluecolor
oplot, bias(ch2t29r355g368d355_0p4s ), noisy4x(ch2t29r355g368d355_0p4s ), psym = 2, color =bluecolor
oplot, bias(ch2t29r355g368d355_0p4s ), noisy4x(ch2t29r355g368d355_0p4s ), color = bluecolor


;plot, bias(ch1t28r355g368d355_0p4s), latent_strength(ch1t28r355g368d355_0p4s), psym = 4, xtitle = 'bias', ytitle = 'latent strength', title = '0.4s r355g368_361d355' , charsize=1.5, /ylog
;oplot, bias(ch1t28r355g368d355_0p4s), latent_strength(ch1t28r355g368d355_0p4s)

;oplot, bias(ch2t28r355g368d355_0p4s ), latent_strength(ch2t28r355g368d355_0p4s ), psym = 4, color = bluecolor
;oplot, bias(ch2t28r355g368d355_0p4s ), latent_strength(ch2t28r355g368d355_0p4s ), color = bluecolor
;oplot, bias(ch2t29r355g368d355_0p4s ), latent_strength(ch2t29r355g368d355_0p4s ), psym = 2, color =bluecolor
;oplot, bias(ch2t29r355g368d355_0p4s ), latent_strength(ch2t29r355g368d355_0p4s ), color = bluecolor

;12s
plot, bias(ch1t28r355g368d355_12s), meannoise(ch1t28r355g368d355_12s), psym = 4, xtitle = 'bias', ytitle = 'mean noise(MJy/sr)', title = '12s r355g368_361d355', yrange=[.03,.06], charsize=1.5
oplot, bias(ch1t28r355g368d355_12s), meannoise(ch1t28r355g368d355_12s)

oplot, bias(ch2t28r355g368d355_12s ), meannoise(ch2t28r355g368d355_12s ), psym = 4, color = bluecolor
oplot, bias(ch2t28r355g368d355_12s ), meannoise(ch2t28r355g368d355_12s ), color = bluecolor
oplot, bias(ch2t29r355g368d355_12s ), meannoise(ch2t29r355g368d355_12s ), psym = 2, color =bluecolor
oplot, bias(ch2t29r355g368d355_12s ), meannoise(ch2t29r355g368d355_12s ), color = bluecolor
;-----
plot, bias(ch1t28r355g368d355_12s), noise_darkcal(ch1t28r355g368d355_12s), psym = 4, xtitle = 'bias', ytitle = 'mean noise(MJy/sr)', title = '12s r355g368_361d355', yrange=[.03,.06], charsize=1.5
oplot, bias(ch1t28r355g368d355_12s), noise_darkcal(ch1t28r355g368d355_12s)

oplot, bias(ch2t28r355g368d355_12s ), noise_darkcal(ch2t28r355g368d355_12s ), psym = 4, color = bluecolor
oplot, bias(ch2t28r355g368d355_12s ), noise_darkcal(ch2t28r355g368d355_12s ), color = bluecolor
oplot, bias(ch2t29r355g368d355_12s ), noise_darkcal(ch2t29r355g368d355_12s ), psym = 2, color =bluecolor
oplot, bias(ch2t29r355g368d355_12s ), noise_darkcal(ch2t29r355g368d355_12s ), color = bluecolor


plot, bias(ch1t28r355g368d355_12s), noisy4x(ch1t28r355g368d355_12s), psym = 4, xtitle = 'bias', ytitle = 'fraction of noisy pix 4x zodi', title = '12s r355g368_361d355',  charsize=1.5, /ylog, yrange=[.001, .01]
oplot, bias(ch1t28r355g368d355_12s), noisy4x(ch1t28r355g368d355_12s)

oplot, bias(ch2t28r355g368d355_12s ), noisy4x(ch2t28r355g368d355_12s ), psym = 4, color = bluecolor
oplot, bias(ch2t28r355g368d355_12s ), noisy4x(ch2t28r355g368d355_12s ), color = bluecolor
oplot, bias(ch2t29r355g368d355_12s ), noisy4x(ch2t29r355g368d355_12s ), psym = 2, color =bluecolor
oplot, bias(ch2t29r355g368d355_12s ), noisy4x(ch2t29r355g368d355_12s ), color = bluecolor
legend, ['28k', '29k'], psym=[4,2],/center,/bottom
legend, ['ch1', 'ch2'], linestyle = [0,0], color = [blackcolor,bluecolor], /center,/top

;plot, bias(ch1t28r355g368d355_12s), latent_strength(ch1t28r355g368d355_12s), psym = 4, xtitle = 'bias', ytitle = 'latent strength', title = '12s r355g368_361d355',  charsize=1.5, /ylog
;oplot, bias(ch1t28r355g368d355_12s), latent_strength(ch1t28r355g368d355_12s)

;oplot, bias(ch2t28r355g368d355_12s ), latent_strength(ch2t28r355g368d355_12s ), psym = 4, color = bluecolor
;oplot, bias(ch2t28r355g368d355_12s ), latent_strength(ch2t28r355g368d355_12s ), color = bluecolor
;oplot, bias(ch2t29r355g368d355_12s ), latent_strength(ch2t29r355g368d355_12s ), psym = 2, color =bluecolor
;oplot, bias(ch2t29r355g368d355_12s ), latent_strength(ch2t29r355g368d355_12s ), color = bluecolor
;xyouts, 450, 0.007,'CR?'

;-------
;100s
plot, bias(ch1t28r355g368d355_100s), meannoise(ch1t28r355g368d355_100s), psym = 4, xtitle = 'bias', ytitle = 'mean noise(MJy/sr)', title = '100s r355g368_361d355', yrange=[.006,.012], charsize=1.5
oplot, bias(ch1t28r355g368d355_100s), meannoise(ch1t28r355g368d355_100s)

oplot, bias(ch2t28r355g368d355_100s ), meannoise(ch2t28r355g368d355_100s ), psym = 4, color = bluecolor
oplot, bias(ch2t28r355g368d355_100s ), meannoise(ch2t28r355g368d355_100s ), color = bluecolor
oplot, bias(ch2t29r355g368d355_100s ), meannoise(ch2t29r355g368d355_100s ), psym = 2, color = bluecolor
oplot, bias(ch2t29r355g368d355_100s ), meannoise(ch2t29r355g368d355_100s ), color = bluecolor
;----
plot, bias(ch1t28r355g368d355_100s), meannoise(ch1t28r355g368d355_100s), psym = 4, xtitle = 'bias', ytitle = 'mean noise(MJy/sr)', title = '100s r355g368_361d355', yrange=[.006,.012], charsize=1.5
oplot, bias(ch1t28r355g368d355_100s), meannoise(ch1t28r355g368d355_100s)

oplot, bias(ch2t28r355g368d355_100s ), meannoise(ch2t28r355g368d355_100s ), psym = 4, color = bluecolor
oplot, bias(ch2t28r355g368d355_100s ), meannoise(ch2t28r355g368d355_100s ), color = bluecolor
oplot, bias(ch2t29r355g368d355_100s ), meannoise(ch2t29r355g368d355_100s ), psym = 2, color = bluecolor
oplot, bias(ch2t29r355g368d355_100s ), meannoise(ch2t29r355g368d355_100s ), color = bluecolor

plot, bias(ch1t28r355g368d355_100s), noisy4x(ch1t28r355g368d355_100s), psym = 4, xtitle = 'bias', ytitle = 'fraction of noisy pix 4x zodi', title = '100s r355g368_361d355',  charsize=1.5, /ylog, yrange=[.0001, .01]
oplot, bias(ch1t28r355g368d355_100s), noisy4x(ch1t28r355g368d355_100s)

oplot, bias(ch2t28r355g368d355_100s ), noisy4x(ch2t28r355g368d355_100s ), psym = 4, color = bluecolor
oplot, bias(ch2t28r355g368d355_100s ), noisy4x(ch2t28r355g368d355_100s ), color = bluecolor
oplot, bias(ch2t29r355g368d355_100s ), noisy4x(ch2t29r355g368d355_100s ), psym = 2, color = bluecolor
oplot, bias(ch2t29r355g368d355_100s ), noisy4x(ch2t29r355g368d355_100s ), color = bluecolor

;plot, bias(ch1t28r355g368d355_100s), latent_strength(ch1t28r355g368d355_100s), psym = 4, xtitle = 'bias', ytitle = 'latent strength', title = '100s r355g368_361d355',  charsize=1.5, /ylog
;oplot, bias(ch1t28r355g368d355_100s), latent_strength(ch1t28r355g368d355_100s)

;oplot, bias(ch2t28r355g368d355_100s ), latent_strength(ch2t28r355g368d355_100s ), psym = 4, color = bluecolor
;oplot, bias(ch2t28r355g368d355_100s ), latent_strength(ch2t28r355g368d355_100s ), color = bluecolor
;oplot, bias(ch2t29r355g368d355_100s ), latent_strength(ch2t29r355g368d355_100s ), psym = 2, color = bluecolor
;oplot, bias(ch2t29r355g368d355_100s ), latent_strength(ch2t29r355g368d355_100s ), color = bluecolor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;looking at vreset and vdduc
;ch1

ch1t28b55g368_0p4s = where(temp ge 28.0 and temp le 28.2 and  vgg1 eq -3.68 and exptime eq 0.4 and channel eq 1 and bias eq 550)
ch1t28b55g368_12s = where(temp ge 28.0 and temp le 28.2 and  vgg1 eq -3.68 and exptime eq 12 and channel eq 1 and bias eq 550)
ch1t28b55g368_100s = where(temp ge 28.0 and temp le 28.2 and vgg1 eq -3.68 and exptime eq 100 and channel eq 1 and bias eq 550)


;-----
;ch2
ch2t28b55g368_0p4s = where(temp ge 28.0 and temp le 28.2 and  vgg1 eq -3.61 and exptime eq 0.4 and channel eq 2 and bias eq 550)
ch2t28b55g368_12s = where(temp ge 28.0 and temp le 28.2 and vgg1 eq -3.61 and exptime eq 12 and channel eq 2 and bias eq 550)
ch2t28b55g368_100s = where(temp ge 28.0 and temp le 28.2 and  vgg1 eq -3.61 and exptime eq 100 and channel eq 2 and bias eq 550)


;---------
;plots
;0.4s
plot, vdduc(ch1t28b55g368_0p4s), meannoise(ch1t28b55g368_0p4s), psym = 4, xtitle = 'vdduc & vreset', ytitle = 'mean noise(MJy/sr)', title = '0.4s b55g368',yrange=[1.4,2.0], charsize=1.5,xrange=[-3.65,-3.50]
oplot, vdduc(ch1t28b55g368_0p4s), meannoise(ch1t28b55g368_0p4s)

oplot, vdduc(ch2t28b55g368_0p4s ), meannoise(ch2t28b55g368_0p4s ), psym = 4, color = bluecolor
oplot, vdduc(ch2t28b55g368_0p4s ), meannoise(ch2t28b55g368_0p4s ), color = bluecolor

plot, vdduc(ch1t28b55g368_0p4s), noisy4x(ch1t28b55g368_0p4s), psym = 4, xtitle = 'vdduc & vreset', ytitle = 'fraction of noisy pix 4x zodi', title = '0.4s b55g368', charsize=1.5, /ylog, yrange=[.001, .01],xrange=[-3.65,-3.50]
oplot, vdduc(ch1t28b55g368_0p4s), noisy4x(ch1t28b55g368_0p4s)

oplot, vdduc(ch2t28b55g368_0p4s ), noisy4x(ch2t28b55g368_0p4s ), psym = 4, color = bluecolor
oplot, vdduc(ch2t28b55g368_0p4s ), noisy4x(ch2t28b55g368_0p4s ), color = bluecolor

plot, vdduc(ch1t28b55g368_0p4s), latent_strength(ch1t28b55g368_0p4s), psym = 4, xtitle = 'vdduc & vreset', ytitle = 'latent strength', title = '0.4s b55g368', charsize=1.5, /ylog,xrange=[-3.65,-3.50]
oplot, vdduc(ch1t28b55g368_0p4s), latent_strength(ch1t28b55g368_0p4s)

oplot, vdduc(ch2t28b55g368_0p4s ), latent_strength(ch2t28b55g368_0p4s ), psym = 4, color = bluecolor
oplot, vdduc(ch2t28b55g368_0p4s ), latent_strength(ch2t28b55g368_0p4s ), color = bluecolor

;12s
plot, vdduc(ch1t28b55g368_12s), meannoise(ch1t28b55g368_12s), psym = 4, xtitle = 'vdduc & vreset', ytitle = 'mean noise(MJy/sr)', title = '12s b55g368', yrange=[.03,.06], charsize=1.5,xrange=[-3.65,-3.50]
oplot, vdduc(ch1t28b55g368_12s), meannoise(ch1t28b55g368_12s)

oplot, vdduc(ch2t28b55g368_12s ), meannoise(ch2t28b55g368_12s ), psym = 4, color = bluecolor
oplot, vdduc(ch2t28b55g368_12s ), meannoise(ch2t28b55g368_12s ), color = bluecolor

plot, vdduc(ch1t28b55g368_12s), noisy4x(ch1t28b55g368_12s), psym = 4, xtitle = 'vdduc & vreset', ytitle = 'fraction of noisy pix 4x zodi', title = '12s b55g368',  charsize=1.5, /ylog, yrange=[.001, .01],xrange=[-3.65,-3.50]
oplot, vdduc(ch1t28b55g368_12s), noisy4x(ch1t28b55g368_12s)

oplot, vdduc(ch2t28b55g368_12s ), noisy4x(ch2t28b55g368_12s ), psym = 4, color = bluecolor
oplot, vdduc(ch2t28b55g368_12s ), noisy4x(ch2t28b55g368_12s ), color = bluecolor
legend, ['28k', '29k'], psym=[4,2],/center,/bottom
legend, ['ch1', 'ch2'], linestyle = [0,0], color = [blackcolor,bluecolor], /center,/top

plot, vdduc(ch1t28b55g368_12s), latent_strength(ch1t28b55g368_12s), psym = 4, xtitle = 'vdduc & vreset', ytitle = 'latent strength', title = '12s b55g368',  charsize=1.5, /ylog,xrange=[-3.65,-3.50]
oplot, vdduc(ch1t28b55g368_12s), latent_strength(ch1t28b55g368_12s)

oplot, vdduc(ch2t28b55g368_12s ), latent_strength(ch2t28b55g368_12s ), psym = 4, color = bluecolor
oplot, vdduc(ch2t28b55g368_12s ), latent_strength(ch2t28b55g368_12s ), color = bluecolor
legend, ['28k', '29k'], psym=[4,2],/center,/bottom
legend, ['ch1', 'ch2'], linestyle = [0,0], color = [blackcolor,bluecolor], /center,/top

;-------
;100s
plot, vdduc(ch1t28b55g368_100s), meannoise(ch1t28b55g368_100s), psym = 4, xtitle = 'vdduc & vreset', ytitle = 'mean noise(MJy/sr)', title = '100s b55g368', yrange=[.006,.012], charsize=1.5,xrange=[-3.65,-3.50]
oplot, vdduc(ch1t28b55g368_100s), meannoise(ch1t28b55g368_100s)

oplot, vdduc(ch2t28b55g368_100s ), meannoise(ch2t28b55g368_100s ), psym = 4, color = bluecolor
oplot, vdduc(ch2t28b55g368_100s ), meannoise(ch2t28b55g368_100s ), color = bluecolor

plot, vdduc(ch1t28b55g368_100s), noisy4x(ch1t28b55g368_100s), psym = 4, xtitle = 'vdduc & vreset', ytitle = 'fraction of noisy pix 4x zodi', title = '100s b55g368',  charsize=1.5, /ylog, yrange=[.0001, .01],xrange=[-3.65,-3.50]
oplot, vdduc(ch1t28b55g368_100s), noisy4x(ch1t28b55g368_100s)

oplot, vdduc(ch2t28b55g368_100s ), noisy4x(ch2t28b55g368_100s ), psym = 4, color = bluecolor
oplot, vdduc(ch2t28b55g368_100s ), noisy4x(ch2t28b55g368_100s ), color = bluecolor

plot, vdduc(ch1t28b55g368_100s), latent_strength(ch1t28b55g368_100s), psym = 4, xtitle = 'vdduc & vreset', ytitle = 'latent strength', title = '100s b55g368',  charsize=1.5, /ylog,xrange=[-3.65,-3.50],/nodata
oplot, vdduc(ch1t28b55g368_100s), latent_strength(ch1t28b55g368_100s)

oplot, vdduc(ch2t28b55g368_100s ), latent_strength(ch2t28b55g368_100s ), psym = 4, color = bluecolor
oplot, vdduc(ch2t28b55g368_100s ), latent_strength(ch2t28b55g368_100s ), color = bluecolor



ps_end


end
