pro make_bc03_filter
;ch1
f = [0.000490,0.000586,0.000453,0.000458,0.000439,0.000531,0.000641,0.000442,0.000535,0.000655,0.000525,0.000631,0.000618,0.000496,0.000700,0.000895,0.000873,0.001045,0.001294,0.001390,0.001683,0.002132,0.002413,0.002782,0.003343,0.003846,0.004596,0.005760,0.006889,0.008313,0.010170,0.012180,0.014820,0.018160,0.021810,0.026460,0.031970,0.038570,0.046500,0.055220,0.065330,0.077400,0.091470,0.106100,0.122000,0.140600,0.158900,0.178800,0.200200,0.220500,0.241200,0.261800,0.282200,0.298600,0.314000,0.330400,0.343100,0.352900,0.360600,0.368500,0.372300,0.373700,0.375000,0.374400,0.371800,0.370300,0.368600,0.365300,0.361900,0.358000,0.354100,0.351300,0.349500,0.348400,0.345200,0.343700,0.344600,0.343600,0.343600,0.344300,0.344800,0.347400,0.349600,0.351200,0.354100,0.357800,0.362300,0.366100,0.369800,0.374800,0.379500,0.384100,0.389500,0.393600,0.397400,0.400700,0.403400,0.406500,0.410400,0.412500,0.413700,0.415800,0.415800,0.416300,0.416500,0.415300,0.414500,0.414300,0.414500,0.413400,0.412500,0.411700,0.409100,0.408600,0.409100,0.406900,0.406000,0.405600,0.405000,0.405500,0.404400,0.402400,0.402700,0.403600,0.404200,0.404700,0.404300,0.404700,0.406600,0.406800,0.406600,0.407000,0.406700,0.406500,0.406700,0.407200,0.407000,0.406600,0.406300,0.405100,0.405100,0.404900,0.403900,0.404800,0.404100,0.401600,0.403300,0.405600,0.405600,0.406300,0.407100,0.408000,0.410900,0.413500,0.416000,0.418700,0.420800,0.424500,0.427300,0.429100,0.432100,0.433300,0.433500,0.435300,0.436700,0.435400,0.434400,0.435600,0.435100,0.434300,0.434000,0.430300,0.426500,0.425400,0.424300,0.422100,0.420200,0.417900,0.416300,0.415700,0.414300,0.413900,0.414000,0.413100,0.413400,0.414500,0.415800,0.417600,0.418600,0.419400,0.421800,0.424800,0.427100,0.430500,0.434300,0.436900,0.440200,0.443400,0.445100,0.447700,0.450600,0.452900,0.455700,0.457200,0.456400,0.456400,0.458300,0.459900,0.460500,0.461400,0.461500,0.461500,0.461600,0.461500,0.462100,0.462900,0.462000,0.459600,0.458000,0.458200,0.457100,0.455500,0.455500,0.454700,0.453700,0.452800,0.451800,0.452100,0.452600,0.452600,0.451500,0.450500,0.450300,0.450800,0.452200,0.452800,0.452800,0.454400,0.456200,0.457400,0.458600,0.459500,0.460200,0.460200,0.459600,0.459500,0.458800,0.459100,0.459700,0.458400,0.457600,0.456800,0.455800,0.455300,0.453600,0.450900,0.448700,0.446900,0.445600,0.444800,0.442500,0.440600,0.440100,0.438100,0.436900,0.437600,0.436400,0.435200,0.435500,0.435700,0.436700,0.438200,0.438900,0.439000,0.440500,0.444100,0.446000,0.446800,0.449500,0.452500,0.455100,0.457200,0.459000,0.460200,0.462300,0.465800,0.468800,0.470200,0.470500,0.470500,0.471000,0.472100,0.473200,0.473400,0.472700,0.473000,0.472800,0.472800,0.473200,0.472300,0.471000,0.469700,0.469400,0.469700,0.467500,0.464700,0.463700,0.463500,0.462200,0.460400,0.460100,0.459600,0.458400,0.456700,0.456800,0.457500,0.456800,0.457100,0.457100,0.455600,0.456400,0.456900,0.454100,0.453000,0.452900,0.451100,0.449500,0.447000,0.444300,0.443000,0.441300,0.438500,0.437400,0.435200,0.431700,0.430200,0.429400,0.428000,0.428000,0.428000,0.428200,0.430100,0.433000,0.436400,0.439600,0.442300,0.444500,0.446100,0.447100,0.448300,0.447300,0.443700,0.434700,0.421100,0.405000,0.383500,0.355900,0.328600,0.293700,0.256300,0.221500,0.186900,0.152200,0.123700,0.098900,0.077500,0.058690,0.044850,0.033730,0.025150,0.018770,0.013760,0.009610,0.006888,0.004779,0.002983,0.002031,0.001633,0.001354,0.001214,0.001186,0.001077,0.001128,0.001287,0.001269,0.001177,0.001125,0.001074,0.000936,0.000751,0.000590]


l =[3.081060,3.082890,3.084730,3.086560,3.088400,3.090240,3.092080,3.093930,3.095780,3.097620,3.099480,3.101330,3.103190,3.105040,3.106900,3.108770,3.110630,3.112500,3.114370,3.116240,3.118110,3.119990,3.121870,3.123750,3.125630,3.127520,3.129400,3.131290,3.133190,3.135080,3.136980,3.138880,3.140780,3.142680,3.144590,3.146490,3.148410,3.150320,3.152230,3.154150,3.156070,3.157990,3.159920,3.161840,3.163770,3.165700,3.167640,3.169570,3.171510,3.173450,3.175400,3.177340,3.179290,3.181240,3.183190,3.185150,3.187110,3.189070,3.191030,3.192990,3.194960,3.196930,3.198900,3.200880,3.202860,3.204840,3.206820,3.208800,3.210790,3.212780,3.214770,3.216760,3.218760,3.220760,3.222760,3.224770,3.226770,3.228780,3.230790,3.232810,3.234820,3.236840,3.238870,3.240890,3.242920,3.244950,3.246980,3.249010,3.251050,3.253090,3.255130,3.257180,3.259220,3.261270,3.263320,3.265380,3.267440,3.269500,3.271560,3.273630,3.275690,3.277760,3.279840,3.281910,3.283990,3.286070,3.288160,3.290240,3.292330,3.294420,3.296520,3.298620,3.300710,3.302820,3.304920,3.307030,3.309140,3.311250,3.313370,3.315490,3.317610,3.319730,3.321860,3.323990,3.326120,3.328260,3.330390,3.332530,3.334680,3.336820,3.338970,3.341120,3.343280,3.345430,3.347590,3.349760,3.351920,3.354090,3.356260,3.358430,3.360610,3.362790,3.364970,3.367160,3.369350,3.371540,3.373730,3.375930,3.378130,3.380330,3.382530,3.384740,3.386950,3.389170,3.391380,3.393600,3.395820,3.398050,3.400280,3.402510,3.404740,3.406980,3.409220,3.411460,3.413710,3.415960,3.418210,3.420460,3.422720,3.424980,3.427250,3.429510,3.431780,3.434060,3.436330,3.438610,3.440890,3.443180,3.445460,3.447750,3.450050,3.452350,3.454650,3.456950,3.459260,3.461560,3.463880,3.466190,3.468510,3.470830,3.473160,3.475480,3.477820,3.480150,3.482490,3.484830,3.487170,3.489520,3.491870,3.494220,3.496580,3.498940,3.501300,3.503660,3.506030,3.508410,3.510780,3.513160,3.515540,3.517930,3.520310,3.522710,3.525100,3.527500,3.529900,3.532300,3.534710,3.537120,3.539540,3.541960,3.544380,3.546800,3.549230,3.551660,3.554090,3.556530,3.558970,3.561420,3.563860,3.566320,3.568770,3.571230,3.573690,3.576150,3.578620,3.581090,3.583570,3.586050,3.588530,3.591010,3.593500,3.595990,3.598490,3.600990,3.603490,3.606000,3.608510,3.611020,3.613540,3.616060,3.618580,3.621110,3.623640,3.626170,3.628710,3.631250,3.633790,3.636340,3.638890,3.641450,3.644010,3.646570,3.649140,3.651710,3.654280,3.656860,3.659440,3.662020,3.664610,3.667200,3.669800,3.672400,3.675000,3.677610,3.680220,3.682830,3.685450,3.688070,3.690690,3.693320,3.695950,3.698590,3.701230,3.703870,3.706520,3.709170,3.711830,3.714490,3.717150,3.719820,3.722490,3.725160,3.727840,3.730520,3.733210,3.735900,3.738590,3.741290,3.743990,3.746690,3.749400,3.752120,3.754830,3.757550,3.760280,3.763010,3.765740,3.768480,3.771220,3.773960,3.776710,3.779460,3.782220,3.784980,3.787750,3.790510,3.793290,3.796060,3.798840,3.801630,3.804420,3.807210,3.810010,3.812810,3.815620,3.818430,3.821240,3.824060,3.826880,3.829710,3.832540,3.835370,3.838210,3.841050,3.843900,3.846750,3.849610,3.852470,3.855330,3.858200,3.861070,3.863950,3.866830,3.869720,3.872610,3.875500,3.878400,3.881300,3.884210,3.887120,3.890040,3.892960,3.895890,3.898810,3.901750,3.904690,3.907630,3.910580,3.913530,3.916480,3.919440,3.922410,3.925380,3.928350,3.931330,3.934310,3.937300,3.940290,3.943290,3.946290,3.949290,3.952300,3.955320,3.958340,3.961360,3.964390,3.967420,3.970460,3.973500,3.976550,3.979600,3.982660,3.985720,3.988790,3.991860,3.994930,3.998010,4.001100,4.004190,4.007280,4.010380]

for i = 0, n_elements(l) -1 do print, l(i)*1E4,f(i) 

;ch2

f = [0.001201,0.001099,0.001117,0.001213,0.001211,0.001253,0.001231,0.001209,0.001237,0.001197,0.001128,0.001091,0.001037,0.000956,0.000858,0.000717,0.000602,0.000578,0.000587,0.000474,0.000377,0.000317,0.000233,0.000230,0.000157,0.000228,0.000424,0.000474,0.000563,0.000633,0.000451,0.000464,0.000850,0.000942,0.000812,0.001039,0.001136,0.001040,0.001322,0.001361,0.001194,0.001184,0.001108,0.001127,0.001191,0.001142,0.001138,0.000946,0.000972,0.001021,0.000773,0.000693,0.000647,0.000523,0.000451,0.000366,0.000398,0.000391,0.000223,0.000226,0.000307,0.000401,0.000664,0.001206,0.001624,0.001920,0.002697,0.003559,0.004236,0.005163,0.006273,0.007327,0.008786,0.010680,0.012370,0.014320,0.016680,0.019370,0.022320,0.025810,0.029940,0.034790,0.040240,0.046900,0.054810,0.063810,0.074360,0.086690,0.101000,0.117500,0.135800,0.157300,0.179300,0.202300,0.228100,0.252200,0.277500,0.300400,0.320700,0.339000,0.355000,0.366900,0.375500,0.380500,0.383700,0.385400,0.384600,0.382400,0.379800,0.378100,0.374300,0.369500,0.365700,0.360700,0.356800,0.353100,0.348800,0.344300,0.340200,0.336700,0.333700,0.330800,0.328300,0.326800,0.326200,0.325400,0.324900,0.325100,0.326000,0.328300,0.330300,0.332400,0.335100,0.338600,0.344800,0.350700,0.355200,0.361700,0.368900,0.376400,0.383800,0.391600,0.399500,0.408500,0.418100,0.426900,0.434900,0.441100,0.447300,0.454500,0.460500,0.466000,0.470900,0.475500,0.479800,0.484100,0.488500,0.490800,0.492200,0.494800,0.496700,0.497200,0.497400,0.497000,0.497200,0.496300,0.494200,0.492300,0.490400,0.488300,0.487300,0.484700,0.480500,0.477300,0.474900,0.472800,0.470400,0.468000,0.467400,0.466400,0.464600,0.465100,0.465900,0.465700,0.466200,0.467700,0.471400,0.474400,0.477300,0.480700,0.483200,0.486400,0.490200,0.492900,0.495400,0.498300,0.501400,0.504600,0.506900,0.509000,0.510800,0.512800,0.514800,0.515600,0.516200,0.519400,0.522800,0.523700,0.525200,0.527900,0.529900,0.531300,0.533600,0.536800,0.539100,0.540700,0.542300,0.543800,0.545900,0.547300,0.547000,0.547100,0.547900,0.548000,0.548100,0.548100,0.546400,0.546200,0.546700,0.544300,0.541900,0.541700,0.541800,0.540000,0.537500,0.536500,0.535100,0.532300,0.531000,0.529700,0.526600,0.523500,0.521800,0.519100,0.516500,0.515500,0.513500,0.511000,0.509900,0.508200,0.506700,0.506900,0.506900,0.505200,0.504000,0.503700,0.503200,0.503400,0.503200,0.501800,0.500600,0.499800,0.498100,0.495500,0.493900,0.493000,0.491100,0.488900,0.486000,0.483600,0.481200,0.477700,0.475200,0.472500,0.469700,0.470500,0.471700,0.472600,0.474400,0.475100,0.476500,0.479700,0.482600,0.484500,0.486300,0.488900,0.491000,0.492300,0.494100,0.495500,0.497000,0.498300,0.498500,0.497900,0.496600,0.494600,0.493800,0.492300,0.489100,0.485300,0.483300,0.480500,0.476700,0.473800,0.470700,0.467100,0.464800,0.462200,0.458200,0.454300,0.451300,0.449600,0.448500,0.446700,0.445200,0.443600,0.442100,0.442100,0.442100,0.440500,0.439900,0.441000,0.441100,0.440500,0.440900,0.441000,0.440600,0.441100,0.441000,0.441100,0.441400,0.441800,0.443500,0.444100,0.444600,0.445400,0.446100,0.446000,0.443800,0.441900,0.439100,0.434400,0.427100,0.418700,0.410800,0.402300,0.393600,0.386100,0.380000,0.376400,0.374300,0.372600,0.371600,0.372800,0.368600,0.356200,0.330300,0.293200,0.247800,0.192700,0.139000,0.096840,0.065350,0.042360,0.027350,0.017490,0.010810,0.006694,0.003758,0.001993,0.000901,0.000591,0.000486,0.000303,0.000300,0.000369,0.000303,0.000205,0.000180,0.000096,0.000108,0.000118,0.000190,0.000274,0.000265,0.000305,0.000315,0.000283,0.000371,0.000458,0.000466,0.000460,0.000500,0.000489,0.000411,0.000421,0.000477,0.000431,0.000384,0.000402,0.000386]


l = [3.722490,3.725160,3.727840,3.730520,3.733210,3.735900,3.738590,3.741290,3.743990,3.746690,3.749400,3.752120,3.754830,3.757550,3.760280,3.763010,3.765740,3.768480,3.771220,3.773960,3.776710,3.779460,3.782220,3.784980,3.787750,3.790510,3.793290,3.796060,3.798840,3.801630,3.804420,3.807210,3.810010,3.812810,3.815620,3.818430,3.821240,3.824060,3.826880,3.829710,3.832540,3.835370,3.838210,3.841050,3.843900,3.846750,3.849610,3.852470,3.855330,3.858200,3.861070,3.863950,3.866830,3.869720,3.872610,3.875500,3.878400,3.881300,3.884210,3.887120,3.890040,3.892960,3.895890,3.898810,3.901750,3.904690,3.907630,3.910580,3.913530,3.916480,3.919440,3.922410,3.925380,3.928350,3.931330,3.934310,3.937300,3.940290,3.943290,3.946290,3.949290,3.952300,3.955320,3.958340,3.961360,3.964390,3.967420,3.970460,3.973500,3.976550,3.979600,3.982660,3.985720,3.988790,3.991860,3.994930,3.998010,4.001100,4.004190,4.007280,4.010380,4.013490,4.016590,4.019710,4.022830,4.025950,4.029080,4.032210,4.035350,4.038490,4.041640,4.044790,4.047950,4.051110,4.054280,4.057450,4.060630,4.063810,4.067000,4.070190,4.073390,4.076590,4.079800,4.083010,4.086230,4.089450,4.092680,4.095910,4.099150,4.102390,4.105640,4.108890,4.112150,4.115420,4.118680,4.121960,4.125240,4.128520,4.131810,4.135110,4.138410,4.141710,4.145020,4.148340,4.151660,4.154990,4.158320,4.161660,4.165000,4.168350,4.171700,4.175060,4.178420,4.181790,4.185170,4.188550,4.191930,4.195330,4.198720,4.202130,4.205530,4.208950,4.212370,4.215790,4.219220,4.222660,4.226100,4.229550,4.233000,4.236460,4.239920,4.243390,4.246870,4.250350,4.253830,4.257330,4.260820,4.264330,4.267840,4.271350,4.274870,4.278400,4.281930,4.285470,4.289020,4.292570,4.296130,4.299690,4.303260,4.306830,4.310410,4.314000,4.317590,4.321190,4.324790,4.328400,4.332020,4.335640,4.339270,4.342900,4.346540,4.350190,4.353840,4.357500,4.361160,4.364830,4.368510,4.372190,4.375880,4.379580,4.383280,4.386990,4.390710,4.394430,4.398150,4.401890,4.405630,4.409370,4.413130,4.416890,4.420650,4.424420,4.428200,4.431990,4.435780,4.439570,4.443380,4.447190,4.451010,4.454830,4.458660,4.462500,4.466340,4.470190,4.474050,4.477910,4.481780,4.485660,4.489540,4.493430,4.497330,4.501240,4.505150,4.509060,4.512990,4.516920,4.520860,4.524800,4.528750,4.532710,4.536680,4.540650,4.544630,4.548620,4.552610,4.556610,4.560620,4.564630,4.568650,4.572680,4.576720,4.580760,4.584810,4.588870,4.592930,4.597010,4.601090,4.605170,4.609270,4.613370,4.617470,4.621590,4.625710,4.629840,4.633980,4.638120,4.642280,4.646440,4.650600,4.654780,4.658960,4.663150,4.667350,4.671550,4.675760,4.679980,4.684210,4.688450,4.692690,4.696940,4.701200,4.705470,4.709740,4.714020,4.718310,4.722610,4.726910,4.731230,4.735550,4.739870,4.744210,4.748560,4.752910,4.757270,4.761640,4.766010,4.770400,4.774790,4.779190,4.783600,4.788020,4.792440,4.796880,4.801320,4.805770,4.810230,4.814690,4.819170,4.823650,4.828140,4.832640,4.837150,4.841660,4.846190,4.850720,4.855270,4.859820,4.864370,4.868940,4.873520,4.878100,4.882700,4.887300,4.891910,4.896530,4.901160,4.905790,4.910440,4.915090,4.919760,4.924430,4.929110,4.933800,4.938500,4.943210,4.947920,4.952650,4.957380,4.962130,4.966880,4.971640,4.976410,4.981190,4.985980,4.990780,4.995590,5.000410,5.005230,5.010070,5.014920,5.019770,5.024630,5.029510,5.034390,5.039280,5.044190,5.049100,5.054020,5.058950,5.063890,5.068840,5.073800,5.078770,5.083750,5.088740,5.093740,5.098740,5.103760,5.108790,5.113830,5.118880,5.123940,5.129000,5.134080,5.139170,5.144270,5.149380,5.154500,5.159630,5.164760,5.169910,5.175070,5.180240,5.185420,5.190610,5.195810,5.201030,5.206250,5.211480,5.216720,5.221980]

for i = 0, n_elements(l) -1 do print, l(i)*1E4,f(i) 

end