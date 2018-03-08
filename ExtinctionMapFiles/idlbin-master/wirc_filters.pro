pro wirc_filters
!P.multi = [0,1,2]

;k_cont


	
l = [	2159.121557,	2160.119698,	2161.11784 ,	2162.115982,	2163.114124,	2164.112265,	2165.110407,	2166.108549,	2167.106691,	2168.104832,	2169.102974,	2170.101116,	2171.099258,	2172.097399,	2173.095541,	2174.093683,	2175.091825,	2176.089966,	2177.088108,	2178.08625 ,	2179.084392,	2180.082533,	2181.080675,	2182.078817,	2183.076959,	2184.075101,	2185.073242,	2186.071384,	2187.069526,	2188.067668,	2189.065809,	2190.063951,	2191.062093,	2192.060235,	2193.058376,	2194.056518,	2195.05466 ,	2196.052802,	2197.050943,	2198.049085,	2199.047227,	2200.045369,	2201.04351 ,	2202.041652,	2203.039794,	2204.037936,	2205.036077,	2206.034219,	2207.032361,	2208.030503,	2209.028644,	2210.026786,	2211.024928,	2212.02307 ,	2213.021211,	2214.019353,	2215.017495,	2216.015637,	2217.013778,	2218.01192 ,	2219.010062,	2220.008204,	2221.006345,	2222.004487,	2223.002629,	2224.000771,	2224.998913,	2225.997054,	2226.995196,	2227.993338,	2228.99148 ,	2229.989621,	2230.987763,	2231.985905,	2232.984047,	2233.982188,	2234.98033 ,	2235.978472,	2236.976614,	2237.974755,	2238.972897,	2239.971039,	2240.969181,	2241.967322,	2242.965464,	2243.963606,	2244.961748,	2245.959889,	2246.958031,	2247.956173,	2248.954315,	2249.952456,	2250.950598,	2251.94874 ,	2252.946882,	2253.945023,	2254.943165,	2255.941307,	2256.939449,	2257.93759 ,	2258.935732,	2259.933874,	2260.932016,	2261.930157,	2262.928299,	2263.926441,	2264.924583,	2265.922725,	2266.920866,	2267.919008,	2268.91715 ,	2269.915292,	2270.913433,	2271.911575,	2272.909717,	2273.907859,	2274.906   ,	2275.904142,	2276.902284,	2277.900426,	2278.898567,	2279.896709,	2280.894851,	2281.892993,	2282.891134,	2283.889276,	2284.887418,	2285.88556 ,	2286.883701,	2287.881843,	2288.879985,	2289.878127,	2290.876268,	2291.87441 ,	2292.872552,	2293.870694,	2294.868835,	2295.866977,	2296.865119,	2297.863261,	2298.861402,	2299.859544,	2300.857686,	2301.855828,	2302.85397 ,	2303.852111,	2304.850253,	2305.848395,	2306.846537,	2307.844678,	2308.84282 ,	2309.840962,	2310.839104,	2311.837245,	2312.835387,	2313.833529,	2314.831671,	2315.829812,	2316.827954,	2317.826096,	2318.824238,	2319.822379,	2320.820521,	2321.818663,	2322.816805,	2323.814946,	2324.813088,	2325.81123 ,	2326.809372,	2327.807513,	2328.805655,	2329.803797,	2330.801939,	2331.80008 ,	2332.798222,	2333.796364,	2334.794506,	2335.792647,	2336.790789,	2337.788931,	2338.787073,	2339.785214,	2340.783356,	2341.781498,	2342.77964 ,	2343.777782,	2344.775923,	2345.774065,	2346.772207,	2347.770349,	2348.76849 ,	2349.766632,	2350.764774,	2351.762916,	2352.761057,	2353.759199,	2354.757341,	2355.755483,	2356.753624,	2357.751766,	2358.749908]

t = [0.0007,0.00025,0,0,0.00005,0.00035,0.0006,0.0003,0.00075,0,0,0.0003,0.00005,0.0003,0.00025,0.00075,0.0002,0.00025,0,0.0004,0.0007,0.0005,0,0.0003,0.00035,0.00025,0.00055,0.00055,0,0,0.00045,0.00045,0,0,0.0002,0.00005,0.0003,0.00085,0.00045,0,0.00015,0.0003,0,0.0004,0.001,0.0006,0,0.0001,0.00005,0.0002,0.00055,0.0003,0.00075,0.0005,0.0004,0.0005,0.00115,0.00035,0.00095,0.00035,0.001,0.00105,0.00125,0.00195,0.0018,0.00175,0.00235,0.00245,0.00265,0.00375,0.0043,0.00545,0.0063,0.0076,0.009,0.0113,0.01395,0.0168,0.0211,0.0255,0.0318,0.0399,0.0504,0.06275,0.07855,0.09795,0.1229,0.1537,0.1854,0.2236,0.2637,0.30585,0.34965,0.3932,0.43405,0.47455,0.50975,0.54645,0.5784,0.6105,0.6396,0.66565,0.6851,0.699,0.7067,0.70645,0.6984,0.68655,0.67065,0.6531,0.63695,0.6219,0.61065,0.60415,0.6042,0.6089,0.62015,0.63565,0.6548,0.67195,0.6794,0.6655,0.6251,0.55715,0.4672,0.3729,0.28735,0.2159,0.15865,0.12015,0.08955,0.06845,0.0514,0.04045,0.03215,0.0256,0.01995,0.01635,0.01325,0.0107,0.00925,0.0074,0.00625,0.00525,0.00475,0.00425,0.00435,0.00365,0.00275,0.0018,0.0023,0.00185,0.00155,0.0016,0.00155,0.00125,0.001,0.0009,0.00095,0.00115,0.00035,0.0006,0.0003,0.00065,0.00035,0.0006,0.00045,0.0003,0.0003,0.0006,0.00045,0.0005,0.00025,0,0.0007,0.00025,0.0006,0.00015,0.0003,0,0.00045,0.00035,0.0004,0.00045,0.00045,0.00025,0,0.00015,0.00015,0,0,0,0,0,0,0,0,0,0.0001,0.00045,0.0001]

plot, l, t, title = "kcont", xrange = [2200, 2350], charsize = 1.5


;co
l = [	2155.238985	,	2156.237127	,	2157.235268	,	2158.23341	,	2159.231552	,	2160.229694	,	2161.227835	,	2162.225977	,	2163.224119	,	2164.222261	,	2165.220402	,	2166.218544	,	2167.216686	,	2168.214828	,	2169.212969	,	2170.211111	,	2171.209253	,	2172.207395	,	2173.205536	,	2174.203678	,	2175.20182	,	2176.199962	,	2177.198103	,	2178.196245	,	2179.194387	,	2180.192529	,	2181.19067	,	2182.188812	,	2183.186954	,	2184.185096	,	2185.183237	,	2186.181379	,	2187.179521	,	2188.177663	,	2189.175805	,	2190.173946	,	2191.172088	,	2192.17023	,	2193.168372	,	2194.166513	,	2195.164655	,	2196.162797	,	2197.160939	,	2198.15908	,	2199.157222	,	2200.155364	,	2201.153506	,	2202.151647	,	2203.149789	,	2204.147931	,	2205.146073	,	2206.144214	,	2207.142356	,	2208.140498	,	2209.13864	,	2210.136781	,	2211.134923	,	2212.133065	,	2213.131207	,	2214.129348	,	2215.12749	,	2216.125632	,	2217.123774	,	2218.121915	,	2219.120057	,	2220.118199	,	2221.116341	,	2222.114482	,	2223.112624	,	2224.110766	,	2225.108908	,	2226.107049	,	2227.105191	,	2228.103333	,	2229.101475	,	2230.099617	,	2231.097758	,	2232.0959	,	2233.094042	,	2234.092184	,	2235.090325	,	2236.088467	,	2237.086609	,	2238.084751	,	2239.082892	,	2240.081034	,	2241.079176	,	2242.077318	,	2243.075459	,	2244.073601	,	2245.071743	,	2246.069885	,	2247.068026	,	2248.066168	,	2249.06431	,	2250.062452	,	2251.060593	,	2252.058735	,	2253.056877	,	2254.055019	,	2255.05316	,	2256.051302	,	2257.049444	,	2258.047586	,	2259.045727	,	2260.043869	,	2261.042011	,	2262.040153	,	2263.038294	,	2264.036436	,	2265.034578	,	2266.03272	,	2267.030862	,	2268.029003	,	2269.027145	,	2270.025287	,	2271.023429	,	2272.02157	,	2273.019712	,	2274.017854	,	2275.015996	,	2276.014137	,	2277.012279	,	2278.010421	,	2279.008563	,	2280.006704	,	2281.004846	,	2282.002988	,	2283.00113	,	2283.999271	,	2284.997413	,	2285.995555	,	2286.993697	,	2287.991838	,	2288.98998	,	2289.988122	,	2290.986264	,	2291.984405	,	2292.982547	,	2293.980689	,	2294.978831	,	2295.976972	,	2296.975114	,	2297.973256	,	2298.971398	,	2299.969539	,	2300.967681	,	2301.965823	,	2302.963965	,	2303.962106	,	2304.960248	,	2305.95839	,	2306.956532	,	2307.954674	,	2308.952815	,	2309.950957	,	2310.949099	,	2311.947241	,	2312.945382	,	2313.943524	,	2314.941666	,	2315.939808	,	2316.937949	,	2317.936091	,	2318.934233	,	2319.932375	,	2320.930516	,	2321.928658	,	2322.9268	,	2323.924942	,	2324.923083	,	2325.921225	,	2326.919367	,	2327.917509	,	2328.91565	,	2329.913792	,	2330.911934	,	2331.910076	,	2332.908217	,	2333.906359	,	2334.904501	,	2335.902643	,	2336.900784	,	2337.898926	,	2338.897068	,	2339.89521	,	2340.893351	,	2341.891493	,	2342.889635	,	2343.887777	,	2344.885918	,	2345.88406	,	2346.882202	,	2347.880344	,	2348.878486	,	2349.876627	,	2350.874769	,	2351.872911	,	2352.871053	,	2353.869194	,	2354.867336	]


t = [0.0002,0,0.0002,0.0002,0.00065,0.0004,0.00015,0.00015,0,0,0.00005,0.00035,0,0.00005,0.00005,0.00015,0.00035,0,0,0.0001,0.00045,0,0.00035,0.0004,0,0.0004,0.00025,0.0002,0,0.00015,0.00005,0,0,0.0002,0.0002,0,0.0003,0,0,0,0,0,0,0.0003,0,0.00055,0.0007,0.00115,0.00045,0,0,0,0.0002,0,0,0,0.0001,0.00005,0.0001,0.0004,0.0005,0.0003,0,0,0,0.00005,0.0001,0,0.0002,0,0.00005,0,0.00035,0,0,0.0005,0.00055,0.0007,0.0001,0.00045,0,0.0001,0.0002,0.0004,0.0007,0,0.0002,0.0005,0.001,0.00085,0.00105,0.00105,0.0011,0.0016,0.00145,0.0024,0.00255,0.00305,0.0034,0.0044,0.00535,0.0062,0.00725,0.009,0.0115,0.01345,0.0169,0.0205,0.02515,0.0313,0.03875,0.0488,0.05995,0.07425,0.08965,0.1085,0.12975,0.1544,0.1816,0.21015,0.24235,0.2766,0.3105,0.347,0.3849,0.42455,0.46505,0.50435,0.5397,0.5727,0.59605,0.6122,0.6167,0.61285,0.5999,0.581,0.5596,0.53675,0.5164,0.49995,0.4896,0.48245,0.4814,0.48645,0.49735,0.5135,0.5328,0.54735,0.5504,0.5348,0.49245,0.42595,0.34605,0.2663,0.19705,0.14525,0.1057,0.07775,0.05775,0.04315,0.0333,0.02525,0.01935,0.01585,0.0118,0.01015,0.00805,0.0066,0.0056,0.0048,0.00435,0.00315,0.003,0.00225,0.00225,0.0013,0.0013,0.001,0.0008,0.0009,0.00095,0.00065,0.00065,0.0013,0.00065,0.0002,0.0001,0.00015,0.0005,0.0003,0.001,0.00055,0.00055,0.00035,0.00055,0.0001,0,0.0001,0.0001,0.0002,0.00045]

plot, l, t, title = 'CO', xrange = [2200, 2350], charsize = 1.5

end