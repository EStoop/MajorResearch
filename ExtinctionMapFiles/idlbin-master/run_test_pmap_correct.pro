pro run_test_pmap_correct
;ganz normal
;test_pmap_correct, 'HD7924b', /run_nn, /plot_pmap, /plot_nn; /run_pmap,

;without the first four frames in the pmap dataset
;test_pmap_correct, 'HD7924b',  /run_pmap,  /plot_pmap, /pmapff

;with half weights on pixphasecorr np & xyfwhm and all 64 frames
;test_pmap_correct, 'HD7924b', 0, 0, /run_nn,  /plot_nn, /halfweights
;----------------


;ganz normal
test_pmap_correct, 'HAT-P-22', 1, 1, /run_pmap, /run_nn, /plot_pmap, /plot_nn

;without the first four frames in the pmap dataset
test_pmap_correct, 'HAT-P-22',  1, 1, /run_pmap,  /plot_pmap,/pmapff

;with half weights on pixphasecorr np & xyfwhm and all 64 frames
test_pmap_correct, 'HAT-P-22', 1, 1, /run_nn,  /plot_nn, /halfweights

end
