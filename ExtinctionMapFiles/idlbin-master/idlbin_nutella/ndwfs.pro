pro ndwfs

readlargecol, '/Users/jkrick/nep/clusters/flamex/NDWFS_Bw_34_35.cat.txt', v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12,v13,v14,v15, $
            v16,v17,v18,v19,v20,v21,v22,v23,v24,v25, v26,v27,v28,v29,v30,v31,v32,$
                 v33,v34,v35,v36,v37,v38,v39,v40,v41,v42,v43,v44,v45,v46,v47,v48,$
                 v49,v50,v51,v52,v53,v54,v55,v56,v57,v58,v59,v60,v61,v62,v63,v64,$
                 v65,v66,v67,v68,v69,v70,v71,v72,v73,v74,v75,v76,v77,v78,v79,v80,$
                 v81,v82,v83,v84,v85,v86,v87,v88,v89,v90,v91,v92,v93,v94,v95,v96,$
                 v97,v98,v99,v100,v101, format='A'


a = where(v1 gt 0. and v68 lt 35)
print, n_elements(a)
openw, outlun, '/Users/jkrick/nep/clusters/flamex/NDWFS_Bw.cat.txt', /get_lun
for i = 0, n_elements(a) - 1 do begin
   printf, outlun, v1(a),v2(a),v3(a),v4(a),v5(a),v6(a),v7(a),v8(a),v9(a),v10(a),v11(a),v12(a),v13(a),v14(a),v15(a), $
            v16(a),v17(a),v18(a),v19(a),v20(a),v21(a),v22(a),v23(a),v24(a),v25(a), v26(a),v27(a),v28(a),v29(a),v30(a),v31(a),v32(a),$
                 v33(a),v34(a),v35(a),v36(a),v37(a),v38(a),v39(a),v40(a),v41(a),v42(a),v43(a),v44(a),v45(a),v46(a),v47(a),v48(a),$
                 v49(a),v50(a),v51(a),v52(a),v53(a),v54(a),v55(a),v56(a),v57(a),v58(a),v59(a),v60(a),v61(a),v62(a),v63(a),v64(a),$
                 v65(a),v66(a),v67(a),v68(a),v69(a),v70(a),v71(a),v72(a),v73(a),v74(a),v75(a),v76(a),v77(a),v78(a),v79(a),v80(a),$
                 v81(a),v82(a),v83(a),v84(a),v85(a),v86(a),v87(a),v88(a),v89(a),v90(a),v91(a),v92(a),v93(a),v94(a),v95(a),v96(a),$
                 v97(a),v98(a),v99(a),v100(a), v101(a)
endfor
free_lun, outlun
end
