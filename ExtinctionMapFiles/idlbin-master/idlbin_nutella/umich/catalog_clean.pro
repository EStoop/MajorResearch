IDL Version 6.2, Mac OS X (darwin ppc m32). (c) 2005, Research Systems, Inc.
Installation number: 91593-1.
Licensed for use by: Jet Propulsion Laboratory

IDL> .run "/Users/jkrick/idlbin/umich/histogram.pro"

colors = GetColor(/load, Start=1)
                   ^
% Syntax error.
  At: /Users/jkrick/idlbin/umich/histogram.pro, Line 3

plothist, ,magarr, xhist, yhist, bin=0.10;, /noprint,xrange=[0,480]
           ^
% Syntax error.
  At: /Users/jkrick/idlbin/umich/histogram.pro, Line 28
% 2 Compilation error(s) in module HISTOGRAM.
IDL> .comp GetColor
% Compiled module: COLOR24.
% Compiled module: GETCOLOR.
IDL> .run "/Users/jkrick/idlbin/umich/histogram.pro"

plothist, ,magarr, xhist, yhist, bin=0.10;, /noprint,xrange=[0,480]
           ^
% Syntax error.
  At: /Users/jkrick/idlbin/umich/histogram.pro, Line 28
% 1 Compilation error(s) in module HISTOGRAM.
IDL> .run "/Users/jkrick/idlbin/umich/histogram.pro"
% Compiled module: HISTOGRAM.
IDL> histogram
% CLOSE: Variable is undefined: LUN.
% Execution halted at: HISTOGRAM          21
  /Users/jkrick/idlbin/umich/histogram.pro
%                      $MAIN$          
IDL> .run "/Users/jkrick/idlbin/umich/histogram.pro"
% Procedure was compiled while active: HISTOGRAM. Returning.
% Compiled module: HISTOGRAM.
IDL> histogram
% Compiled module: PS_OPEN.
% DEVICE: Error opening file. Unit: 100
          File: Users/jkrick/palomar/lfc/catalog/rhist.ps
  No such file or directory
% Execution halted at: PS_OPEN           114
  /Users/jkrick/idlbin/umich/others/ps_open.pro
%                      HISTOGRAM          25
  /Users/jkrick/idlbin/umich/histogram.pro
%                      $MAIN$          
IDL> .run "/Users/jkrick/idlbin/umich/histogram.pro"
% Procedure was compiled while active: HISTOGRAM. Returning.
% Compiled module: HISTOGRAM.
IDL> histogram
% PS_OPEN: WARNING: device already opened to PS, closing it first.
% PS_OPEN:          any plot in progress will be lost.
% PS_OPEN: output redirect to PostScript file
           /Users/jkrick/palomar/LFC/catalog/rhist.ps
% Compiled module: PLOTHIST.
% Compiled module: MPFITFUN.
% Compiled module: MPFIT.
% Compiled module: GAUSS.
Iter      1   CHI-SQUARE =       17513346.          DOF = 206
    P(0) =              132.000  
    P(1) =              10.0000  
    P(2) =              5000.00  
Iter      1   CHI-SQUARE =       17513346.          DOF = 206
    P(0) =              132.000  
    P(1) =              10.0000  
    P(2) =              5000.00  
% Compiled module: PS_CLOSE.
% PS_CLOSE: plotting device restored as X
% Program caused arithmetic error: Floating underflow
% Program caused arithmetic error: Floating overflow
% Program caused arithmetic error: Floating illegal operand
IDL> .run "/Users/jkrick/idlbin/umich/histogram.pro"
% Compiled module: HISTOGRAM.
IDL> histogram
% PS_OPEN: output redirect to PostScript file
           /Users/jkrick/palomar/LFC/catalog/rhist.ps
Iter      1   CHI-SQUARE =   2.4751245E+08          DOF = 206
    P(0) =              25.0000  
    P(1) =              5.00000  
    P(2) =              5000.00  
Iter      2   CHI-SQUARE =       6769086.5          DOF = 206
    P(0) =              24.7988  
    P(1) =              4.17396  
    P(2) =              1224.91  
Iter      3   CHI-SQUARE =       2473440.5          DOF = 206
    P(0) =              24.5022  
    P(1) =              1.52605  
    P(2) =              1363.91  
Iter      4   CHI-SQUARE =       591341.88          DOF = 206
    P(0) =              25.0075  
    P(1) =              1.69507  
    P(2) =              1887.82  
Iter      5   CHI-SQUARE =       462679.53          DOF = 206
    P(0) =              24.8367  
    P(1) =              1.55158  
    P(2) =              1960.83  
Iter      6   CHI-SQUARE =       458281.97          DOF = 206
    P(0) =              24.8649  
    P(1) =              1.57100  
    P(2) =              1961.55  
Iter      7   CHI-SQUARE =       458224.22          DOF = 206
    P(0) =              24.8630  
    P(1) =              1.56710  
    P(2) =              1964.17  
Iter      8   CHI-SQUARE =       458223.31          DOF = 206
    P(0) =              24.8634  
    P(1) =              1.56744  
    P(2) =              1963.98  
Iter      8   CHI-SQUARE =       458223.31          DOF = 206
    P(0) =              24.8634  
    P(1) =              1.56744  
    P(2) =              1963.98  
% PS_CLOSE: plotting device restored as X
% Program caused arithmetic error: Floating underflow
IDL> 

  Process idl finished
IDL Version 6.2, Mac OS X (darwin ppc m32). (c) 2005, Research Systems, Inc.
Installation number: 91593-1.
Licensed for use by: Jet Propulsion Laboratory

IDL> .run "/Users/jkrick/idlbin/umich/histogram.pro"

colors = GetColor(/load, Start=1)
                   ^
% Syntax error.
  At: /Users/jkrick/idlbin/umich/histogram.pro, Line 3
% 1 Compilation error(s) in module HISTOGRAM.
IDL> .comp GetColor
% Compiled module: COLOR24.
% Compiled module: GETCOLOR.
IDL> .run "/Users/jkrick/idlbin/umich/histogram.pro"
% Compiled module: HISTOGRAM.
IDL> histogram
% Compiled module: PS_OPEN.
% PS_OPEN: output redirect to PostScript file
           /Users/jkrick/palomar/LFC/catalog/rhist.ps
% Compiled module: PLOTHIST.
% Compiled module: MPFITFUN.
% Compiled module: MPFIT.
% Compiled module: GAUSS.
Iter      1   CHI-SQUARE =   2.4751245E+08          DOF = 206
    P(0) =              25.0000  
    P(1) =              5.00000  
    P(2) =              5000.00  
Iter      2   CHI-SQUARE =       6769086.5          DOF = 206
    P(0) =              24.7988  
    P(1) =              4.17396  
    P(2) =              1224.91  
Iter      3   CHI-SQUARE =       2473440.5          DOF = 206
    P(0) =              24.5022  
    P(1) =              1.52605  
    P(2) =              1363.91  
Iter      4   CHI-SQUARE =       591341.88          DOF = 206
    P(0) =              25.0075  
    P(1) =              1.69507  
    P(2) =              1887.82  
Iter      5   CHI-SQUARE =       462679.53          DOF = 206
    P(0) =              24.8367  
    P(1) =              1.55158  
    P(2) =              1960.83  
Iter      6   CHI-SQUARE =       458281.97          DOF = 206
    P(0) =              24.8649  
    P(1) =              1.57100  
    P(2) =              1961.55  
Iter      7   CHI-SQUARE =       458224.22          DOF = 206
    P(0) =              24.8630  
    P(1) =              1.56710  
    P(2) =              1964.17  
Iter      8   CHI-SQUARE =       458223.31          DOF = 206
    P(0) =              24.8634  
    P(1) =              1.56744  
    P(2) =              1963.98  
Iter      8   CHI-SQUARE =       458223.31          DOF = 206
    P(0) =              24.8634  
    P(1) =              1.56744  
    P(2) =              1963.98  
% Compiled module: PS_CLOSE.
% PS_CLOSE: plotting device restored as X
% Program caused arithmetic error: Floating underflow
IDL> .run "/Users/jkrick/idlbin/umich/histogram.pro"
% Compiled module: HISTOGRAM.
IDL> histogram
% PS_OPEN: output redirect to PostScript file
           /Users/jkrick/palomar/LFC/catalog/rhist.ps
Iter      1   CHI-SQUARE =   2.4751245E+08          DOF = 206
    P(0) =              25.0000  
    P(1) =              5.00000  
    P(2) =              5000.00  
Iter      2   CHI-SQUARE =       6769086.5          DOF = 206
    P(0) =              24.7988  
    P(1) =              4.17396  
    P(2) =              1224.91  
Iter      3   CHI-SQUARE =       2473440.5          DOF = 206
    P(0) =              24.5022  
    P(1) =              1.52605  
    P(2) =              1363.91  
Iter      4   CHI-SQUARE =       591341.88          DOF = 206
    P(0) =              25.0075  
    P(1) =              1.69507  
    P(2) =              1887.82  
Iter      5   CHI-SQUARE =       462679.53          DOF = 206
    P(0) =              24.8367  
    P(1) =              1.55158  
    P(2) =              1960.83  
Iter      6   CHI-SQUARE =       458281.97          DOF = 206
    P(0) =              24.8649  
    P(1) =              1.57100  
    P(2) =              1961.55  
Iter      7   CHI-SQUARE =       458224.22          DOF = 206
    P(0) =              24.8630  
    P(1) =              1.56710  
    P(2) =              1964.17  
Iter      8   CHI-SQUARE =       458223.31          DOF = 206
    P(0) =              24.8634  
    P(1) =              1.56744  
    P(2) =              1963.98  
Iter      8   CHI-SQUARE =       458223.31          DOF = 206
    P(0) =              24.8634  
    P(1) =              1.56744  
    P(2) =              1963.98  
% PS_CLOSE: plotting device restored as X
% Program caused arithmetic error: Floating underflow
IDL> .run "/Users/jkrick/idlbin/umich/histogram.pro"
% Compiled module: HISTOGRAM.
IDL> histogram
% PS_OPEN: output redirect to PostScript file
           /Users/jkrick/palomar/LFC/catalog/rhist.ps
Iter      1   CHI-SQUARE =   2.4751245E+08          DOF = 206
    P(0) =              25.0000  
    P(1) =              5.00000  
    P(2) =              5000.00  
Iter      2   CHI-SQUARE =       6769086.5          DOF = 206
    P(0) =              24.7988  
    P(1) =              4.17396  
    P(2) =              1224.91  
Iter      3   CHI-SQUARE =       2473440.5          DOF = 206
    P(0) =              24.5022  
    P(1) =              1.52605  
    P(2) =              1363.91  
Iter      4   CHI-SQUARE =       591341.88          DOF = 206
    P(0) =              25.0075  
    P(1) =              1.69507  
    P(2) =              1887.82  
Iter      5   CHI-SQUARE =       462679.53          DOF = 206
    P(0) =              24.8367  
    P(1) =              1.55158  
    P(2) =              1960.83  
Iter      6   CHI-SQUARE =       458281.97          DOF = 206
    P(0) =              24.8649  
    P(1) =              1.57100  
    P(2) =              1961.55  
Iter      7   CHI-SQUARE =       458224.22          DOF = 206
    P(0) =              24.8630  
    P(1) =              1.56710  
    P(2) =              1964.17  
Iter      8   CHI-SQUARE =       458223.31          DOF = 206
    P(0) =              24.8634  
    P(1) =              1.56744  
    P(2) =              1963.98  
Iter      8   CHI-SQUARE =       458223.31          DOF = 206
    P(0) =              24.8634  
    P(1) =              1.56744  
    P(2) =              1963.98  
% PS_CLOSE: plotting device restored as X
% Program caused arithmetic error: Floating underflow
IDL> .run "/Users/jkrick/idlbin/umich/histogram.pro"
% Compiled module: HISTOGRAM.
IDL> histogram
% PS_OPEN: output redirect to PostScript file
           /Users/jkrick/palomar/LFC/catalog/rhist.ps
Iter      1   CHI-SQUARE =   2.4751245E+08          DOF = 206
    P(0) =              25.0000  
    P(1) =              5.00000  
    P(2) =              5000.00  
Iter      2   CHI-SQUARE =       6769086.5          DOF = 206
    P(0) =              24.7988  
    P(1) =              4.17396  
    P(2) =              1224.91  
Iter      3   CHI-SQUARE =       2473440.5          DOF = 206
    P(0) =              24.5022  
    P(1) =              1.52605  
    P(2) =              1363.91  
Iter      4   CHI-SQUARE =       591341.88          DOF = 206
    P(0) =              25.0075  
    P(1) =              1.69507  
    P(2) =              1887.82  
Iter      5   CHI-SQUARE =       462679.53          DOF = 206
    P(0) =              24.8367  
    P(1) =              1.55158  
    P(2) =              1960.83  
Iter      6   CHI-SQUARE =       458281.97          DOF = 206
    P(0) =              24.8649  
    P(1) =              1.57100  
    P(2) =              1961.55  
Iter      7   CHI-SQUARE =       458224.22          DOF = 206
    P(0) =              24.8630  
    P(1) =              1.56710  
    P(2) =              1964.17  
Iter      8   CHI-SQUARE =       458223.31          DOF = 206
    P(0) =              24.8634  
    P(1) =              1.56744  
    P(2) =              1963.98  
Iter      8   CHI-SQUARE =       458223.31          DOF = 206
    P(0) =              24.8634  
    P(1) =              1.56744  
    P(2) =              1963.98  
% PS_CLOSE: plotting device restored as X
% Program caused arithmetic error: Floating underflow
IDL> .run "/Users/jkrick/idlbin/catalog_merge.pro"
% Compiled module: CATALOG_MERGE.
IDL> catalog_merge
