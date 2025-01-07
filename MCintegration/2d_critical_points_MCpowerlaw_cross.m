LaunchKernels[16];

W[k_, r_] := Exp[-k^2*r^2/2];
R = 15/60/180*Pi;

If[Length[$CommandLine] > 3,
    arg = $CommandLine[[4]];
    i = Quiet[ToExpression[arg], {ToExpression::sntx}],
    (
        Print["Error: No input argument provided."];
        Exit[];
    )
];

If[!NumericQ[i],
    Print["Error: Invalid input argument. Must be numeric."];
    Exit[];
];

m = -0.95047263101408`;
sigma2[n_] := NIntegrate[(1/2/Pi)*ell^(2*n + m + 1)*W[ell, R]^2, {ell, 0, Infinity}];
sigma20 = sigma2[0];
sigma21 = sigma2[1];
sigma22 = sigma2[2];
gamma = sigma21/Sqrt[sigma20*sigma22];

c11[r_] := (1/sigma20)*(1/2/Pi)*NIntegrate[k*BesselJ[0, k*r]*W[k, R]^2*k^m, {k, 0, Infinity}, MaxRecursion -> 200, WorkingPrecision -> 32, PrecisionGoal -> 10, AccuracyGoal -> Infinity];
c12[r_] := (-1/Sqrt[sigma20*sigma21])*(1/2/Pi)*NIntegrate[k^2*BesselJ[1, k*r]*W[k, R]^2*k^m, {k, 0, Infinity}, MaxRecursion -> 200, WorkingPrecision -> 32, PrecisionGoal -> 10, AccuracyGoal -> Infinity];
c14[r_] := (-1/Sqrt[sigma20*sigma22])*(1/2/Pi)*NIntegrate[k^3/2*(Hypergeometric0F1Regularized[2, (-1/4)* k^2* r^2] - 2*BesselJ[2, k*r])*W[k, R]^2*k^m, {k, 0, Infinity}, MaxRecursion -> 200, WorkingPrecision -> 32, PrecisionGoal -> 10, AccuracyGoal -> Infinity];
c16[r_] := (-1/Sqrt[sigma20*sigma22])*(1/2/Pi)*NIntegrate[k^3/2*Hypergeometric0F1Regularized[2, (-1/4)* k^2* r^2]*W[k, R]^2*k^m, {k, 0, Infinity}, MaxRecursion -> 200, WorkingPrecision -> 32, PrecisionGoal -> 10, AccuracyGoal -> Infinity];
c21[r_] := -c12[r];
c22[r_] := (1/sigma21)*(1/2/Pi)*NIntegrate[k^3/2*(Hypergeometric0F1Regularized[2, (-1/4)* k^2* r^2] - 2*BesselJ[2, k*r])*W[k, R]^2*k^m, {k, 0, Infinity}, MaxRecursion -> 200, WorkingPrecision -> 32, PrecisionGoal -> 10, AccuracyGoal -> Infinity];
c24[r_] := (-1/Sqrt[sigma21*sigma22])*(1/2/Pi)*NIntegrate[k^4*(3*BesselJ[2, k*r]/(k*r) - BesselJ[3, k*r])*W[k, R]^2*k^m, {k, 0, Infinity}, MaxRecursion -> 200, WorkingPrecision -> 32, PrecisionGoal -> 10, AccuracyGoal -> Infinity];
c26[r_] := (-1/Sqrt[sigma21*sigma22])*(1/2/Pi)*NIntegrate[k^4*(BesselJ[2, k*r]/(k*r))*W[k, R]^2*k^m, {k, 0, Infinity}, MaxRecursion -> 200, WorkingPrecision -> 32, PrecisionGoal -> 10, AccuracyGoal -> Infinity];
c33[r_] := (1/sigma21)*(1/2/Pi)*NIntegrate[k^3/2*Hypergeometric0F1Regularized[2, (-1/4)* k^2* r^2]*W[k, R]^2*k^m, {k, 0, Infinity}, MaxRecursion -> 200, WorkingPrecision -> 32, PrecisionGoal -> 10, AccuracyGoal -> Infinity];
c35[r_] := (-1/Sqrt[sigma21*sigma22])*(1/2/Pi)*NIntegrate[k^4*(BesselJ[2, k*r]/(k*r))*W[k, R]^2*k^m, {k, 0, Infinity}, MaxRecursion -> 200, WorkingPrecision -> 32, PrecisionGoal -> 10, AccuracyGoal -> Infinity];
c41[r_] := c14[r];
c42[r_] := -c24[r];
c44[r_] := (1/sigma22)*(1/2/Pi)*NIntegrate[k^5*(3/(k*r)^2 - 1)*BesselJ[2, k*r]*W[k, R]^2*k^m, {k, 0, Infinity}, MaxRecursion -> 200, WorkingPrecision -> 32, PrecisionGoal -> 10, AccuracyGoal -> Infinity];
c46[r_] := (1/sigma22)*(1/2/Pi)*NIntegrate[k^5/(k*r)*(BesselJ[1, k*r] - 3*BesselJ[2, k*r]/(k*r))*W[k, R]^2*k^m, {k, 0, Infinity}, MaxRecursion -> 200, WorkingPrecision -> 32, PrecisionGoal -> 10, AccuracyGoal -> Infinity];
c53[r_] := -c35[r];
c55[r_] := c46[r];
c61[r_] := c16[r];
c62[r_] := -c26[r];
c64[r_] := c46[r];
c66[r_] := (1/sigma22)*(1/2/Pi)*NIntegrate[k^5/(k*r)^2*3*BesselJ[2, k*r]*W[k, R]^2*k^m, {k, 0, Infinity}, MaxRecursion -> 200, WorkingPrecision -> 32, PrecisionGoal -> 10, AccuracyGoal -> Infinity];

cov2pt[r_] := {{1, 0, 0, -gamma/2, 0, -gamma/2, c11[r], c12[r], 0, c14[r], 0, c16[r]},{0, 1/2, 0, 0, 0, 0, c21[r], c22[r], 0, c24[r], 0, c26[r]},{0, 0, 1/2, 0, 0, 0, 0, 0, c33[r], 0, c35[r], 0},{-gamma/2, 0, 0, 3/8, 0, 1/8, c41[r], c42[r], 0, c44[r], 0, c46[r]},{0, 0, 0, 0, 1/8, 0, 0, 0, c53[r], 0, c55[r], 0},{-gamma/2, 0, 0, 1/8, 0, 3/8, c61[r], c62[r], 0, c64[r], 0, c66[r]},{c11[r], c21[r], 0, c41[r], 0, c61[r], 1, 0, 0, -gamma/2, 0, -gamma/2},{c12[r], c22[r], 0, c42[r], 0, c62[r], 0, 1/2, 0, 0, 0, 0},{0, 0, c33[r], 0, c53[r], 0, 0, 0, 1/2, 0, 0, 0},{c14[r], c24[r], 0, c44[r], 0, c64[r], -gamma/2, 0, 0, 3/8, 0, 1/8},{0, 0, c35[r], 0, c55[r], 0, 0, 0, 0, 0, 1/8, 0},{c16[r], c26[r], 0, c46[r], 0, c66[r], -gamma/2, 0, 0, 1/8, 0, 3/8}};
covauto = {{1, 0, 0, -gamma/2, 0, -gamma/2},{0, 1/2, 0, 0, 0, 0},{0, 0, 1/2, 0, 0, 0},{-gamma/2, 0, 0, 3/8, 0, 1/8},{0, 0, 0, 0, 1/8, 0},{-gamma/2, 0, 0, 1/8, 0, 3/8}};

ConditionnalMultinormalDistribution[pdf_, val_, idx_] := Module[{S = pdf[[2]], m = pdf[[1]], odx, Sigmaa, Sigmab, Sigmac,Mu2, S2, idx2, val2}, odx = Flatten[{Complement[Range[Length[S]], Flatten[{idx}]]}]; Sigmaa = (S[[odx]] // Transpose)[[odx]];idx2 = Flatten[{idx}]; val2 = Flatten[{val}];Sigmac = (S[[odx]] // Transpose)[[idx2]] // Transpose; Sigmab = (S[[idx2]] // Transpose)[[idx2]]; Mu2 = m[[odx]] + Sigmac . Inverse[Sigmab] . (val2 - m[[idx]]); S2 =Sigmaa -  Sigmac . Inverse[Sigmab] . Transpose[Sigmac]; S2 = 1/2 (S2 + Transpose[S2]); If[Length[Mu2] == 1, NormalDistribution[Mu2 // First, Sqrt[S2 // First // First]], MultinormalDistribution[Mu2, S2]]] /; Head[pdf] == MultinormalDistribution;

jointMCcum[rr_, nunu_, npt_] := (Cov = cov2pt[rr]; pdf = MultinormalDistribution[Table[0, {Length[Cov]}], Cov]; pdfcond =ConditionnalMultinormalDistribution[pdf, {0, 0, 0, 0}, {2, 3, 8, 9}]; Clear[pcrit]; pcrit = PDF[MarginalDistribution[pdf, {2, 3, 8, 9}], {0, 0, 0, 0}]; data = RandomVariate[pdfcond, npt]; Clear[ff]; ff[x_, x11_, x12_, x22_, y_, y11_, y12_, y22_] := Boole[x11 + x22 > 0 && x11 x22 - x12^2 > 0 && y11 y22 - y12^2 < 0 && x > nunu && y > nunu]*(x11 x22 - x12^2)*Abs[y11 y22 - y12^2]; (Plus@@Map[ff@@#&, data])/npt*pcrit);

peakMCcum1[nunu_, npt_] := (Cov = covauto; pdf = MultinormalDistribution[Table[0, {Length[Cov]}], Cov]; pdfcond = ConditionnalMultinormalDistribution[pdf, {0, 0}, {2, 3}]; Clear[pcrit]; pcrit = PDF[MarginalDistribution[pdf, {2, 3}], {0, 0}]; data = RandomVariate[pdfcond, npt]; Clear[ff]; ff[x_, x11_, x12_, x22_] := Boole[x11 + x22 > 0 && x11 x22 - x12^2 > 0 && x > nunu]*(x11 x22 - x12^2); (Plus@@Map[ff@@#&, data])/npt*pcrit);

peakMCcum2[nunu_, npt_] := (Cov = covauto; pdf = MultinormalDistribution[Table[0, {Length[Cov]}], Cov]; pdfcond = ConditionnalMultinormalDistribution[pdf, {0, 0}, {2, 3}]; Clear[pcrit]; pcrit = PDF[MarginalDistribution[pdf, {2, 3}], {0, 0}]; data = RandomVariate[pdfcond, npt]; Clear[ff]; ff[x_, x11_, x12_, x22_] := Boole[x11 x22 - x12^2 < 0 && x > nunu]*Abs[x11 x22 - x12^2]; (Plus@@Map[ff@@#&, data])/npt*pcrit);

np = peakMCcum1[0.3, 1000000]*peakMCcum2[0.3, 1000000];
tt = ParallelTable[{r, jointMCcum[r, 0.3, 20000000]/np}, {r, 0.000290888, 0.101811, 0.00145444}];

Export[StringJoin["/home/barthele/Laurence/peak_correlation/MCintegration/MC_integration_void_saddle_powerlaw", ToString[i], ".m"], tt];

Exit[];