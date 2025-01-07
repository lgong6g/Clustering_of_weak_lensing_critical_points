LaunchKernels[32];
l=ReadList["/home/barthele/Laurence/peak_correlation/2d_multipole_array_0_to_2e4_121bins.dat", Number, RecordLists->True]//Flatten;
lmin = l[[1]];
lmax = l[[-1]];

pks=ReadList["/home/barthele/Laurence/peak_correlation/zslice0.75_2d_kappa_power_spectrum_ell0_to_2e4_121bins.dat",Number,RecordLists->True]//Flatten;
Pk = {l, pks} // Transpose // Interpolation;
  
W[k_, r_] := Exp[-k^2*r^2/2];
smoothangle = 15/60/180*Pi;

list = ParallelTable[{{r, n, m}, NIntegrate[BesselJ[n, r*k]*k^m*W[k, smoothangle]^2*Pk[k], {k, lmin, lmax}, Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0, "MaxErrorIncreases" -> 2000}, MaxRecursion -> 40, WorkingPrecision -> 16, PrecisionGoal-> 5, AccuracyGoal -> Infinity]}, {r, 0, 1.0, .0001}, {n, 0, 4, 1}, {m, 1, 5, 1}];
Export["/home/barthele/Laurence/peak_correlation/X_table_r_nbessel_mk_planck18_cosmology_zslice0.75_r1.0_dr_0.0001_ell0_to_2e4_15arcmin_smooth.m", list];
Exit[];
