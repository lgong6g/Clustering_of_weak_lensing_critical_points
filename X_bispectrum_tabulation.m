LaunchKernels[40];
l=ReadList["background_quantity/2d_multipole_array_0_to_2e4_121bins.dat", Number, RecordLists->True]//Flatten;
lmin = l[[1]];
lmax = l[[-1]];
zlospk=ReadList["background_quantity/z_los_0.05_until_1.5_100bins.dat",Number,RecordLists->True]//Flatten;
pklzgrid=Import["background_quantity/P_chi_z_grid_until_z0.05_to_1.5_100bins_l0_to_2e4_121bins_enhanced_BAO.dat"];
list={};
For[i=1,i<=Length[l],i++,nl={zlospk,pklzgrid[[i]]}//Transpose;
nl=Replace[nl,{x_,y_}:>{l[[i]],x,y},{1}];AppendTo[list,nl]];
pkgrid=list//Flatten[#,1]&;
Pkgrid=pkgrid//Interpolation;

W[k_, r_] := Exp[-k^2*r^2/2];
smoothangle = 15/60/180*Pi;

list1=ParallelTable[{{r,z},NIntegrate[BesselJ[2,r*k]*k*W[k,smoothangle]*Pkgrid[k,z],{k,lmin,lmax}, Method->{"GlobalAdaptive", "SymbolicProcessing"->0, "MaxErrorIncreases"-> 2000}, MaxRecursion->40, WorkingPrecision -> 16, PrecisionGoal->5, AccuracyGoal->Infinity]},{r,0.0,0.3,.0001},{z,0.05,1.5,0.01}];
Export["X_chi_21_table_r_nbessel_mk_planck18_cosmology_dirac_delta_z_r0.3_dr_0.0001_dz_0.01_ell0_to_2e4_15arcmin_smooth_enhanced_BAO.m", list1];

list2=ParallelTable[{{r,z},NIntegrate[BesselJ[2,r*k]*k^3*W[k,smoothangle]*Pkgrid[k,z],{k,lmin,lmax}, Method->{"GlobalAdaptive", "SymbolicProcessing"->0, "MaxErrorIncreases"-> 2000}, MaxRecursion->40, WorkingPrecision -> 16, PrecisionGoal->5, AccuracyGoal->Infinity]},{r,0.0,0.3,.0001},{z,0.05,1.5,0.01}];
Export["X_chi_23_table_r_nbessel_mk_planck18_cosmology_dirac_delta_z_r0.3_dr_0.0001_dz_0.01_ell0_to_2e4_15arcmin_smooth_enhanced_BAO.m", list2];

list3=ParallelTable[{{r,z},NIntegrate[BesselJ[3,r*k]*k^2*W[k,smoothangle]*Pkgrid[k,z],{k,lmin,lmax}, Method->{"GlobalAdaptive", "SymbolicProcessing"->0, "MaxErrorIncreases"-> 2000}, MaxRecursion->40, WorkingPrecision -> 16, PrecisionGoal->5, AccuracyGoal->Infinity]},{r,0.0,0.3,.0001},{z,0.05,1.5,0.01}];
Export["X_chi_32_table_r_nbessel_mk_planck18_cosmology_dirac_delta_z_r0.3_dr_0.0001_dz_0.01_ell0_to_2e4_15arcmin_smooth_enhanced_BAO.m", list3];

list4=ParallelTable[{{r,z},NIntegrate[BesselJ[3,r*k]*k^4*W[k,smoothangle]*Pkgrid[k,z],{k,lmin,lmax}, Method->{"GlobalAdaptive", "SymbolicProcessing"->0, "MaxErrorIncreases"-> 2000}, MaxRecursion->40, WorkingPrecision -> 16, PrecisionGoal->5, AccuracyGoal->Infinity]},{r,0.0,0.3,.0001},{z,0.05,1.5,0.01}];
Export["X_chi_34_table_r_nbessel_mk_planck18_cosmology_dirac_delta_z_r0.3_dr_0.0001_dz_0.01_ell0_to_2e4_15arcmin_smooth_enhanced_BAO.m", list4];

list5=ParallelTable[{{r,z},NIntegrate[BesselJ[4,r*k]*k^3*W[k,smoothangle]*Pkgrid[k,z],{k,lmin,lmax}, Method->{"GlobalAdaptive", "SymbolicProcessing"->0, "MaxErrorIncreases"-> 2000}, MaxRecursion->40, WorkingPrecision -> 16, PrecisionGoal->5, AccuracyGoal->Infinity]},{r,0.0,0.3,.0001},{z,0.05,1.5,0.01}];
Export["X_chi_43_table_r_nbessel_mk_planck18_cosmology_dirac_delta_z_r0.3_dr_0.0001_dz_0.01_ell0_to_2e4_15arcmin_smooth_enhanced_BAO.m", list5];

list6=ParallelTable[{{r,z},NIntegrate[BesselJ[0,r*k]*k^1*W[k,smoothangle]*Pkgrid[k,z],{k,lmin,lmax}, Method->{"GlobalAdaptive", "SymbolicProcessing"->0, "MaxErrorIncreases"-> 2000}, MaxRecursion->40, WorkingPrecision -> 16, PrecisionGoal->5, AccuracyGoal->Infinity]},{r,0.0,0.3,.0001},{z,0.05,1.5,0.01}];
Export["X_chi_01_table_r_nbessel_mk_planck18_cosmology_dirac_delta_z_r0.3_dr_0.0001_dz_0.01_ell0_to_2e4_15arcmin_smooth_enhanced_BAO.m", list6];

list7=ParallelTable[{{r,z},NIntegrate[BesselJ[0,r*k]*k^3*W[k,smoothangle]*Pkgrid[k,z],{k,lmin,lmax}, Method->{"GlobalAdaptive", "SymbolicProcessing"->0, "MaxErrorIncreases"-> 2000}, MaxRecursion->40, WorkingPrecision -> 16, PrecisionGoal->5, AccuracyGoal->Infinity]},{r,0.0,0.3,.0001},{z,0.05,1.5,0.01}];
Export["X_chi_03_table_r_nbessel_mk_planck18_cosmology_dirac_delta_z_r0.3_dr_0.0001_dz_0.01_ell0_to_2e4_15arcmin_smooth_enhanced_BAO.m", list7];

list8=ParallelTable[{{r,z},NIntegrate[BesselJ[1,r*k]*W[k,smoothangle]*Pkgrid[k,z],{k,lmin,lmax}, Method->{"GlobalAdaptive", "SymbolicProcessing"->0, "MaxErrorIncreases"-> 2000}, MaxRecursion->40, WorkingPrecision -> 16, PrecisionGoal->5, AccuracyGoal->Infinity]},{r,0.0,0.3,.0001},{z,0.05,1.5,0.01}];
Export["X_chi_10_table_r_nbessel_mk_planck18_cosmology_dirac_delta_z_r0.3_dr_0.0001_dz_0.01_ell0_to_2e4_15arcmin_smooth_enhanced_BAO.m", list8];

list9=ParallelTable[{{r,z},NIntegrate[BesselJ[1,r*k]*k^2*W[k,smoothangle]*Pkgrid[k,z],{k,lmin,lmax}, Method->{"GlobalAdaptive", "SymbolicProcessing"->0, "MaxErrorIncreases"-> 2000}, MaxRecursion->40, WorkingPrecision -> 16, PrecisionGoal->5, AccuracyGoal->Infinity]},{r,0.0,0.3,.0001},{z,0.05,1.5,0.01}];
Export["X_chi_12_table_r_nbessel_mk_planck18_cosmology_dirac_delta_z_r0.3_dr_0.0001_dz_0.01_ell0_to_2e4_15arcmin_smooth_enhanced_BAO.m", list9];

list10=ParallelTable[{{r,z},NIntegrate[BesselJ[1,r*k]*k^4*W[k,smoothangle]*Pkgrid[k,z],{k,lmin,lmax}, Method->{"GlobalAdaptive", "SymbolicProcessing"->0, "MaxErrorIncreases"-> 2000}, MaxRecursion->40, WorkingPrecision -> 16, PrecisionGoal->5, AccuracyGoal->Infinity]},{r,0.0,0.3,.0001},{z,0.05,1.5,0.01}];
Export["X_chi_14_table_r_nbessel_mk_planck18_cosmology_dirac_delta_z_r0.3_dr_0.0001_dz_0.01_ell0_to_2e4_15arcmin_smooth_enhanced_BAO.m", list10];

Exit[];
