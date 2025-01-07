LaunchKernels[32];
Om0 = 0.3143878052130123;
X0 = 1
X1 = -1;

g10000 = Import["renormalized_bias_functions_pvs/gp10000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.m"];
g10000 = g10000//Interpolation;
g1000 = g10000[0.3];

g01000 = Import["renormalized_bias_functions_pvs/gp01000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.m"];
g01000 = g01000//Interpolation;
g0100 = g01000[0.3];

g20000 = Import["renormalized_bias_functions_pvs/gp20000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.m"];
g20000 = g20000//Interpolation;
g2000 = g20000[0.3];

g11000 = Import["renormalized_bias_functions_pvs/gp11000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.m"];
g11000 = g11000//Interpolation;
g1100 = g11000[0.3];

g02000 = Import["renormalized_bias_functions_pvs/gp02000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.m"];
g02000 = g02000//Interpolation;
g0200 = g02000[0.3];

g00100 = Import["renormalized_bias_functions_pvs/gp00100_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.m"];
g00100 = g00100//Interpolation;
g0010 = g00100[0.3];

g00010 = Import["renormalized_bias_functions_pvs/gp00010_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.m"];
g00010 = g00010//Interpolation;
g0001 = g00010[0.3];

l=ReadList["background_quantity/2d_multipole_array_0_to_2e4_121bins.dat", Number, RecordLists->True]//Flatten;
lmin = l[[1]];
lmax = l[[-1]];

pks = ReadList["background_quantity/projected_2d_kappa_power_spectrum_ell0_to_2e4_121bins_enhanced_BAO.dat", Number, RecordLists -> True] // Flatten;
Pk = {l, pks} // Transpose // Interpolation;

zlospk=ReadList["background_quantity/z_los_0.05_until_1.5_100bins.dat",Number,RecordLists->True]//Flatten;
pklzgrid=Import["background_quantity/P_chi_z_grid_until_z0.05_to_1.5_100bins_l0_to_2e4_121bins_enhanced_BAO.dat"];
list={};
For[i=1,i<=Length[l],i++,nl={zlospk,pklzgrid[[i]]}//Transpose;
nl=Replace[nl,{x_,y_}:>{l[[i]],x,y},{1}];AppendTo[list,nl]];
pkgrid=list//Flatten[#,1]&;
Pkgrid=pkgrid//Interpolation;

zlos = ReadList["background_quantity/z_los_until_1.5_100bins.dat", Number, RecordLists -> True] // Flatten;
zmin = zlos[[1]];
zmax = zlos[[-1]];
chiz = ReadList["background_quantity/chi_los_until_z1.5_100bins.dat", Number, RecordLists -> True] // Flatten;
Hz = ReadList["background_quantitity/H_los_until_z1.5_100bins.dat", Number, RecordLists -> True] // Flatten;
chi = {zlos, chiz} // Transpose // Interpolation;
H = {zlos, Hz} // Transpose // Interpolation;

Xtable = Import["X_table_r_nbessel_mk_planck18_cosmology_dirac_delta_z_r1.0_dr_0.0001_ell0_to_2e4_15arcmin_smooth.m"];
Xtable = Xtable // Flatten[#, 2] &;
X = Xtable // Interpolation;

X01table = Import["X_chi_01_table_r_nbessel_mk_planck18_cosmology_dirac_delta_z_r0.3_dr_0.0001_dz_0.01_ell0_to_2e4_15arcmin_smooth.m"];
X01table = X01table // Flatten[#, 1] &;
X01 = X01table // Interpolation;
X03table = Import["X_chi_03_table_r_nbessel_mk_planck18_cosmology_dirac_delta_z_r0.3_dr_0.0001_dz_0.01_ell0_to_2e4_15arcmin_smooth.m"];
X03table = X03table // Flatten[#, 1] &;
X03 = X03table // Interpolation;
X10table = Import["X_chi_10_table_r_nbessel_mk_planck18_cosmology_dirac_delta_z_r0.3_dr_0.0001_dz_0.01_ell0_to_2e4_15arcmin_smooth.m"];
X10table = X10table // Flatten[#, 1] &;
X10 = X10table // Interpolation;
X12table = Import["X_chi_12_table_r_nbessel_mk_planck18_cosmology_dirac_delta_z_r0.3_dr_0.0001_dz_0.01_ell0_to_2e4_15arcmin_smooth.m"];
X12table = X12table // Flatten[#, 1] &;
X12 = X12table // Interpolation;
X14table = Import["X_chi_14_table_r_nbessel_mk_planck18_cosmology_dirac_delta_z_r0.3_dr_0.0001_dz_0.01_ell0_to_2e4_15arcmin_smooth.m"];
X14table = X14table // Flatten[#, 1] &;
X14 = X14table // Interpolation;
X21table = Import["X_chi_21_table_r_nbessel_mk_planck18_cosmology_dirac_delta_z_r0.3_dr_0.0001_dz_0.01_ell0_to_2e4_15arcmin_smooth.m"];
X21table = X21table // Flatten[#, 1] &;
X21 = X21table // Interpolation;
X23table = Import["X_chi_23_table_r_nbessel_mk_planck18_cosmology_dirac_delta_z_r0.3_dr_0.0001_dz_0.01_ell0_to_2e4_15arcmin_smooth.m"];
X23table = X23table // Flatten[#, 1] &;
X23 = X23table // Interpolation;
X32table = Import["X_chi_32_table_r_nbessel_mk_planck18_cosmology_dirac_delta_z_r0.3_dr_0.0001_dz_0.01_ell0_to_2e4_15arcmin_smooth.m"];
X32table = X32table // Flatten[#, 1] &;
X32 = X32table // Interpolation;
X34table = Import["X_chi_34_table_r_nbessel_mk_planck18_cosmology_dirac_delta_z_r0.3_dr_0.0001_dz_0.01_ell0_to_2e4_15arcmin_smooth.m"];
X34table = X34table // Flatten[#, 1] &;
X34 = X34table // Interpolation;
X43table = Import["X_chi_43_table_r_nbessel_mk_planck18_cosmology_dirac_delta_z_r0.3_dr_0.0001_dz_0.01_ell0_to_2e4_15arcmin_smooth.m"];
X43table = X43table // Flatten[#, 1] &;
X43 = X43table // Interpolation;

W[k_, r_] := Exp[-k^2*r^2/2];
angle = 15/60/180*Pi;

Pklinear = ParallelTable[{k, (g1000 + g0100*k^2)^2*W[k, angle]^2*Pk[k]}, {k, l}];
Export["PkLO_peak_nu0.3_correction_planck18_cosmology_dirac_delta_z_lmin0_lmax2e4_15arcmin_smooth.m", Pklinear];

Pksquare[k_]:=1/(4*Pi)*NIntegrate[r*BesselJ[0, k*r]*(g2000^2*X[r, 0, 1]^2 + 4*g1100*g2000*X[r, 0, 1]*X[r, 0, 3] + 2*g1100^2*X[r, 0, 1]*X[r, 0, 5] + (2*g1100^2 + 2*g0200*g2000 - 4*g0001*g2000)*X[r, 0, 3]^2 + (4*g0200*g1100 - 8*g0001*g1100)*X[r, 0, 5]*X[r, 0, 3] + (4*g0001^2 - 4*g0001*g0200 + g0200^2)*X[r, 0, 5]^2 + 4*g0010*g2000*(X[r, 1, 2]^2) + 8*g0010*g1100*(X[r, 1, 4]*X[r, 1, 2]) - (8*g0001*g0010 - 4*g0010*g0200)*(X[r, 1, 4]^2) + (g0010^2 + 2*g0001*g2000)*(2*X[r, 2, 3]^2 + 2*X[r, 0, 3]^2) + 4*g0001*g1100*(2*X[r, 2, 5]*X[r, 2, 3] + 2*X[r, 0, 5]*X[r, 0, 3]) - (4*g0001^2 - 2*g0001*g0200)*(2*X[r, 2, 5]^2 + 2*X[r, 0, 5]^2) + 2*g0001*g0010*(2*X[r, 3, 4]^2 + 6*X[r, 1, 4]^2) + g0001^2*(2*X[r, 4, 5]^2 + 8*X[r, 2, 5]^2 + 6*X[r, 0, 5]^2)), {r, 0.0, 1.0}, Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0, "MaxErrorIncreases" -> 2000}, MaxRecursion -> 40, WorkingPrecision -> 16, PrecisionGoal -> 5, AccuracyGoal -> Infinity];

PkBsmooth1st[k_] := (g1000 + g0100*k^2)*W[k, angle]*NIntegrate[1/H[z]*(3/2)^3*H[0.0]^6*Om0^3*(1 + z)^3/chi[z]*((chi[1.5] - chi[z])/chi[1.5])^3/2/Pi*r*BesselJ[0,k*r]*(2/7)*((-7*g0010 + 12*g1100)*X01[r, z]*X03[r, z] + (g0001 + 6*g0200)*X03[r, z]^2 + (13*g0010 - 7*g1100)*X12[r, z]^2 - 7*g1100*X10[r, z]*X14[r, z] - 7*(g0001 + g0200)*X12[r, z]*X14[r, z] + g2000*(6*X01[r, z]^2 - 7*X10[r, z]*X12[r, z] + X21[r, z]^2) - (7*g0010 - 2*g1100)*X21[r, z]*X23[r, z] + (12*g0001 + g0200)*X23[r, z]^2 + g0010*X32[r, z]^2 - 7*g0001*X32[r, z]*X34[r, z] + g0001*X43[r, z]^2), {z, 0.05, 1.5}, {r, 0.0, 0.3}, Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0, "MaxErrorIncreases" -> 2000}, MaxRecursion -> 40, WorkingPrecision -> 16, PrecisionGoal -> 5, AccuracyGoal -> Infinity];

PkBsmooth2nd[k_] := (g1000 + g0100*k^2)*W[k, angle]^2*NIntegrate[1/H[z]*(3/2)^3*H[0.0]^6*Om0^3*(1 + z)^3/chi[z]*((chi[1.5] - chi[z])/chi[1.5])^3*Pkgrid[k, z]/Pi*(10/7*k2*(BesselI[0, k*k2*angle^2]*(g2000 + 2*(g1100 + g0010)*k2^2 + g1100*k^2 + (g0200 - 2*g0001)*k2^2*k^2 + (g0200 + 2*g0001)*k2^4) - BesselI[1, k*k2*angle^2]*(2*(g1100 + g0010 - 2/angle^2*g0001)*k*k2 + 2*(g0200 + 2*g0001)*k*k2^3) + BesselI[2, k*k2*angle^2]*4*g0001*k^2*k2^2) - (k2^2/k + k)*(BesselI[1, k*k2*angle^2]*(g2000 + 2*(g1100 + g0010)*k2^2 + g1100*k^2 + (g0200 - 2*g0001)*k2^2*k^2 + (g0200 + 2*g0001)*k2^4 - 1/angle^2*(2*(g1100 + g0010) + 2*(g0200 + 2*g0001)*k2^2)) + BesselI[2, k*k2*angle^2]*(12/angle^2*g0001*k*k2 -2*(g1100 + g0010)*k*k2 - 2*(g0200 + 2*g0001)*k*k2^3) + BesselI[3, k*k2*angle^2]*4*g0001*k^2*k2^2) + 4/7*k2*(1/angle^2/k/k2*BesselI[1, k*k2*angle^2]*(g2000 + 2*(g1100 + g0010)*k2^2 + g1100*k^2 + (g0200 - 2*g0001)*k2^2*k^2 + (g0200 + 2*g0001)*k2^4) + BesselI[2, k*k2*angle^2]*(g2000 + 2*(g1100 + g0010)*k2^2 + g1100*k^2 + (g0200 - 2*g0001)*k2^2*k^2 + (g0200 + 2*g0001)*k2^4 - 6/angle^2*(g1100 + g0010) - 6/angle^2*(g0200 + 2*g0001)*k2^2 + 4*g0001*(3/angle^4 + k^2*k2^2)) - BesselI[3, k*k2*angle^2]*(2*(g1100 + g0010)*k*k2 + 2*(g0200 + 2*g0001)*k*k2^3)))*W[k2, angle]^2*Pkgrid[k2, z], {z, 0.05, 1.5}, {k2, lmin, lmax}, Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0, "MaxErrorIncreases" -> 2000}, MaxRecursion -> 40, WorkingPrecision -> 16, PrecisionGoal -> 5, AccuracyGoal -> Infinity];

PkB[k_] := PkBsmooth1st[k] + PkBsmooth2nd[k];

Pksquaretable = ParallelTable[{k, Pksquare[k]}, {k, l}];
Export["PkNLOsquare_peak_nu0.3_correction_planck18_cosmology_dirac_delta_z_lmin0_lmax2e4_r1.0_dr_0.0001_15arcmin_smooth.m", Pksquaretable];

PkBtable = ParallelTable[{k, PkB[k]}, {k, l}];
Export["PkNLOB_peak_nu0.3_correction_planck18_cosmology_dirac_delta_z_lmin0_lmax2e4_r0.3_dr_0001_z_1.5_dz_0.01_15arcmin_smooth.m", PkBtable];

Exit[];
