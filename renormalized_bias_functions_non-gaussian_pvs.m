Om0 = 0.3143878052130123;
W[k_, r_] := Exp[-k^2*r^2/2];
angle = 15/60/180*Pi;
l = ReadList["background_quantity/2d_multipole_array_0_to_2e4_121bins.dat", Number, RecordLists -> True] //Flatten;
lmin = l[[1]];
lmax = l[[-1]];
pks = ReadList["background_quantity/projected_2d_kappa_power_spectrum_ell0_to_2e4_121bins.dat", Number,RecordLists -> True] // Flatten;
Pk = {l, pks} // Transpose // Interpolation;

zlos = ReadList["background_quantity/peak_correlation/z_los_until_1.5_100bins.dat", Number, RecordLists -> True] // Flatten;
zmin = zlos[[1]];
zmax = zlos[[-1]];
chiz = ReadList["background_quantity/chi_los_until_z1.5_100bins.dat", Number, RecordLists -> True] // Flatten;
Hz = ReadList["background_quantity/H_los_until_z1.5_100bins.dat", Number, RecordLists -> True] // Flatten;
chi = {zlos, chiz} // Transpose // Interpolation;
H = {zlos, Hz} // Transpose // Interpolation;

zlospk = ReadList["background_quantity/z_los_0.05_until_1.5_100bins.dat", Number, RecordLists -> True] // Flatten;
pklzgrid = Import["background_quantity/P_chi_z_grid_until_z0.05_to_1.5_100bins_l0_to_2e4_121bins.dat"];
list = {};
For[i = 1, i <= Length[l], i++,
 nl = {zlospk, pklzgrid[[i]]} // Transpose;
 nl = Replace[nl, {x_, y_} :> {l[[i]], x, y}, {1}]; AppendTo[list, nl]]
pkgrid = list // Flatten[#, 1] &;
Pkgrid = pkgrid // Interpolation;

sigma2[n_] := NIntegrate[1/H[z]*(3/2)^2*H[0.0]^4*Om0^2*(1 + z)^2*((chi[1.5] - chi[z])/chi[1.5])^2*ell^(2*n + 1)/(2*Pi)*W[ell, angle]^2*Pkgrid[ell, z], {z, 0.05, 1.5}, {ell, lmin, lmax}];

sigma20 = sigma2[0];
sigma21 = sigma2[1];
sigma22 = sigma2[2];

gamma = sigma21/Sqrt[sigma20*sigma22];
X0 = 1;
X1 = -1;

gpdenom[nu_, i_, j_, k_, l_, m_] := Power[sigma20, i/2]*Power[sigma21, k]*Power[sigma22, j/2 + l + 3*m/2]*NIntegrate[1/(2*Pi*Sqrt[1 - gamma^2])*Exp[-(y^2 + x^2 - 2*gamma*x*y)/2/(1 - gamma^2)]*(Exp[-x^2] + x^2 - 1), {x, 0, Infinity}, {y, nu, Infinity}];
gp10000[nu_] := NIntegrate[X0*1/(2*Pi*Sqrt[1 - gamma^2])*Exp[-(nu^2 + x^2 - 2*gamma*x*nu)/2/(1 - gamma^2)]*(Exp[-x^2] + x^2 - 1), {x, 0, Infinity}]/gpdenom[nu, 1, 0, 0, 0, 0];
gp01000[nu_] := NIntegrate[X0*1/(2*Pi*Sqrt[1 - gamma^2])*Exp[-(y^2 + x^2 - 2*gamma*x*y)/2/(1 - gamma^2)]*(x - y*gamma)/(1 - gamma^2)*(Exp[-x^2] + x^2 - 1), {x, 0, Infinity}, {y, nu, Infinity}]/gpdenom[nu, 0, 1, 0, 0, 0];
gp20000[nu_] := NIntegrate[X0*1/(2*Pi*Sqrt[1 - gamma^2])*Exp[-(nu^2 + x^2 - 2*gamma*x*nu)/2/(1 - gamma^2)]*(nu - x*gamma)/(1 - gamma^2)*(Exp[-x^2] + x^2 - 1), {x, 0, Infinity}]/gpdenom[nu, 2, 0, 0, 0, 0];
gp11000[nu_] := NIntegrate[X0*1/(2*Pi*Sqrt[1 - gamma^2])*Exp[-(nu^2 + x^2 - 2*gamma*x*nu)/2/(1 - gamma^2)]*(x - nu*gamma)/(1 - gamma^2)*(Exp[-x^2] + x^2 - 1), {x, 0, Infinity}]/gpdenom[nu, 1, 1, 0, 0, 0];
gp02000[nu_] := NIntegrate[X0*1/(2*Pi*Sqrt[1 - gamma^2])*Exp[-(y^2 + x^2 - 2*gamma*x*y)/2/(1 - gamma^2)]*((1 + y^2)*gamma^2 + x^2 - 1 - 2*gamma*y*x)/(1 - gamma^2)^2*(Exp[-x^2] + x^2 - 1), {x, 0, Infinity}, {y, nu, Infinity}]/gpdenom[nu, 0, 2, 0, 0, 0];
gp00100[nu_] := NIntegrate[X1*1/(2*Pi*Sqrt[1 - gamma^2])*Exp[-(y^2 + x^2 - 2*gamma*x*y)/2/(1 - gamma^2)]*(Exp[-x^2] + x^2 - 1), {x, 0, Infinity}, {y, nu, Infinity}]/gpdenom[nu, 0, 0, 1, 0, 0];
gp00010[nu_] := NIntegrate[X0*1/(2*Pi*Sqrt[1 - gamma^2])*Exp[-(y^2 + x^2 - 2*gamma*x*y)/2/(1 - gamma^2)]*((1 + x^2)*Exp[-x^2] - 1), {x, 0, Infinity}, {y, nu, Infinity}]/gpdenom[nu, 0, 0, 0, 1, 0];

Nu = Range[-3, 6, 0.1];

gp1000 = ParallelTable[{nu, gp10000[nu]}, {nu, Nu}];
gp0100 = ParallelTable[{nu, gp01000[nu]}, {nu, Nu}];
gp2000 = ParallelTable[{nu, gp20000[nu]}, {nu, Nu}];
gp1100 = ParallelTable[{nu, gp11000[nu]}, {nu, Nu}];
gp0200 = ParallelTable[{nu, gp02000[nu]}, {nu, Nu}];
gp0010 = ParallelTable[{nu, gp00100[nu]}, {nu, Nu}];
gp0001 = ParallelTable[{nu, gp00010[nu]}, {nu, Nu}];
Export["renormalized_bias_functions_pvs/gp10000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.m", gp1000];
Export["renormalized_bias_functions_pvs/gp01000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.m", gp0100];
Export["renormalized_bias_functions_pvs/gp20000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.m", gp2000];
Export["renormalized_bias_functions_pvs/gp11000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.m", gp1100];
Export["renormalized_bias_functions_pvs/gp02000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.m", gp0200];
Export["renormalized_bias_functions_pvs/gp00100_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.m", gp0010];
Export["renormalized_bias_functions_pvs/gp00010_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.m", gp0001];
