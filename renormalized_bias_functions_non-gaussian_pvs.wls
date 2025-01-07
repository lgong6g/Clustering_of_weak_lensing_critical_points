#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Input:: *)
(*Om0=0.3143878052130123;*)


(* ::Input:: *)
(*W[k_,r_]:=Exp[-k^2*r^2/2];*)
(*angle=15/60/180*Pi;*)
(*(*W[k_,r_]:=Exp[-k^2*r^2/2];*)
(*R=20;*)*)


(* ::Input:: *)
(*l=ReadList["/home/laurence/Documents/personal_work/peak_correlation/2d_multipole_array_0_to_2e4_121bins.dat", Number, RecordLists->True]//Flatten;*)
(*(*k3d=ReadList["/home/laurence/Documents/personal_work/peak_correlation/3d_wavenumber_k.dat", Number, RecordLists->True]//Flatten;*)*)


(* ::Input:: *)
(*lmin=l[[1]];*)
(*lmax=l[[-1]];*)
(*(*kmax=k3d[[-1]];*)
(*kmin=k3d[[1]];*)*)


(* ::Input:: *)
(*pks=ReadList["/home/laurence/Documents/personal_work/peak_correlation/zslice0.75_2d_kappa_power_spectrum_ell0_to_2e4_121bins.dat",Number,RecordLists->True]//Flatten;*)
(*(*pk3d=ReadList["/home/laurence/Documents/personal_work/peak_correlation/3d_matter_power_spectrum_at_z0.dat",Number,RecordLists->True]//Flatten;*)*)


(* ::Input:: *)
(*Pk={l,pks}//Transpose//Interpolation;*)
(*(*Pk3d={k3d,pk3d}//Transpose//Interpolation;*)*)


(* ::Input:: *)
(*zlos=ReadList["/home/laurence/Documents/personal_work/peak_correlation/z_los_until_1.5_100bins.dat",Number,RecordLists->True]//Flatten;*)
(*zmin=zlos[[1]];*)
(*zmax=zlos[[-1]];*)
(*chiz=ReadList["/home/laurence/Documents/personal_work/peak_correlation/chi_los_until_z1.5_100bins.dat",Number,RecordLists->True]//Flatten;*)
(*Hz=ReadList["/home/laurence/Documents/personal_work/peak_correlation/H_los_until_z1.5_100bins.dat",Number,RecordLists->True]//Flatten;*)
(*chi={zlos,chiz}//Transpose//Interpolation;*)
(*H={zlos,Hz}//Transpose//Interpolation*)


(* ::Input:: *)
(*zlospk=ReadList["/home/laurence/Documents/personal_work/peak_correlation/z_los_0.05_until_1.5_100bins.dat",Number,RecordLists->True]//Flatten;*)
(*pklzgrid=Import["/home/laurence/Documents/personal_work/peak_correlation/P_chi_z_grid_until_z0.05_to_1.5_100bins_l0_to_2e4_121bins.dat"];*)
(*list={};*)
(*For[i=1,i<=Length[l],i++,nl={zlospk,pklzgrid[[i]]}//Transpose;*)
(*nl=Replace[nl,{x_,y_}:>{l[[i]],x,y},{1}];AppendTo[list,nl]]*)
(*pkgrid=list//Flatten[#,1]&;*)
(*Pkgrid=pkgrid//Interpolation*)


(* ::Input:: *)
(*sigma2[n_]:=NIntegrate[(1/2/Pi)*ell^(2*n+1)*W[ell,angle]^2*Pk[ell],{ell,lmin,lmax}];*)


(* ::Input:: *)
(*(*sigma2[n_]:=NIntegrate[1/H[z]*(3/2)^2*H[0.0]^4*Om0^2*(1+z)^2*((chi[1.5]-chi[z])/chi[1.5])^2*ell^(2*n+1)/(2*Pi)*W[ell,angle]^2*Pkgrid[ell,z],{z,0.05,1.5},{ell,lmin,lmax}];*)
(*sigma2ki[n_,i_]:=NIntegrate[(1/2/Pi)*ell^(2*n+i+1)*W[ell,angle]^2,{ell,0,Infinity}];*)*)
(*(*sigma2kiv2[n_,i_]:=NIntegrate[(1/2/Pi)*ell^(2*n+i+1)*W[ell,angle]^2,{ell,1,Infinity}];*)*)
(*(*sigma3d2[n_]:=NIntegrate[4*Pi/(2*Pi)^3*Pk3d[k]*k^(2*n+2)*W[k,R]^2,{k,kmin,kmax}];*)*)


(* ::Input:: *)
(*sigma20=sigma2[0];*)
(*sigma21=sigma2[1];*)
(*sigma22=sigma2[2];*)


(* ::Input:: *)
(*sigma20ki=sigma2ki[0,-1];*)
(*sigma21ki=sigma2ki[1,-1];*)
(*sigma22ki=sigma2ki[2,-1];*)
(*(*sigma20kiv2=sigma2kiv2[0,-3];*)
(*sigma21kiv2=sigma2kiv2[1,-3];*)
(*sigma22kiv2=sigma2kiv2[2,-3];*)*)


(* ::Input:: *)
(*sigma20ki*)


(* ::Input:: *)
(*sigma22*)


(* ::Input:: *)
(*(*sigma22d0=sigma22d[0];*)
(*sigma22d1=sigma22d[1];*)
(*sigma22d2=sigma22d[2];*)*)


(* ::Input:: *)
(*(*sigma3d20=sigma3d2[0];*)
(*sigma3d21=sigma3d2[1];*)
(*sigma3d22=sigma3d2[2];*)*)


(* ::Input:: *)
(*gamma=sigma21/Sqrt[sigma20*sigma22];*)
(*(*gammaki=sigma21ki/Sqrt[sigma20ki*sigma22ki];*)*)
(*(*gamma3d=sigma3d21/Sqrt[sigma3d20*sigma3d22];*)*)


(* ::Input:: *)
(*gammaki*)


(* ::Input:: *)
(*X0=1;*)
(*X1=-1;*)


(* ::Input:: *)
(*gpdenomki[nu_,i_,j_,k_,l_,m_]:=Power[sigma20ki, i/2]*Power[sigma21ki, k]*Power[sigma22ki, j/2+l+3*m/2]*NIntegrate[1/(2*Pi*Sqrt[1-gammaki^2])*Exp[-(y^2+x^2-2*gammaki*x*y)/2/(1-gammaki^2)]*(Exp[-x^2]+x^2-1), {x, 0, Infinity},{y, nu, Infinity}];*)
(*gp10000ki[nu_]:=NIntegrate[X0*1/(2*Pi*Sqrt[1-gammaki^2])*Exp[-(nu^2+x^2-2*gammaki*x*nu)/2/(1-gammaki^2)]*(Exp[-x^2]+x^2-1), {x, 0, Infinity}]/gpdenomki[nu,1,0,0,0,0];*)
(*gp01000ki[nu_]:=NIntegrate[X0*1/(2*Pi*Sqrt[1-gammaki^2])*Exp[-(y^2+x^2-2*gammaki*x*y)/2/(1-gammaki^2)]*(x-y*gammaki)/(1-gammaki^2)*(Exp[-x^2]+x^2-1), {x, 0, Infinity},{y, nu, Infinity}]/gpdenomki[nu,0,1,0,0,0];*)


(* ::Input:: *)
(*g1[\[Nu]_,ell_]:=gp10000ki[\[Nu]]+gp01000ki[\[Nu]]*ell^2;*)


(* ::Input:: *)
(*g1[6,100]*)


(* ::Input:: *)
(*gpdenom[nu_,i_,j_,k_,l_,m_]:=Power[sigma20, i/2]*Power[sigma21, k]*Power[sigma22, j/2+l+3*m/2]*NIntegrate[1/(2*Pi*Sqrt[1-gamma^2])*Exp[-(y^2+x^2-2*gamma*x*y)/2/(1-gamma^2)]*(Exp[-x^2]+x^2-1), {x, 0, Infinity},{y, nu, Infinity}];*)
(*gp10000[nu_]:=NIntegrate[X0*1/(2*Pi*Sqrt[1-gamma^2])*Exp[-(nu^2+x^2-2*gamma*x*nu)/2/(1-gamma^2)]*(Exp[-x^2]+x^2-1), {x, 0, Infinity}]/gpdenom[nu,1,0,0,0,0];*)
(*gp01000[nu_]:=NIntegrate[X0*1/(2*Pi*Sqrt[1-gamma^2])*Exp[-(y^2+x^2-2*gamma*x*y)/2/(1-gamma^2)]*(x-y*gamma)/(1-gamma^2)*(Exp[-x^2]+x^2-1), {x, 0, Infinity},{y, nu, Infinity}]/gpdenom[nu,0,1,0,0,0];*)
(*gp20000[nu_]:=NIntegrate[X0*1/(2*Pi*Sqrt[1-gamma^2])*Exp[-(nu^2+x^2-2*gamma*x*nu)/2/(1-gamma^2)]*(nu-x*gamma)/(1-gamma^2)*(Exp[-x^2]+x^2-1), {x, 0, Infinity}]/gpdenom[nu,2,0,0,0,0];*)
(*gp11000[nu_]:=NIntegrate[X0*1/(2*Pi*Sqrt[1-gamma^2])*Exp[-(nu^2+x^2-2*gamma*x*nu)/2/(1-gamma^2)]*(x-nu*gamma)/(1-gamma^2)*(Exp[-x^2]+x^2-1), {x, 0, Infinity}]/gpdenom[nu,1,1,0,0,0];*)
(*gp02000[nu_]:=NIntegrate[X0*1/(2*Pi*Sqrt[1-gamma^2])*Exp[-(y^2+x^2-2*gamma*x*y)/2/(1-gamma^2)]*((1+y^2)*gamma^2+x^2-1-2*gamma*y*x)/(1-gamma^2)^2*(Exp[-x^2]+x^2-1),{x, 0, Infinity},{y, nu, Infinity}]/gpdenom[nu,0,2,0,0,0];*)
(*gp00100[nu_]:=NIntegrate[X1*1/(2*Pi*Sqrt[1-gamma^2])*Exp[-(y^2+x^2-2*gamma*x*y)/2/(1-gamma^2)]*(Exp[-x^2]+x^2-1),{x, 0, Infinity},{y, nu, Infinity}]/gpdenom[nu,0,0,1,0,0];*)
(*gp00010[nu_]:=NIntegrate[X0*1/(2*Pi*Sqrt[1-gamma^2])*Exp[-(y^2+x^2-2*gamma*x*y)/2/(1-gamma^2)]*((1+x^2)*Exp[-x^2]-1),{x, 0, Infinity},{y, nu, Infinity}]/gpdenom[nu,0,0,0,1,0];*)


(* ::Input:: *)
(*gp10000[0.3]+gp01000[0.3]*100^2*)


(* ::Input:: *)
(*10.0/Sqrt[sigma20]*)


(* ::Input:: *)
(*Nu=Range[-3,6,0.1];*)


(* ::Input:: *)
(*gp1000=ParallelTable[{nu,gp10000[nu]},{nu,Nu}];*)
(*gp0100=ParallelTable[{nu,gp01000[nu]},{nu,Nu}];*)
(*gp2000=ParallelTable[{nu,gp20000[nu]},{nu,Nu}];*)
(*gp1100=ParallelTable[{nu,gp11000[nu]},{nu,Nu}];*)
(*gp0200=ParallelTable[{nu,gp02000[nu]},{nu,Nu}];*)
(*gp0010=ParallelTable[{nu,gp00100[nu]},{nu,Nu}];*)
(*gp0001=ParallelTable[{nu,gp00010[nu]},{nu,Nu}];*)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gp10000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth_zslice0.75.m",gp1000];*)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gp01000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth_zslice0.75.m",gp0100];*)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gp20000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth_zslice0.75.m",gp2000];*)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gp11000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth_zslice0.75.m",gp1100];*)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gp02000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth_zslice0.75.m",gp0200];*)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gp00100_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth_zslice0.75.m",gp0010];*)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gp00010_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth_zslice0.75.m",gp0001];*)


(* ::Input:: *)
(*(*gpdenom3d[nu_,i_,j_,k_,l_,m_]:=Power[sigma3d20, i/2]*Power[sigma3d21, k]*Power[sigma3d22, j/2+l+3*m/2]*NIntegrate[1/(2*Pi*Sqrt[1-gamma3d^2])*Exp[-(y^2+x^2-2*gamma3d*x*y)/2/(1-gamma3d^2)]*((x/2)*(x^2-3)*(Erf[0.5*Sqrt[5/2]*x]+Erf[Sqrt[5/2]*x])+Sqrt[2/5/Pi]*((x^2/2-8/5)*Exp[-5*x^2/2]+(31/4*x^2+8/5)*Exp[-5*x^2/8])), {x, 0, Infinity},{y, nu, Infinity}];*)
(*gp100003d[nu_]:=NIntegrate[X0*1/(2*Pi*Sqrt[1-gamma3d^2])*Exp[-(nu^2+x^2-2*gamma3d*x*nu)/2/(1-gamma3d^2)]*((x/2)*(x^2-3)*(Erf[0.5*Sqrt[5/2]*x]+Erf[Sqrt[5/2]*x])+Sqrt[2/5/Pi]*((x^2/2-8/5)*Exp[-5*x^2/2]+(31/4*x^2+8/5)*Exp[-5*x^2/8])), {x, 0, Infinity}]/gpdenom3d[nu,1,0,0,0,0];*)*)


(* ::Input:: *)
(*gp1000 = Import["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gp10000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.m"];*)
(*gp0100 = Import["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gp01000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.m"];*)
(*gp2000 = Import["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gp20000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.m"];*)
(*gp1100 = Import["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gp11000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.m"];*)
(*gp0200 = Import["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gp02000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.m"];*)
(*gp0010 = Import["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gp00100_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.m"];*)
(*gp0001 = Import["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gp00010_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.m"];*)


(* ::Input:: *)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gp10000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.csv",gp1000];*)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gp01000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.csv",gp0100];*)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gp20000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.csv",gp2000];*)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gp11000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.csv",gp1100];*)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gp02000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.csv",gp0200];*)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gp00100_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.csv",gp0010];*)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gp00010_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.csv",gp0001];*)


(* ::Input:: *)
(*ListPlot[gp0001]*)


(* ::Input:: *)
(*gvdenom[nu_,i_,j_,k_,l_,m_]:=Power[sigma20, i/2]*Power[sigma21, k]*Power[sigma22, j/2+l+3*m/2]*NIntegrate[1/(2*Pi*Sqrt[1-gamma^2])*Exp[-(y^2+x^2-2*gamma*x*y)/2/(1-gamma^2)]*(Exp[-x^2]+x^2-1), {x,-Infinity,0},{y,nu,Infinity}];*)
(*gv10000[nu_]:=NIntegrate[X0*1/(2*Pi*Sqrt[1-gamma^2])*Exp[-(nu^2+x^2-2*gamma*x*nu)/2/(1-gamma^2)]*(Exp[-x^2]+x^2-1), {x, -Infinity,0}]/gvdenom[nu,1,0,0,0,0];*)
(*gv01000[nu_]:=NIntegrate[X0*1/(2*Pi*Sqrt[1-gamma^2])*Exp[-(y^2+x^2-2*gamma*x*y)/2/(1-gamma^2)]*(x-y*gamma)/(1-gamma^2)*(Exp[-x^2]+x^2-1), {x, -Infinity,0},{y,nu,Infinity}]/gvdenom[nu,0,1,0,0,0];*)
(*gv20000[nu_]:=NIntegrate[X0*1/(2*Pi*Sqrt[1-gamma^2])*Exp[-(nu^2+x^2-2*gamma*x*nu)/2/(1-gamma^2)]*(nu-x*gamma)/(1-gamma^2)*(Exp[-x^2]+x^2-1), {x, -Infinity,0}]/gvdenom[nu,2,0,0,0,0];*)
(*gv11000[nu_]:=NIntegrate[X0*1/(2*Pi*Sqrt[1-gamma^2])*Exp[-(nu^2+x^2-2*gamma*x*nu)/2/(1-gamma^2)]*(x-nu*gamma)/(1-gamma^2)*(Exp[-x^2]+x^2-1), {x, -Infinity,0}]/gvdenom[nu,1,1,0,0,0];*)
(*gv02000[nu_]:=NIntegrate[X0*1/(2*Pi*Sqrt[1-gamma^2])*Exp[-(y^2+x^2-2*gamma*x*y)/2/(1-gamma^2)]*((1+y^2)*gamma^2+x^2-1-2*gamma*y*x)/(1-gamma^2)^2*(Exp[-x^2]+x^2-1),{x, -Infinity,0},{y,nu,Infinity}]/gvdenom[nu,0,2,0,0,0];*)
(*gv00100[nu_]:=NIntegrate[X1*1/(2*Pi*Sqrt[1-gamma^2])*Exp[-(y^2+x^2-2*gamma*x*y)/2/(1-gamma^2)]*(Exp[-x^2]+x^2-1),{x, -Infinity,0},{y,nu,Infinity}]/gvdenom[nu,0,0,1,0,0];*)
(*gv00010[nu_]:=NIntegrate[X0*1/(2*Pi*Sqrt[1-gamma^2])*Exp[-(y^2+x^2-2*gamma*x*y)/2/(1-gamma^2)]*((1+x^2)*Exp[-x^2]-1),{x,-Infinity,0},{y, nu,Infinity}]/gvdenom[nu,0,0,0,1,0];*)


(* ::Input:: *)
(*gv1000=ParallelTable[{nu,gv10000[nu]},{nu,Nu}];*)
(*gv0100=ParallelTable[{nu,gv01000[nu]},{nu,Nu}];*)
(*gv2000=ParallelTable[{nu,gv20000[nu]},{nu,Nu}];*)
(*gv1100=ParallelTable[{nu,gv11000[nu]},{nu,Nu}];*)
(*gv0200=ParallelTable[{nu,gv02000[nu]},{nu,Nu}];*)
(*gv0010=ParallelTable[{nu,gv00100[nu]},{nu,Nu}];*)
(*gv0001=ParallelTable[{nu,gv00010[nu]},{nu,Nu}];*)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gv10000_table_planck18_cosmology_nu-3to6_dnu0.1_20arcmin_smooth.m",gv1000];*)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gv01000_table_planck18_cosmology_nu-3to6_dnu0.1_20arcmin_smooth.m",gv0100];*)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gv20000_table_planck18_cosmology_nu-3to6_dnu0.1_20arcmin_smooth.m",gv2000];*)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gv11000_table_planck18_cosmology_nu-3to6_dnu0.1_20arcmin_smooth.m",gv1100];*)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gv02000_table_planck18_cosmology_nu-3to6_dnu0.1_20arcmin_smooth.m",gv0200];*)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gv00100_table_planck18_cosmology_nu-3to6_dnu0.1_20arcmin_smooth.m",gv0010];*)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gv00010_table_planck18_cosmology_nu-3to6_dnu0.1_20arcmin_smooth.m",gv0001];*)


(* ::Input:: *)
(*gv1000 = Import["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gv10000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.m"];*)
(*gv0100 = Import["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gv01000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.m"];*)
(*gv2000 = Import["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gv20000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.m"];*)
(*gv1100 = Import["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gv11000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.m"];*)
(*gv0200 = Import["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gv02000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.m"];*)
(*gv0010 = Import["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gv00100_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.m"];*)
(*gv0001 = Import["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gv00010_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.m"];*)


(* ::Input:: *)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gv10000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.csv",gv1000];*)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gv01000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.csv",gv0100];*)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gv20000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.csv",gv2000];*)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gv11000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.csv",gv1100];*)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gv02000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.csv",gv0200];*)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gv00100_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.csv",gv0010];*)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gv00010_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.csv",gv0001];*)


(* ::Input:: *)
(*ListPlot[gv2000]*)


(* ::Input:: *)
(*gsdenom[nu_,i_,j_,k_,l_,m_]:=Power[sigma20, i/2]*Power[sigma21, k]*Power[sigma22, j/2+l+3*m/2]*NIntegrate[1/(2*Pi*Sqrt[1-gamma^2])*Exp[-(y^2+x^2-2*gamma*x*y)/2/(1-gamma^2)]*(-Exp[-x^2]), {x,-Infinity,Infinity},{y, nu, Infinity},Method->{"GlobalAdaptive","SymbolicProcessing"->0,"MaxErrorIncreases"->2000},MaxRecursion->40,WorkingPrecision->16,PrecisionGoal->5,AccuracyGoal->Infinity];*)
(*gs10000[nu_]:=NIntegrate[X0*1/(2*Pi*Sqrt[1-gamma^2])*Exp[-(nu^2+x^2-2*gamma*x*nu)/2/(1-gamma^2)]*(-Exp[-x^2]), {x, -Infinity,Infinity},Method->{"GlobalAdaptive","SymbolicProcessing"->0,"MaxErrorIncreases"->2000},MaxRecursion->40,WorkingPrecision->16,PrecisionGoal->5,AccuracyGoal->Infinity]/gsdenom[nu,1,0,0,0,0];*)
(*gs01000[nu_]:=NIntegrate[X0*1/(2*Pi*Sqrt[1-gamma^2])*Exp[-(y^2+x^2-2*gamma*x*y)/2/(1-gamma^2)]*(x-y*gamma)/(1-gamma^2)*(-Exp[-x^2]), {x, -Infinity,Infinity},{y, nu, Infinity},Method->{"GlobalAdaptive","SymbolicProcessing"->0,"MaxErrorIncreases"->2000},MaxRecursion->40,WorkingPrecision->16,PrecisionGoal->5,AccuracyGoal->Infinity]/gsdenom[nu,0,1,0,0,0];*)
(*gs20000[nu_]:=NIntegrate[X0*1/(2*Pi*Sqrt[1-gamma^2])*Exp[-(nu^2+x^2-2*gamma*x*nu)/2/(1-gamma^2)]*(nu-x*gamma)/(1-gamma^2)*(-Exp[-x^2]), {x, -Infinity,Infinity},Method->{"GlobalAdaptive","SymbolicProcessing"->0,"MaxErrorIncreases"->2000},MaxRecursion->40,WorkingPrecision->8,PrecisionGoal->5,AccuracyGoal->Infinity]/gsdenom[nu,2,0,0,0,0];*)
(*gs11000[nu_]:=NIntegrate[X0*1/(2*Pi*Sqrt[1-gamma^2])*Exp[-(nu^2+x^2-2*gamma*x*nu)/2/(1-gamma^2)]*(x-nu*gamma)/(1-gamma^2)*(-Exp[-x^2]), {x, -Infinity,Infinity},Method->{"GlobalAdaptive","SymbolicProcessing"->0,"MaxErrorIncreases"->2000},MaxRecursion->40,WorkingPrecision->16,PrecisionGoal->5,AccuracyGoal->Infinity]/gsdenom[nu,1,1,0,0,0];*)
(*gs02000[nu_]:=NIntegrate[X0*1/(2*Pi*Sqrt[1-gamma^2])*Exp[-(y^2+x^2-2*gamma*x*y)/2/(1-gamma^2)]*((1+y^2)*gamma^2+x^2-1-2*gamma*y*x)/(1-gamma^2)^2*(-Exp[-x^2]),{x, -Infinity,Infinity},{y, nu, Infinity},Method->{"GlobalAdaptive","SymbolicProcessing"->0,"MaxErrorIncreases"->2000},MaxRecursion->40,WorkingPrecision->16,PrecisionGoal->5,AccuracyGoal->Infinity]/gsdenom[nu,0,2,0,0,0];*)
(*gs00100[nu_]:=NIntegrate[X1*1/(2*Pi*Sqrt[1-gamma^2])*Exp[-(y^2+x^2-2*gamma*x*y)/2/(1-gamma^2)]*(-Exp[-x^2]),{x, -Infinity,Infinity},{y, nu, Infinity}]/gsdenom[nu,0,0,1,0,0];*)
(*gs00010[nu_]:=NIntegrate[X0*1/(2*Pi*Sqrt[1-gamma^2])*Exp[-(y^2+x^2-2*gamma*x*y)/2/(1-gamma^2)]*(-(1+x^2)*Exp[-x^2]),{x,-Infinity,Infinity},{y, nu, Infinity},Method->{"GlobalAdaptive","SymbolicProcessing"->0,"MaxErrorIncreases"->2000},MaxRecursion->40,WorkingPrecision->16,PrecisionGoal->5,AccuracyGoal->Infinity]/gsdenom[nu,0,0,0,1,0];*)


(* ::Input:: *)
(*gs1000=ParallelTable[{nu,gs10000[nu]},{nu,Nu}];*)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gs10000_table_planck18_cosmology_nu-3to6_dnu0.1_20arcmin_smooth.m",gs1000];*)
(**)


(* ::Input:: *)
(*gs0100=ParallelTable[{nu,gs01000[nu]},{nu,Nu}];*)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gs01000_table_planck18_cosmology_nu-3to6_dnu0.1_20arcmin_smooth.m",gs0100];*)


(* ::Input:: *)
(*gs0100[[All,2]]=0.0*)


(* ::Input:: *)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gs01000_table_planck18_cosmology_nu-3to6_dnu0.1_H-1j.m",gs0100];*)


(* ::Input:: *)
(*ListPlot[gs2000]*)


(* ::Input:: *)
(*gs2000=ParallelTable[{nu,gs20000[nu]},{nu,Nu}];*)
(*gs2000[[31,2]]=0.0;*)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gs20000_table_planck18_cosmology_nu-3to6_dnu0.1_20arcmin_smooth.m",gs2000];*)
(**)


(* ::Input:: *)
(*gs1100=ParallelTable[{nu,gs11000[nu]},{nu,Nu}];*)
(*gs1100[[31,2]]=0.0;*)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gs11000_table_planck18_cosmology_nu-3to6_dnu0.1_20arcmin_smooth.m",gs1100];*)
(**)


(* ::Input:: *)
(*gs0200=ParallelTable[{nu,gs02000[nu]},{nu,Nu}];*)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gs02000_table_planck18_cosmology_nu-3to6_dnu0.1_20arcmin_smooth.m",gs0200];*)


(* ::Input:: *)
(*gs0010=ParallelTable[{nu,gs00100[nu]},{nu,Nu}];*)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gs00100_table_planck18_cosmology_nu-3to6_dnu0.1_20arcmin_smooth.m",gs0010];*)


(* ::Input:: *)
(*gs0001=ParallelTable[{nu,gs00010[nu]},{nu,Nu}];*)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gs00010_table_planck18_cosmology_nu-3to6_dnu0.1_20arcmin_smooth.m",gs0001];*)


(* ::Input:: *)
(*gs1000 = Import["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gs10000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.m"];*)
(*gs0100 = Import["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gs01000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.m"];*)
(*gs2000 = Import["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gs20000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.m"];*)
(*gs1100 = Import["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gs11000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.m"];*)
(*gs0200 = Import["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gs02000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.m"];*)
(*gs0010 = Import["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gs00100_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.m"];*)
(*gs0001 = Import["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gs00010_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.m"];*)


(* ::Input:: *)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gs10000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.csv",gs1000];*)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gs01000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.csv",gs0100];*)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gs20000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.csv",gs2000];*)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gs11000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.csv",gs1100];*)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gs02000_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.csv",gs0200];*)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gs00100_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.csv",gs0010];*)
(*Export["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gs00010_table_planck18_cosmology_nu-3to6_dnu0.1_15arcmin_smooth.csv",gs0001];*)


(* ::Input:: *)
(*ListPlot[{gp1000ls,gp1000,gv1000,gs1000},AxesLabel->{nu,"g1000"},PlotLegends->{"peak larger smooth","peak","void","saddle"},PlotRange->All]*)


(* ::Input:: *)
(*gp1000ls=Import["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gp10000_table_planck18_cosmology_nu-3to6_dnu0.1_30arcmin_smooth.m"];*)
(*gp1000ls=gp1000ls//Interpolation*)


(* ::Input:: *)
(*gv1000[0.3]*)


(* ::Input:: *)
(*gv0100=Import["/home/laurence/Documents/personal_work/peak_correlation/renormalized_bias_functions_pvs/gv01000_table_planck18_cosmology_nu-3to6_dnu0.1.m"];*)
(*gv0100=gv0100//Interpolation*)


(* ::Input:: *)
(*gv0100[0.3]*)


(* ::Input:: *)
(**)
