clear;
clear all;
%p=parpool(8)
g = MTB.geometry("Kai_2band_C");

g = MTB.read_poscar(g,"data/kai_2band_hr/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4/qe/gamma_low/wannier90_hr_p1.dat','data/TaIrTe4/qe/gamma_low/wannier90_hr_p2.dat');
% [g.ham,g.hopr] = MTB.wannier.read_hr('data/TaIrTe4_2d_tb/wannier90_hr_p1.dat'