clear;
clear all;
%p=parpool(8)
g = MTB.geometry("Gra");
g = MTB.read_poscar(g,"data/SRO/1111/POSCAR");
[g.ham,g.hopr] = MTB.wannier.read_hr('data/SRO/1111/wannier90_hr_p1.dat','data/SRO/1111/wannier90_hr_p2.dat');
%%
filename="data/SRO/1111/wannier90_hr_1111.dat";
MTB.write_hr(g,filename)