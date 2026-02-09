clear;
clear all;
%parpool(8)
g = MTB.geometry("TaIrTe4");
g = MTB.read_poscar(g,"data/TaIrTe4/POSCAR-TaIrTe4");
[g.ham,g.hopr] = MTB.read_hr('data/TaIrTe4/wannier90_hr.dat');
efermi=7.3742;
%% calculate slab plane bands
MillerIndices=[0,0,1];
Umatrix=g.MillerIndicestoumatrix(MillerIndices);
Urot=g.surfab;
knum=21;
nslab=3;
Occ=120*nslab;

kxline=[0,1];
kyline=[0,1];
[Kx,Ky] = get_Slab2Dkmesh(g,kxline,kyline,knum);

[~,Enk]=MTB.ham.get_slab_plane_bands(g,Kx,Ky,nslab);

Enk=Enk-efermi;
%%
filename='TaIrTe4_Enk-400x400.dat'
writeEnk(Enk,Kx,Ky,Occ,filename)

save("TaIrTe4_Enk-400x400.mat","Enk","-v7.3")

%% Enk(nk,nk,nband) Kx(nk,nk) Ky(nk,nk)
function writeEnk(Enk,Kx,Ky,Occ,filename)
    file=fopen(filename, 'w');
    fprintf(file, '%s\n', '# kx     ky     kz     Ev5     Ev4     Ev3     Ev2     Ev1     Ec1     Ec2     Ec3     Ec4     Ec5');
    for i=1:size(Enk,1)
        for j=1:size(Enk,2)
            fprintf(file, [repmat('%12.6f',1,13) '\n'],...
                Kx(i,j), Ky(i,j), 0.0, ...
                Enk(i,j,Occ-4), Enk(i,j,Occ-3), Enk(i,j,Occ-2), Enk(i,j,Occ-1), Enk(i,j,Occ), ...
                Enk(i,j,Occ+1), Enk(i,j,Occ+2), Enk(i,j,Occ+3), Enk(i,j,Occ+4), Enk(i,j,Occ+5));
        end
        fprintf(file,'\n');
    end
    fclose(file);
end