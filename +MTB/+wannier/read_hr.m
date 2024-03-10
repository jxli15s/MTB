function [hamiltonian,hopping_r] = read_hr(filename)
    %read wannier90_hr.dat
    %Return:
    %hamiltonian: real space hopping parameters
    %hopping_r: real space hopping coordinates
    wannier90_hrdata=importdata(filename);
    textData=wannier90_hrdata.textdata;
    numData=wannier90_hrdata.data;
    nbands=numData(1);
    nrpts=numData(2);
    degrpts=numData(3:2+nrpts);
    numData(1:nrpts+2)=[];
    hoppingdata=reshape(numData,7,[])';
    hopping_r=zeros(3,nrpts);
    for i=1:nrpts
        hopping_r(1:3,i)=hoppingdata(1+nbands^2*(i-1),1:3)';
    end

    hamiltonian=zeros(nbands,nbands,nrpts);
    for i=1:nrpts
        for j=1:nbands^2
            hamiltonian(hoppingdata(j+nbands^2*(i-1),4),hoppingdata(j+nbands^2*(i-1),5),i)=...
            hoppingdata(j+nbands^2*(i-1),6)+1j*hoppingdata(j+nbands^2*(i-1),7)/degrpts(i);
        end
    end
end
