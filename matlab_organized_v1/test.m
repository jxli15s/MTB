%%
tic
[hamiltonian,hopping_r] = read_hr("data/MnBiTe-xu/wannier90_hr_MnBiTe.dat");
toc
%%
tic
[hamiltonian,hopping_r] = read_hr_v2("data/MnBiTe-xu/wannier90_hr_p1.dat","data/MnBiTe-xu/wannier90_hr_p2.dat");
toc
%%
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
    hopping_r=hopping_r'
    hamiltonian=zeros(nbands,nbands,nrpts);
    for i=1:nrpts
        for j=1:nbands^2
            hamiltonian(hoppingdata(j+nbands^2*(i-1),4),hoppingdata(j+nbands^2*(i-1),5),i)=...
            hoppingdata(j+nbands^2*(i-1),6)+1j*hoppingdata(j+nbands^2*(i-1),7)/degrpts(i);
        end
    end
    
    % parfor num = nrpts*nbands^2
    %     [j,i]=ind2sub([nbands^2,nrpts],num);
    %      hamiltonian(hoppingdata(j+nbands^2*(i-1),4),hoppingdata(j+nbands^2*(i-1),5),i)=...
    %      hoppingdata(j+nbands^2*(i-1),6)+1j*hoppingdata(j+nbands^2*(i-1),7)/degrpts(i);
    % end
end


function [hamiltonian,hopping_r] = read_hr_v2(filename1,filename2)
    %read wannier90_hr.dat
    %Return:
    %hamiltonian: real space hopping parameters
    %hopping_r: real space hopping coordinates
    data1=importdata(filename1);
    numData=data1.data;
    nbands=numData(1);
    nrpts=numData(2);
    degrpts=numData(3:2+nrpts);
    data2=readmatrix(filename2);
    hopping_r=data2(1:nbands*nbands:end,1:3);
    hamiltonian=data2(:,6)+1j*data2(:,7);
    degrpts=kron(degrpts,ones(nbands*nbands,1));
    hamiltonian=hamiltonian./degrpts;
    hamiltonian=reshape(hamiltonian,[nbands,nbands,nrpts]);
    % for i=1:nrpts
    %     hamiltonian(:,:,i)=hamiltonian(:,:,i)/degrpts(i);
    % end
end