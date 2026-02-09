clc;
clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Read the structure and hop from fplo   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=read_fplo("RhG-15s-fplo");
% g=roate_geometry([1,0,0],g);
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Write the wannier90_hr.dat from fplo   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filename="data/WTe2/fplo/wannier90_hr.dat";
% g.a=g.a*0.529177249;
% g.b=2*pi*inv(g.a');
% g.wpos=g.wpos*0.529177249;
% 
% for i=1:size(g.wpos,1)
%     if g.wpos(i,1)<0
%         g.wpos(i,:)=g.wpos(i,:)+g.a(1,:);
%     end
%     if g.wpos(i,2)<0
%         g.wpos(i,:)=g.wpos(i,:)+g.a(2,:);
%     end
%     if g.wpos(i,3)<0
%         g.wpos(i,:)=g.wpos(i,:)+g.a(3,:);
%     end
% end
% 
% Electric_field_in_evpA=0.0000;
% g=add_elec(g,Electric_field_in_evpA);
% g.atoms=g.wpos*inv(g.a);
%%
% filename="data/Graphene/15s/fplo/encut_25/wannier90_formula/wannier90_hr.dat";
% MTB.write_hr(g,filename)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Read the structure and hop from wannier90   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
g=get_g_from_wannier_encut20("Rhg-15s-wannier");
Electric_field_in_evpA=0.00100;
g=add_elec(g,Electric_field_in_evpA);
tic;
% efermi=get_ef(g);
toc;
% clc;
% clear;
% g_fplo=get_g_from_wannier_encut20("Rhg-15s-wannier");
% Electric_field_in_evpA=0.0000;
% g_fplo=add_elec(g_fplo,Electric_field_in_evpA);
% tic;
% efermi=get_ef(g_fplo);
% toc;
% efermi=0;
%%
kpoint=[1/3,1/3,0]
[Energy,Psik,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%               Check the 2pi periodictivity        %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nbands,~,nrpts]=size(g.ham);
kpoint=[0.0,0.0,0.0]
[Energy,Psik,hk1]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);
kpoint=[1.0,-0.0,-0.0]
[Energy,Psik,hk2]=MTB.ham.get_parity_singleK(g.ham,g.hopr,nbands,nrpts,kpoint,g.a,g.b);

c=hk1-hk2;
max(c,[],"all")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Calculate Bulk Band structure          %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[nbands,~,nrpts]=size(g.ham);
% labels={'\Gamma','K','M','\Gamma'}; % labels for k
% hkpoints={[0.0,0.0,0.0],...
%           [1/3,1/3,0.0],...
%           [0.0,0.5,0.0],...
%           [0.0,0.0,0.0]...
%           };% hkpoints-high symmetry k points
labels={'\Gamma','K','M'}; % labels for k
hkpoints={[1/3,1/3,0.0]*0.9,...
          [1/3,1/3,0.0],...
          [1/3,1/3,0.0]+([0.0,0.5,0.0]-[1/3,1/3,0.0])*0.2,...
          };% hkpoints-high symmetry k points
nk=251;

efermi=0;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
% [Energy,kpath,kindex]=MTB.ham.get_bulk_bands_atom_gauge(g.ham,g.hopr,g.wpos,nbands,nrpts,hkpoints,nk,g.a,g.b);

MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"SrSnO-bulk")
hold on;
plot(kpath,Energy(15,:)-efermi,'Color','blue','LineWidth',2);
plot(kpath,Energy(16,:)-efermi,"Color",'red','LineWidth',2);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Calculate Band Structure in 2D             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
knum=101;
kxline=[1/3-0.04,1/3+0.04];
kyline=[1/3-0.04,1/3+0.04];
[Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,knum);
[~,Enk]=MTB.ham.get_bulk_plane_bands(g,Kx,Ky,Kz);
Enk=Enk-efermi;
%%
figure()
% surf(Enk(:,:,14),'EdgeColor','none')
surf(Enk(:,:,15),'EdgeColor','none')
hold on;
surf(Enk(:,:,16),'EdgeColor','none')
% surf(Enk(:,:,17),'EdgeColor','none')
% zlim([-0.05,0.05])
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                Check the 2pi periodictivity        %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
E_list=linspace(0,3e-3,1201);% 4100
Enum=1001; %1001
Dos_all = zeros(length(E_list), Enum);
TDos_all = zeros(length(E_list), Enum);
parfor E_idx=1:length(E_list)
    %get the model
    fprintf('Current working on E_idx:%6f , E_idx: %d',E_list(E_idx),E_idx)
    tic;
    g=get_g_from_wannier("Rhg-15s-wannier");
    Electric_field_in_evpA=E_list(E_idx);
    g=add_elec(g,Electric_field_in_evpA);
    efermi=get_ef(g);
    %get the eigenvalue on plane
    knum=501; % 1001
    kxline=[2/3-0.04,2/3+0.04];
    kyline=[1/3-0.04,1/3+0.04];
    [Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,knum);
    [~,Enk]=MTB.ham.get_bulk_plane_bands(g,Kx,Ky,Kz);
    Enk=Enk-efermi;
    %get the DOS and TDos
    plottap=2;
    Emin=-0.12;
    Emax=0.12;
    eps=(Emax-Emin)/Enum;
    Nband=size(Enk,3);
    [Eaxis,Dos,TDos]=MTB.ham.get_dos(Enk,eps,Enum,Emin,Emax,knum,plottap);
    Dos=Dos*(Emax-Emin)/Enum;
    Dos_all(E_idx,:) = Dos(:).';
    TDos_all(E_idx,:) = TDos(:).';
    toc;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                        Functions                  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function fplo=roate_geometry(MillerIndices,g)
        Umatrix = g.MillerIndicestoumatrix(MillerIndices);
        Urot = g.surfab;
        fplo=g;
end

function fplo=read_fplo(name)
        g = MTB.geometry(name);
        % g = MTB.read_poscar(g,"data/Graphene/15s/fplo/POSCAR");
        % pos=textread("data/Graphene/15s/fplo/wpos");
        % g.wpos=pos;
        % ham=textread("data/Graphene/15s/fplo/mydata-p1");
        g = MTB.read_poscar(g,"data/Graphene/15s/fplo/encut_25/POSCAR");
        pos=textread("data/Graphene/15s/fplo/encut_25/wpos");
        g.wpos=pos;
        ham=textread("data/Graphene/15s/fplo/encut_25/mydata-p1");
        tol = 1e-9;  % 容差，避免浮点误差
        is_int1 = abs(ham(:,1) - round(ham(:,1))) < tol;
        is_int2 = abs(ham(:,2) - round(ham(:,2))) < tol;
        orbital_index_logical = is_int1 & is_int2 & all(ham(:,3:5) == 0, 2);

        % orbital_index_logical= ham(:,3)==0&ham(:,4)==0&ham(:,5)==0;
        orbital_index=find(orbital_index_logical);
        orbital=ham(orbital_index_logical,1:2);
        orbital_num=sqrt(size(orbital_index,1));

        ham2=ham;
        tic;
        for i=1:size(orbital_index,1)-1
            fprintf("%d\n",i)
            if (orbital_index(i+1)-orbital_index(i))>1
                ham2(orbital_index(i)+1:orbital_index(i+1)-1,1:3)=ham2(orbital_index(i)+1:orbital_index(i+1)-1,1:3)+(g.wpos(orbital(i,1),:)-g.wpos(orbital(i,2),:));
                ham2(orbital_index(i)+1:orbital_index(i+1)-1,1:3)=ham2(orbital_index(i)+1:orbital_index(i+1)-1,1:3)/g.a;
            end
        end
        i=size(orbital_index,1);
        ham2(orbital_index(i)+1:end,1:3)=ham2(orbital_index(i)+1:end,1:3)+(g.wpos(orbital(i,1),:)-g.wpos(orbital(i,2),:));
        ham2(orbital_index(i)+1:end,1:3)=ham2(orbital_index(i)+1:end,1:3)/g.a;
        toc;

        %
        tic;
        dim=max(ham2(~orbital_index_logical,1:3))-min(ham2(~orbital_index_logical,1:3))+1;
        hambac=zeros(orbital_num,orbital_num,round(prod(dim)));
        [x,y,z]=meshgrid(round(min(ham2(~orbital_index_logical,1))):round(max(ham2(~orbital_index_logical,1))), ...
                         round(min(ham2(~orbital_index_logical,2))):round(max(ham2(~orbital_index_logical,2))), ...
                         round(min(ham2(~orbital_index_logical,3))):round(max(ham2(~orbital_index_logical,3))));
        hopr=round([x(:),y(:),z(:)]);
        toc;
        %
        tic;
        for i=1:size(orbital_index,1)-1
            fprintf("%d\n",i)
            if (orbital_index(i+1)-orbital_index(i))>1
                for j=orbital_index(i)+1:orbital_index(i+1)-1
                    index=find(ismember(hopr,round(ham2(j,1:3)),'rows'));
                    hambac(orbital(i,1),orbital(i,2),index)=ham2(j,4)+1j*ham2(j,5);
                end
            end
        end


        i=size(orbital_index,1)-1
        for j=orbital_index(i+1)+1:size(ham2,1)
            index=find(ismember(hopr,round(ham2(j,1:3)),'rows'));
            hambac(orbital(i+1,1),orbital(i+1,2),index)=ham2(j,4)+1j*ham2(j,5);
        end

        toc;
        g.ham=hambac;
        g.hopr=hopr;
        fplo=g;
end

function g=get_g_from_wannier(name)
%Create geometry object
g = MTB.geometry(name);
% Read geometry from POSCAR file
g = MTB.read_poscar(g, fullfile("data/Graphene/15s/fplo/wannier90_formula/", "POSCAR"));
% Read Wannier Hamiltonian data
[g.ham, g.hopr] = MTB.wannier.read_hr(...
    fullfile("data/Graphene/15s/fplo/wannier90_formula/", "wannier90_hr_p1.dat"), ...
    fullfile("data/Graphene/15s/fplo/wannier90_formula/", "wannier90_hr_p2.dat"));
g.wpos=g.atoms*g.a;
end

function g=get_g_from_wannier_encut40(name)
%Create geometry object
g = MTB.geometry(name);
% Read geometry from POSCAR file
g = MTB.read_poscar(g, fullfile("data/Graphene/15s/fplo/encut_40/", "POSCAR"));
% Read Wannier Hamiltonian data
[g.ham, g.hopr] = MTB.wannier.read_hr(...
    fullfile("data/Graphene/15s/fplo/encut_40/", "wannier90_hr_p1.dat"), ...
    fullfile("data/Graphene/15s/fplo/encut_40/", "wannier90_hr_p2.dat"));
g.wpos=g.atoms*g.a;
end

function g=get_g_from_wannier_encut20(name)
%Create geometry object
g = MTB.geometry(name);
% Read geometry from POSCAR file
g = MTB.read_poscar(g, fullfile("data/Graphene/15s/fplo/encut_25/wannier90_formula/", "POSCAR"));
% Read Wannier Hamiltonian data
[g.ham, g.hopr] = MTB.wannier.read_hr(...
    fullfile("data/Graphene/15s/fplo/encut_25/wannier90_formula/", "wannier90_hr_p1.dat"), ...
    fullfile("data/Graphene/15s/fplo/encut_25/wannier90_formula/", "wannier90_hr_p2.dat"));
g.wpos=g.atoms*g.a;
end

function writeoutput(filename,list)
    file=fopen(filename,'a+');
    cloumns=size(list,2);
    raws=size(list,1);
    % fprintf(file,'Time on %s\n',datetime('today'));
    % fprintf(file,'raws %d cloums %d\n',raws,cloumns);
    for i=1:raws
        for j=1:cloumns
            fprintf(file,'%12.6f \t',list(i,j));
        end
        fprintf(file,'\n');
    end
    fprintf(file,'\n');
 end


 function write_hr(g,filename)
    % open a file to write wannier hoppings
    fileID=fopen(filename,'w');
    
    [numBands,~,numRpts]=size(g.ham);

    % Get the current date and time
    currentTime = datetime('now');
    formattedTime = datestr(currentTime, 'mm/dd/yyyy at HH:MM:SS');

    % Write the header with the current data and time to the file
    fprintf(fileID,' write on %s\n', formattedTime);
    % Write the orbital numbers to the file
    fprintf(fileID,'\t %d\n', numBands);
    % Write the sites numbers to the file
    fprintf(fileID,'\t %d\n', numRpts);
    
    % numbers per line
    numsPerLine = 15;
    
    % get the site numbers of hoppings
    len=size(g.hopr,1);

    degeneracy = ones(1,len);

    % Write the data into the file with 15 numbers per line
    for i = 1:numsPerLine:len
        % Determine the index range of data for the current line
        endIdx = min(i+numsPerLine-1,len);
        fprintf(fileID, '%5d', degeneracy(i:endIdx));
        fprintf(fileID, '\n');
    end

    for i = 1:numRpts
        for j = 1:numBands
            for k = 1:numBands
                fprintf(fileID, '%5d %5d %5d %5d %5d %12.6f %12.6f\n',g.hopr(i,:),j,k,real(g.ham(j,k,i)),imag(g.ham(j,k,i)));
            end
        end
    end

    
 end

function obj=add_elec(obj,Electric_field_in_evpA)
    dim_H=size(obj.ham,1);
    minrz=min(obj.wpos(:,3));
    maxrz=max(obj.wpos(:,3));
    rz=(minrz+maxrz)/2.0;
    obj.wpos(:,3)=obj.wpos(:,3)-rz;
    ham_index=find(ismember(obj.hopr,[0,0,0],'rows'));
    for i = 1:dim_H
        obj.ham(i,i,ham_index)=obj.ham(i,i,ham_index)+obj.wpos(i,3)*Electric_field_in_evpA;
        obj.wpos(i,3)*Electric_field_in_evpA;
    end 
end

function ef=get_ef(g)
    knum=301;%501
    kxline=[-0.5,0.5];
    kyline=[-0.5,0.5];
    u=0.5;
    [Kx,Ky,Kz] = g.get_Bulk2Dkmesh(kxline,kyline,knum);
    [~,Enk]=MTB.ham.get_bulk_plane_bands(g,Kx,Ky,Kz);
    ef=calculate_ef(Enk(:), u);
end

function efermi = calculate_ef(Enk, u)
    % 计算费米能级，基于填充因子 u
    % 输入:
    % Enk: 本征值矩阵，维度 (knum^2, nbands)
    % u: 填充因子 (0 <= u <= 1)
    %
    % 输出:
    % efermi: 费米能级

    % 将能量值展平并取前 u*N 个最低能量值的最大值
    total_states = numel(Enk);                % 总的能量态数
    occupied_states = ceil(total_states * u); % 填充的态数
    efermi = max(mink(Enk(:), occupied_states));
end