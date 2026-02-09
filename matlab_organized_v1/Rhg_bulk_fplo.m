clc;
clear;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Read the structure and hop from fplo   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% g=read_fplo("RhG-fplo-bulk");
% g=roate_geometry([1,0,0],g);
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Write the wannier90_hr.dat from fplo   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% filename="data/WTe2/fplo/wannier90_hr.dat";
% g.a=g.a*0.529177249;
% g.b=2*pi*inv(g.a');
% g.wpos=g.wpos*0.529177249;

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

% Electric_field_in_evpA=0.0000;
% g=add_elec(g,Electric_field_in_evpA);
% g.atoms=g.wpos*inv(g.a);
% %%
% filename="data/Graphene/bulk/fplo/wannier_formula/wannier90_hr.dat";
% MTB.write_hr(g,filename)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       Read the structure and hop from wannier90   %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
g=get_g_from_wannier_bulk("Rhg-bulk-wannier");
Electric_field_in_evpA=0.0000;
g=add_elec(g,Electric_field_in_evpA);
% tic;
% efermi=get_ef(g);
% toc;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Plot the BZ and get the nodal line in BZ    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ===== 计算第一布里渊区
[bzV, bzRidges, bzFacets] = g.get_brillouin_zone_3d;

% 画一下
figure('Color','white'); hold on; axis equal vis3d; grid off; box on;view(35,25);
% for i = 1:numel(bzFacets)
%     F = bzFacets{i};
%     if size(F,1) >= 3
%         % Fi = convhulln(F);
%         Fi = convhulln(F, {'QJ'});
%         patch('Vertices',F,'Faces',Fi,'FaceColor','#b286d1',...
%               'FaceAlpha',0.05,'EdgeColor','none');
%     end
% end

for i = 1:numel(bzRidges)
    R = bzRidges{i};
    plot3([R(:,1); R(1,1)], [R(:,2); R(1,2)], [R(:,3); R(1,3)], 'k-', 'LineWidth',1.2);
end

plot3(bzV(:,1), bzV(:,2), bzV(:,3), 'ko', 'MarkerSize',4, 'MarkerFaceColor','k');
xlabel('k_x'); ylabel('k_y'); zlabel('k_z');
title('Helical nodal line in Rhg');
hold on;

o=[0,0,0];
b1=g.b(1,:)/2;
b2=g.b(2,:)/2;
b3=g.b(3,:)/2;
lw=1;
cols={'#1f77b4','#2ca02c','#d62728'};
cols={'black','black','black'};
quiver3(o(1), o(2), o(3), b1(1), b1(2), b1(3), 0, 'LineWidth', lw, 'Color', cols{1}, 'MaxHeadSize', 0.25);
    quiver3(o(1), o(2), o(3), b2(1), b2(2), b2(3), 0, 'LineWidth', lw, 'Color', cols{2}, 'MaxHeadSize', 0.25);
    quiver3(o(1), o(2), o(3), b3(1), b3(2), b3(3), 0, 'LineWidth', lw, 'Color', cols{3}, 'MaxHeadSize', 0.25);

% read the nodal line points from nodes.dat
hold on;
data=textread("/Users/jxli/work/dft/fplo/gra/bulk/findnodes/Nodes.dat");
kx=data(:,1);
ky=data(:,2);
kz=data(:,3);
% k1=data(:,6);
% k2=data(:,7);
% k3=data(:,8);
e=data(:,5)+0.338;
% % kx=k1*g.b(1,1)+k2*g.b(2,1)
% % ky=
% % kz=
% k=[k1,k2,k3]*g.b;
% kx=k(:,1)
% ky=k(:,2)
% kz=k(:,3)
% gap=data(:,4);
% figure('Color','White')
% sc = scatter3(kx, ky, kz, 36, E, 'filled');   % 第5个参数就是颜色
% scatter3(kx, ky, kz, 36,'filled','MarkerFaceColor','#036EB8');
% scatter3(kx+g.b(1,1), ky+g.b(1,2), kz+g.b(1,3), 18, 'filled');
% scatter3(kx-g.b(1,1), ky-g.b(1,2), kz-g.b(1,3), 18, 'filled');
% scatter3(kx+g.b(2,1), ky+g.b(2,2), kz+g.b(2,3), 18, 'filled');
% scatter3(kx-g.b(2,1), ky-g.b(2,2), kz-g.b(2,3), 18, 'filled');
% scatter3(kx+g.b(3,1), ky+g.b(3,2), kz+g.b(3,3), 18, 'filled');
% scatter3(kx-g.b(3,1), ky-g.b(3,2), kz-g.b(3,3), 18, 'filled');

%
scatter3(kx, ky, kz, 18, e,'filled');
scatter3(kx+g.b(1,1), ky+g.b(1,2), kz+g.b(1,3), 18,e, 'filled');
scatter3(kx-g.b(1,1), ky-g.b(1,2), kz-g.b(1,3), 18,e, 'filled');
scatter3(kx+g.b(2,1), ky+g.b(2,2), kz+g.b(2,3), 18,e, 'filled');
scatter3(kx-g.b(2,1), ky-g.b(2,2), kz-g.b(2,3), 18,e, 'filled');
scatter3(kx+g.b(3,1), ky+g.b(3,2), kz+g.b(3,3), 18,e, 'filled');
scatter3(kx-g.b(3,1), ky-g.b(3,2), kz-g.b(3,3), 18,e, 'filled');
colorbar
% colormap(slanCM('RdBu'));
colormap(slanCM('Blues'));
colormap(flipud(colormap));
shading interp;
% scatter3(kx, ky, kz, 18, e,'filled');
axis equal;
xlim([-2,2])
ylim([-2,2])
zlim([-2,2])
% clim([-0.0001,1])
% xlim([0.5,1.25])
% ylim([1.2,1.8])

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
% labels={'\Gamma','T','H_2','H_0','L','\Gamma','S'}; % labels for k
% hkpoints={[0,0,0.0],...
%           [0.5,0.5,0.5],...
%           [0.8060155999,0.1939844001,0.5],...
%           [0.5000000000, -0.1939844001, 0.1939844001],...
%           [0.5000000000,0.0000000000,0.0000000000],...
%           [0,0,0.0],...
%           [0.3469922000, -0.3469922000, 0.0000000000]};% hkpoints-high symmetry k points
labels={'\Gamma','T','H_0','L','\Gamma','S_0','S_2','F','\Gamma'}; % labels for k
hkpoints={[0,0,0.0],...
          [0.5,0.5,0.5],...
          [0.5000000000, -0.1939844001, 0.1939844001],...
          [0.5000000000,0.0000000000,0.0000000000],...
          [0,0,0.0],...
          [1-0.3469922000 ,0.3469922000,0],... %[0.3469922000, -0.3469922000, 0.0000000000],...[0, 0.3469922000 , -0.3469922000]
          [0.5000000000, -0.1939844001, 0.1939844001],...
          [0.5000000000,   0.0000000000,   0.5000000000],...
          [0,0,0.0]};% hkpoints-high symmetry k points
% labels={'L','H_0','S_0','\Gamma'}; % labels for k
% hkpoints={[0.5,0,0.0],...
%           [0.5000000000, -0.1939844001, 0.1939844001],...
%           [0.3469922000 ,0,-0.3469922000],... %[0.3469922000, -0.3469922000, 0.0000000000],...[0, 0.3469922000 , -0.3469922000]
%           [0,0,0.0],... 
%           };% hkpoints-high symmetry k points

nk=301;

efermi=0;
[Energy,kpath,kindex]=MTB.ham.get_bulk_bands(g.ham,g.hopr,nbands,nrpts,hkpoints,nk,g.a,g.b);
% [Energy,kpath,kindex]=MTB.ham.get_bulk_bands_atom_gauge(g.ham,g.hopr,g.wpos,nbands,nrpts,hkpoints,nk,g.a,g.b);

MTB.plot.plot_bands(Energy,nbands,efermi,kpath,labels,kindex,"SrSnO-bulk")
hold on;
plot(kpath,Energy(4,:)-efermi,'Color','blue','LineWidth',2);
plot(kpath,Energy(5,:)-efermi,"Color",'red','LineWidth',2);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                    Calculate slab                 %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
g=get_g_from_wannier_bulk("Rhg-bulk-wannier");
Electric_field_in_evpA=0.0000;
g=add_elec(g,Electric_field_in_evpA);

MillerIndices=[1,1,1];
Umatrix = g.MillerIndicestoumatrix(MillerIndices);
Urot = g.surfab;
[nbands,~,nrpts]=size(g.ham); %nbands-number of bands; nrpts-number of r points
%
labels={'\Gamma','K','M'}; % labels for k
hkpoints={[1/3,2/3]*0.9,...
          [1/3,2/3],...
          [1/3,2/3]+([0.0,0.5000000000]-[1/3,2/3])*0.2};% hkpoints-high symmetry k points
nk=51;
nslab=8;
vb_idx=4*(nslab-1);
%
[Energy,kpath,kindex]=MTB.ham.get_slab_bands(g.ham,g.hopr2,nslab,nbands,nrpts,hkpoints,nk,g.a2,g.b2);
efermi=Energy(vb_idx+1,nk)/2+Energy(vb_idx,nk)/2;
MTB.plot.plot_bands(Energy,nbands*nslab,efermi,kpath,labels,kindex,"WS2-slab")
hold on;
plot(kpath,Energy(vb_idx,:)-efermi,'Color','blue','LineWidth',2);
plot(kpath,Energy(vb_idx+1,:)-efermi,"Color",'red','LineWidth',2);
ylim([-0.4 0.4])
%%
kxline=[1/3-0.01,1/3+0.01];
kyline=[2/3-0.01,2/3+0.01];
knum=201;
[Kx,Ky] = get_Slab2Dkmesh(g,kxline,kyline,knum);
[~,Enk]=MTB.ham.get_slab_plane_bands(g,Kx,Ky,nslab);

%%
figure('Color','White')
Ecb=Enk(:,:,vb_idx+1); % for double layer n=3
Evb=Enk(:,:,vb_idx);
surf(Kx,Ky,Ecb-efermi)
hold on;
surf(Kx,Ky,Evb-efermi)
colormap(slanCM('RdBu'))
shading interp
zlim([-0.025,0.04])
%
figure('Color','White')
pcolor(Kx,Ky,Ecb-Evb)
colormap(slanCM('RdBu'))
shading interp
clim([0,0.001])
%%
data=textread("/Users/jxli/work/dft/fplo/gra/bulk/findnodes/Nodes.dat");
kx=data(:,1);
ky=data(:,2);
kz=data(:,3);
% gap=data(:,4);
figure('Color','White')
scatter3(kx, ky, kz, 36, 'filled');
%%


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
        g = MTB.read_poscar(g,"data/Graphene/bulk/fplo/POSCAR");
        pos=textread("data/Graphene/bulk/fplo/wpos");
        g.wpos=pos;
        ham=textread("data/Graphene/bulk/fplo/mydata-p1");
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

function g=get_g_from_wannier_bulk(name)
%Create geometry object
g = MTB.geometry(name);
% Read geometry from POSCAR file
g = MTB.read_poscar(g, fullfile("data/Graphene/bulk/fplo/wannier_formula/", "POSCAR"));
% Read Wannier Hamiltonian data
[g.ham, g.hopr] = MTB.wannier.read_hr(...
    fullfile("data/Graphene/bulk/fplo/wannier_formula/", "wannier90_hr_p1.dat"), ...
    fullfile("data/Graphene/bulk/fplo/wannier_formula/", "wannier90_hr_p2.dat"));
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