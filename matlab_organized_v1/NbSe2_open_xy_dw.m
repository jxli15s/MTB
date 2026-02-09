clear;
clear all;
%parpool('local',48)
% parpool("Threads",48)
% read fplo and set the atom positions
g=read_fplo("NbSe2");
g.atoms=[0,0,0];
g.wpos=g.atoms*g.a;
g.sublattice=1:size(g.atoms,1);
g.orbnum_list=[6,];
g.get_suborbidx;
% build the sqrt(3)x1 supercell
Umatrix=[2,1,0;0,1,0;0,0,1];
shift=[0.0,0,0]; %
gs=MTB.ham.get_supercell_wannier_3d_general(g,Umatrix,shift,"fill");
gs.sublattice=1:size(gs.atoms,1);
gs.wpos=[];
for i=1:length(gs.orbnum_list)
        wpos=repmat(gs.atoms(i,:)*gs.a,gs.orbnum_list(i),1);
        gs.wpos=[gs.wpos;wpos];
end

mu=0;
delta=0.03;
numEigs=12;

% Build n1xn2 open boundary
n1=101;
n2=151;
tic;
gs_xy = MTB.ham.get_xy_open_wannier(gs, n1, n2);
gs_xy.ham = gs_xy.ham{1};
toc;

%删除靠近 x=0 或 y=0 的边界原子（保留 mx 和 my 对称）
del_index = find(gs_xy.wpos(:,1) < 1e-4 | gs_xy.wpos(:,2) < 1e-4);
gs_xy.ham(del_index, :) = [];
gs_xy.ham(:, del_index) = [];
gs_xy.wpos(del_index, :) = [];
gs_xy.atoms=gs_xy.wpos(1:6:end,:)*inv(gs_xy.a);
%%
[E, V] = get_open_xy_pi_junction_2d(gs_xy, mu, delta, numEigs);
%%
save("open_xy_E_V_101_151.mat","E","V","gs_xy","-v7.3");
%%
% load("open_xy_E_V_101_151.mat")
load("data/NbSe2/open_xy/101x151/open_xy_E_V_101_151_small_mem.mat")
% load("data/NbSe2/open_xy/101x151/open_xy_E_V_301_151_small_mem.mat")
figure
plot(real(E),'ro')
%%
figure()
x=gs_xy.wpos(1:end,1);
y=gs_xy.wpos(1:end,2);
p=abs(V).^2;
p=p(1:end/2,:)+p(end/2+1:end,:);
band=5:5;
scatter3(x, y, sum(p(:, band), 2), 100, sum(p(:,band), 2), 'filled');
%%
figure()
slec_indx=find(gs_xy.wpos(:,1)<153.0 & gs_xy.wpos(:,1)>133.0);
selectx=gs_xy.wpos(slec_indx,1);
selecty=gs_xy.wpos(slec_indx,2);
selecty=reshape(selecty,6,[]);
selecty=sum(selecty,1)/6;

strength=p(slec_indx,:);
band=5:5;
strength=sum(strength(:,band),2);
strength=reshape(strength,6,[]);
strength=sum(strength,1);

plot(selecty,strength,'ro');


%%
function [E, V] = get_open_xy_pi_junction_2d(gs_xy, mu, delta, numEigs)
% 构造并求解具有π-junction的BdG Hamiltonian
% 输入参数：
%   mu      : 化学势
%   delta   : 超导 pairing 强度
%   numEigs : 求解的低能本征态数量
%
% 输出参数：
%   E : numEigs 个本征值（升序排列）
%   V : 对应的本征态（列向量）

    % 构造xy开边界体系


    % 构造 Pi junction 相位（y方向）
    ypos = gs_xy.wpos(:,2);
    ymean = mean(ypos);
    val = zeros(size(ypos));
    val(ypos < ymean - 1) = 1;
    val(ypos > ymean + 1) = -1;
    val = val(1:2:end);  % 只取一半轨道（可能是去掉自旋重复）

    % 设置 pairing 和 BdG 结构
    sigma_y = sparse([0, -1j; 1j, 0]);
    orbital = spdiags(val, 0, length(val), length(val));
    nbands = size(gs_xy.ham, 1);

    h_delta = kron(orbital, 1j * sigma_y * delta);
    h_onsite = speye(nbands) * mu;
    hk_e = gs_xy.ham - h_onsite;
    hk_h = h_onsite - gs_xy.ham.';
    hk = [hk_e, h_delta; h_delta', hk_h];
    hk = (hk + hk') / 2;  % 保证厄米
    % 对角化
    tic;
    [V, D] = eigs(hk, numEigs, 'smallestabs');
    toc;
    % 排序
    E = diag(D);
    [E, idx] = sort(E, 'ascend');
    V = V(:, idx);
end


function fplo=read_fplo(name)
        g = MTB.geometry(name);
        g = MTB.read_poscar(g,"data/NbSe2/fplo/POSCAR");
        pos=textread("data/NbSe2/fplo/wpos");
        g.wpos=pos;
        ham=textread("data/NbSe2/fplo/mydata-p1");
        orbital_index_logical= ham(:,3)==0&ham(:,4)==0&ham(:,5)==0;
        orbital_index=find(orbital_index_logical);
        orbital=ham(orbital_index_logical,1:2)
        orbital_num=sqrt(size(orbital_index,1))

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
%%
function C = construct_rotation_operator_periodic(wpos, a, b, R, k)
% wpos: Norb×3，轨道位置（笛卡尔坐标）
% a: 3×3 实空间晶格基矢
% R: 3×3 旋转矩阵
% k: 1×3 动量点
Norb = size(wpos, 1);
C = zeros(Norb);

% 平移搜索范围：-1, 0, 1
range = -2:2;
[X, Y, Z] = ndgrid(range, range, 0);  % 平移方向：xy方向±1
shifts = [X(:), Y(:), Z(:)];  % 共 3×3×1 = 9 个平移
Tlist = shifts * a;  % 将分数坐标转为笛卡尔坐标

wpos_rot = (R * wpos')';  % 所有轨道旋转后的位置
k=k*b;
for i = 1:Norb
    ri_rot = wpos_rot(i, 1:3);  % 旋转后的轨道位置

    min_dist = inf;
    j_best = -1;
    T_best = [0 0 0];

    for j = 1:Norb
        for t = 1:size(Tlist,1)
            T = Tlist(t,:);
            diff = wpos(j,:) - (ri_rot + T);
            d = norm(diff);
            if d < min_dist
                min_dist = d;
                j_best = j;
                T_best = T;
            end
        end
    end

    if min_dist > 1e-3
        error('No matching site found for orbital %d under rotation.', i);
    end


    % R=round((ri_rot-wpos(j_best,:))*inv(a));
    % phase = exp(1j * dot(k, R));
    % C(j_best, i) = phase;
    phase = exp(1j * dot(k, wpos(j_best,:) - ri_rot));
    C(j_best, i) = phase;

end
end