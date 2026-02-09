# MTB (MATLAB Tight-Binding) 工具箱说明文档

## 项目概述

MTB 是一个基于 MATLAB 的凝聚态物理紧束缚模型（Tight-Binding Model）计算工具箱，主要用于计算和分析固体材料的电子结构、拓扑性质和量子几何性质。

**主要功能：**
- 紧束缚模型的能带结构计算（体材料和表面态）
- 拓扑性质分析（Berry曲率、陈数、Wilson loop、Z2不变量）
- 量子几何性质（量子度规、量子几何张量）
- 超导系统的BdG（Bogoliubov-de Gennes）形式计算
- 态密度、反常霍尔电导等输运性质
- 与Wannier90接口，支持第一性原理+紧束缚计算

**适用领域：**
- 拓扑绝缘体和拓扑半金属研究
- 超导材料和拓扑超导体
- 二维材料（石墨烯、过渡金属硫化物等）
- 表面态和边界态研究
- 电子输运性质计算

---

## 目录结构

```
+MTB/
├── geometry.m              # 核心类：晶格几何和哈密顿量管理
├── read_hr.m              # 读取Wannier90的hr.dat文件
├── write_hr.m             # 写出Wannier90格式的hr.dat文件
├── write_wb_hr.m          # 写出WannierBerri格式的hr文件
├── read_poscar.m          # 读取VASP的POSCAR文件
│
├── +ham/                  # 哈密顿量计算和物理性质模块
│   ├── 能带计算
│   │   ├── get_bulk_bands.m           # 体材料能带
│   │   ├── get_slab_bands.m           # 表面/薄膜能带
│   │   ├── get_bulk_plane_bands.m     # 二维平面能带
│   │   ├── get_kline_bands.m          # k路径能带
│   │   └── get_bulk_bands_sparse.m    # 稀疏矩阵能带计算
│   │
│   ├── 超导BdG形式
│   │   ├── get_bulk_bands_BdG.m       # 体材料BdG能带
│   │   ├── get_slab_bands_BdG.m       # 表面BdG能带
│   │   ├── get_bulk_bands_BdG_pwave.m # p波超导
│   │   └── get_bulk_BdG_pwave_wilsonloop.m
│   │
│   ├── 拓扑性质
│   │   ├── get_Berry_curvature.m      # Berry曲率（loop方法）
│   │   ├── get_Berrycurvature_dip.m   # Berry曲率（偶极方法）
│   │   ├── get_wilsonloop.m           # Wilson loop（陈数/拓扑不变量）
│   │   ├── get_wilsonloop_xy.m        # 2D Wilson loop
│   │   ├── get_wilsonloop_mirror.m    # 镜面对称Wilson loop
│   │   ├── get_parity_singleK.m       # 宇称计算（Z2不变量）
│   │   └── get_bcd.m                  # Berry曲率偶极矩
│   │
│   ├── 量子几何
│   │   ├── quantum_geometry_general_k.m      # 单k点量子几何
│   │   ├── quantum_geometry_general_plane.m  # 平面量子几何
│   │   ├── get_quan_metric.m                 # 量子度规
│   │   ├── get_sigma_quantum_metric_dipole.m # 量子度规偶极
│   │   └── get_D_metric_kresolved_from_gv.m  # 从群速度计算D度规
│   │
│   ├── 态密度和输运
│   │   ├── get_dos.m                  # 态密度
│   │   ├── get_dos_FermiDirac.m       # 费米-狄拉克分布态密度
│   │   └── get_ahc.m                  # 反常霍尔电导
│   │
│   ├── 表面态和格林函数
│   │   ├── get_surfstates.m           # 表面态
│   │   ├── get_surfgreen.m            # 表面格林函数
│   │   └── get_anc.m                  # 迭代格林函数
│   │
│   ├── 超胞和能带展开
│   │   ├── get_supercell.m            # 构建超胞
│   │   ├── get_supercell_wannier.m    # Wannier超胞
│   │   └── get_bulk_unfolding_bands.m # 能带展开
│   │
│   └── 辅助函数
│       ├── get_slab_hk.m              # 表面哈密顿量
│       ├── get_slab_h00.m             # 表面onsite哈密顿量
│       └── cal_kmesh.m                # k网格生成
│
├── +plot/                 # 绘图模块
│   ├── plot_bands.m                   # 能带图
│   └── plot_bands_Electric.m          # 含电场的能带图
│
└── +wannier/              # Wannier90接口模块
    ├── read_hr.m                      # 读取hr.dat
    ├── write_hr.m                     # 写入hr.dat
    └── read_poscar.m                  # 读取POSCAR
```

---

## 核心类：geometry

`geometry` 是MTB工具箱的核心类，用于管理晶格结构和哈密顿量。

### 主要属性

```matlab
properties
    name;           % 系统名称
    a;              % 实空间晶格矢量 (3×3矩阵)
    a2;             % 二维晶格矢量 (2×2矩阵，用于表面计算)
    b;              % 倒空间基矢 (3×3矩阵)
    b2;             % 二维倒空间基矢 (2×2矩阵)
    atoms;          % 原子位置（分数坐标）
    hopr;           % 跳跃参数的实空间坐标 R
    hopr2;          % 二维跳跃坐标
    ham;            % 哈密顿量矩阵 H(R)
    wpos;           % Wannier函数中心位置
    iniham;         % 初始哈密顿量（用于恢复）
    sublattice;     % 子晶格标记
    orbnum_list;    % 每个原子的轨道数列表
    suborbidx;      % 子晶格轨道索引
    Rcart;          % 笛卡尔坐标
end
```

### 主要方法

#### 1. 哈密顿量计算
```matlab
hk = obj.get_hk(kpoint)           % 计算H(k)
[E, Psik] = obj.solve_kpoint(k)   % 求解本征值和本征态
```

#### 2. 晶格操作
```matlab
obj.surfab()                                  % 表面化：生成2D晶格
Umatrix = obj.MillerIndicestoumatrix([h,k,l]) % Miller指数转换
obj.lattice_plot()                            % 绘制晶格结构
```

#### 3. 布里渊区
```matlab
[bz_vertices, bz_ridges, bz_facets] = obj.get_brillouin_zone_3d()
[Kx, Ky] = obj.get_Slab2Dkmesh(kxline, kyline, knum)
[Kx, Ky, Kz] = obj.get_Bulk2Dkmesh(kxline, kyline, knum)
```

#### 4. 哈密顿量修改
```matlab
obj.onsite_modify(E)              % 修改在位能（能级移动）
obj.add_zeeman(E)                 % 添加塞曼场
obj.offsite_modify(hopr, ham)     % 修改跳跃项
```

---

## 主要功能模块详解

### 1. 能带结构计算

#### 体材料能带
```matlab
[Energy, kpath, kk] = MTB.ham.get_bulk_bands(...
    hamiltonian, hopping_r, nbands, nrpts, hkpoints, nk, a, b)
```
- **输入：**
  - `hamiltonian`: 哈密顿量矩阵 H(R)
  - `hopping_r`: 跳跃坐标R
  - `nbands`: 能带数
  - `nrpts`: R点数量
  - `hkpoints`: 高对称点路径，如 `{[0,0,0], [0.5,0,0], [0.5,0.5,0]}`
  - `nk`: 每段路径的k点数
  - `a`, `b`: 实空间和倒空间基矢
- **输出：**
  - `Energy`: nbands × nkpts 能量矩阵
  - `kpath`: k点路径坐标
  - `kk`: 高对称点位置

#### 表面态能带
```matlab
[Energy, kpath, kk] = MTB.ham.get_slab_bands(...
    hamiltonian, hopping_r, nslab, nbands, nrpts, hkpoints, nk, a, b)
```
- **新增参数：**
  - `nslab`: 表面层数

#### 二维平面能带
```matlab
[Enk, KX, KY, Unk] = MTB.ham.get_bulk_plane_bands(...
    hamiltonian, hopping_r, nbands, nrpts, kxline, kyline, knum, a, b)
```
- **输出：**
  - `Enk`: knum × knum × nbands 能谱
  - `KX, KY`: k网格坐标
  - `Unk`: 波函数（用于拓扑计算）

### 2. 拓扑性质计算

#### Berry曲率
```matlab
[Omega_k, KX, KY] = MTB.ham.get_Berry_curvature(...
    bandindex, Unk, KX, KY, plottap)
```
- **原理：** 使用Wilson loop方法计算孤立能带的Berry曲率
- **输出：** Ω(k) = Im[log(U₁U₂/(U₃U₄))]/(2π)

#### Wilson Loop（陈数）
```matlab
[wx, unk] = MTB.ham.get_wilsonloop(obj, knum, band1, band2)
```
- **用途：** 计算占据态的陈数，判断拓扑绝缘体
- **输出：** wx为Wilson loop本征值随ky的演化

#### 宇称计算（Z₂不变量）
```matlab
[delta, Parity] = MTB.ham.get_parity_singleK(hamiltonian, ...)
```
- **用途：** 通过时间反演不变动量点（TRIM）的宇称计算Z₂拓扑不变量

#### Berry曲率偶极矩
```matlab
[bcd_x, bcd_y] = MTB.ham.get_bcd(Omega_k, Enk, vx, vy, knum, ef)
```
- **物理意义：** 非线性霍尔效应的来源

### 3. 超导BdG形式

#### BdG哈密顿量
BdG哈密顿量形式：
```
H_BdG = [ H(k) - μ      Δ(k)    ]
        [ Δ†(k)      -H*(-k) + μ ]
```

#### 体材料BdG能带
```matlab
[Energy, kpath, kk] = MTB.ham.get_bulk_bands_BdG(...
    hamiltonian, hopping_r, nbands, nrpts, hkpoints, nk, a, b, mu, delta)
```
- **新增参数：**
  - `mu`: 化学势
  - `delta`: 超导配对强度Δ

#### 表面BdG能带（拓扑超导）
```matlab
[Energy, kpath, kk] = MTB.ham.get_slab_bands_BdG(...)
```
- **用途：** 寻找Majorana零能模

#### p波超导Wilson loop
```matlab
wx = MTB.ham.get_bulk_BdG_pwave_wilsonloop(...)
```
- **用途：** 计算手性拓扑超导的陈数

### 4. 量子几何

#### 量子度规和量子几何张量
```matlab
[gk, Qk, Fk, vk] = MTB.ham.quantum_geometry_general_plane(...
    obj, Kx, Ky, Kz, Unk, Enk, band_list, dk_list, delta)
```
- **物理量：**
  - `gk`: 量子度规 g_ij = Re[Q_ij]
  - `Qk`: 量子几何张量 Q_ij = ⟨∂ᵢu|P⊥|∂ⱼu⟩
  - `Fk`: Berry曲率张量 F_ij = 2·Im[Q_ij]
  - `vk`: 群速度 v_i = ∂E/∂k_i

#### 量子度规偶极
```matlab
[sigma_x, sigma_y] = MTB.ham.get_sigma_quantum_metric_dipole(...)
```
- **物理意义：** 与非线性光学响应相关

### 5. 态密度和输运

#### 态密度（DOS）
```matlab
[Eaxis, Dos, TDos] = MTB.ham.get_dos(En, eps, Enum, Emin, Emax, nk, plottap)
```
- **方法：** 高斯展宽
- **输出：**
  - `Dos`: 态密度
  - `TDos`: 累积态密度

#### 反常霍尔电导（AHC）
```matlab
[Eaxis, sigma] = MTB.ham.get_ahc(Omega_k, Enk, Enum, Emin, Emax, tem)
```
- **公式：** σ_xy = (e²/ℏ) ∫ f(E) Ω(k) d²k
- **参数：** `tem` 为温度（K）

### 6. 表面格林函数

#### 迭代格林函数法
```matlab
[GS_LL, GS_RR, GS_LR] = MTB.ham.get_anc(H00, H01, omega, eta)
```
- **用途：** 计算半无限系统的表面格林函数
- **应用：** LDOS、表面态、输运计算

### 7. 超胞和能带展开

#### 构建超胞
```matlab
[ham_super, hopr_super] = MTB.ham.get_supercell(ham, hopr, nbands, sc)
```
- **输入：** `sc` = [nx, ny, nz] 超胞倍数

#### 能带展开
```matlab
[EKN, A] = MTB.ham.get_bulk_unfolding_bands(...)
```
- **用途：** 将超胞能带投影回原胞布里渊区

---

## 使用示例

### 示例1：读取Wannier90数据并计算能带

```matlab
% 1. 读取wannier90_hr.dat
[ham, hopr] = MTB.read_hr('wannier90_hr.dat');

% 2. 创建geometry对象
g = MTB.geometry('MyMaterial');

% 3. 设置晶格常数（单位：埃）
g.a = [3.0, 0, 0;
       0, 3.0, 0;
       0, 0, 10.0];
g.b = inv(g.a') * 2*pi;

% 4. 赋值哈密顿量
g.ham = ham;
g.hopr = hopr;

% 5. 定义k路径
hkpoints = {[0, 0, 0],      % Gamma
            [0.5, 0, 0],    % X
            [0.5, 0.5, 0],  % M
            [0, 0, 0]};     % Gamma

% 6. 计算能带
nbands = size(ham, 1);
nrpts = size(hopr, 1);
[Energy, kpath, kk] = MTB.ham.get_bulk_bands(...
    ham, hopr, nbands, nrpts, hkpoints, 50, g.a, g.b);

% 7. 绘制能带
labels = {'\Gamma', 'X', 'M', '\Gamma'};
MTB.plot.plot_bands(Energy, nbands, 0, kpath, labels, kk, 'output', 0);
```

### 示例2：计算Berry曲率和陈数

```matlab
% 1. 计算二维k网格的能谱和波函数
knum = 100;
kxline = [0, 1];
kyline = [0, 1];
[Enk, KX, KY, Unk] = MTB.ham.get_bulk_plane_bands(...
    g.ham, g.hopr, nbands, nrpts, kxline, kyline, knum, g.a, g.b);

% 2. 计算Berry曲率（例如第2条能带）
bandindex = 2;
[Omega_k, KX, KY] = MTB.ham.get_Berry_curvature(...
    bandindex, Unk, KX, KY, 1);  % 最后的1表示绘图

% 3. 计算陈数
Chern_number = sum(Omega_k, 'all');
fprintf('Chern number of band %d: %.4f\n', bandindex, Chern_number);

% 4. 计算Wilson loop
[wx, unk] = MTB.ham.get_wilsonloop(g, knum, 1, 2);  % 对能带1-2
```

### 示例3：计算表面态

```matlab
% 1. 首先将晶格表面化（例如沿z方向切表面）
g.surfab();  % 生成a2, b2

% 2. 定义2D k路径
hkpoints_2d = {[0, 0],       % Gamma
               [0.5, 0],     % M
               [0.5, 0.5],   % K
               [0, 0]};      % Gamma

% 3. 计算表面态能带
nslab = 50;  % 表面层数
[Energy_slab, kpath, kk] = MTB.ham.get_slab_bands(...
    g.ham, g.hopr, nslab, nbands, nrpts, hkpoints_2d, 50, g.a, g.b);

% 4. 绘制表面能带
labels = {'\Gamma', 'M', 'K', '\Gamma'};
MTB.plot.plot_bands(Energy_slab, nbands*nslab, 0, kpath, labels, kk, 'slab', 0);
```

### 示例4：超导BdG计算

```matlab
% 1. 设置超导参数
mu = 0.5;      % 化学势 (eV)
delta = 0.1;   % 超导gap (eV)

% 2. 计算BdG能带
[Energy_BdG, kpath, kk] = MTB.ham.get_bulk_bands_BdG(...
    g.ham, g.hopr, nbands, nrpts, hkpoints, 50, g.a, g.b, mu, delta);

% 3. 绘制（能带数变为2*nbands）
MTB.plot.plot_bands(Energy_BdG, 2*nbands, 0, kpath, labels, kk, 'BdG', 0);

% 4. 计算表面态（寻找Majorana零能模）
[Energy_BdG_slab, kpath, kk] = MTB.ham.get_slab_bands_BdG(...
    g.ham, g.hopr, nslab, nbands, nrpts, hkpoints_2d, 50, g.a, g.b, mu, delta);
```

### 示例5：计算量子几何

```matlab
% 1. 获取k网格数据
knum = 50;
[Enk, KX, KY, Unk] = MTB.ham.get_bulk_plane_bands(...
    g.ham, g.hopr, nbands, nrpts, [0,1], [0,1], knum, g.a, g.b);

% 2. 设置微分方向
dk_list = [g.b(1,:)/knum;
           g.b(2,:)/knum];

% 3. 计算量子几何
band_list = [1, 2];  % 计算能带1和2
delta = 1e-6;        % 能隙阈值
Kz = zeros(size(KX));  % 2D系统
[gk, Qk, Fk, vk] = MTB.ham.quantum_geometry_general_plane(...
    g, KX, KY, Kz, Unk, Enk, band_list, dk_list, delta);

% 4. 提取第1条带的量子度规和Berry曲率
gxx_band1 = squeeze(gk(:,:,1,1,1));  % g_xx
gyy_band1 = squeeze(gk(:,:,2,2,1));  % g_yy
gxy_band1 = squeeze(gk(:,:,1,2,1));  % g_xy
Omega_z_band1 = squeeze(Fk(:,:,1,2,1));  % F_xy = Ω_z

% 5. 绘制
figure;
pcolor(KX, KY, gxx_band1);
shading interp;
colorbar;
title('Quantum metric g_{xx}');
```

### 示例6：态密度和反常霍尔电导

```matlab
% 1. 计算态密度
eps = 0.01;     % 高斯展宽
Enum = 1000;    % 能量点数
Emin = -3;
Emax = 3;
[Eaxis, Dos, TDos] = MTB.ham.get_dos(...
    Enk, eps, Enum, Emin, Emax, knum, 1);

% 2. 计算反常霍尔电导
tem = 300;  % 温度 (K)
[Eaxis_ahc, sigma_xy] = MTB.ham.get_ahc(...
    Omega_k, Enk, Enum, Emin, Emax, tem);

% 3. 绘制
figure;
subplot(1,2,1);
plot(Eaxis, Dos, 'LineWidth', 2);
xlabel('E (eV)'); ylabel('DOS');

subplot(1,2,2);
plot(Eaxis_ahc, sigma_xy, 'LineWidth', 2);
xlabel('E (eV)'); ylabel('\sigma_{xy} (e^2/h)');
```

### 示例7：修改哈密顿量（加电场、塞曼场等）

```matlab
% 1. 保存初始哈密顿量
g.iniham = g.ham;

% 2. 移动能级（模拟栅极电压）
E_shift = 0.1;  % eV
g.onsite_modify(E_shift);

% 3. 添加塞曼项（磁场）
% 构造塞曼矩阵（例如自旋劈裂）
nbands = size(g.ham, 1);
B = 1;  % 磁场强度（任意单位）
mu_B = 0.0578;  % Bohr磁子 (meV/T)
Zeeman = zeros(nbands, nbands);
% 假设轨道1-2是自旋上，3-4是自旋下
Zeeman(1:2, 1:2) = eye(2) * mu_B * B;
Zeeman(3:4, 3:4) = -eye(2) * mu_B * B;
g.add_zeeman(Zeeman);

% 4. 修改特定跳跃项
hopr_modify = [1, 0, 0];  % 要修改的R
ham_modify = ones(nbands, nbands) * 0.01;  % 修改量
g.offsite_modify(hopr_modify, ham_modify);

% 5. 重新计算能带
[Energy_modified, kpath, kk] = MTB.ham.get_bulk_bands(...
    g.ham, g.hopr, nbands, nrpts, hkpoints, 50, g.a, g.b);
```

---

## 与Wannier90工作流程

### 典型工作流：DFT → Wannier90 → MTB

```matlab
% 步骤1：从wannier90获取紧束缚模型
[ham, hopr] = MTB.read_hr('wannier90_hr.dat');

% 步骤2：读取晶格结构（从POSCAR或wannier90.win）
% 选项A：从POSCAR读取
lattice = MTB.read_poscar('POSCAR');

% 选项B：手动设置
g = MTB.geometry('Material');
g.a = [...];  % 从wannier90.win获取
g.b = inv(g.a') * 2*pi;

% 步骤3：赋值哈密顿量
g.ham = ham;
g.hopr = hopr;
g.iniham = ham;  % 保存初始值

% 步骤4：后续计算（能带、拓扑等）
...

% 步骤5：（可选）修改后写回hr.dat
MTB.write_hr(g, 'modified_hr.dat');
```

---

## 重要函数列表

### 能带计算（+ham）
| 函数名 | 功能 | 类型 |
|-------|------|------|
| `get_bulk_bands` | 体材料能带 | 基础 |
| `get_slab_bands` | 表面/薄膜能带 | 基础 |
| `get_bulk_plane_bands` | 二维k平面能谱+波函数 | 基础 |
| `get_bulk_plane_bands_with_Ham` | 带哈密顿量输出的平面能带 | 高级 |
| `get_kline_bands` | 任意k路径能带 | 基础 |
| `get_bulk_bands_add_electric` | 含电场的能带 | 扩展 |
| `get_bulk_bands_sparse` | 稀疏矩阵能带（大体系） | 优化 |
| `get_bulk_bands_atom_gauge` | 原子规范能带 | 高级 |
| `get_bulk_bands_full` | 完整信息能带 | 高级 |

### BdG超导（+ham）
| 函数名 | 功能 |
|-------|------|
| `get_bulk_bands_BdG` | 体BdG能带 |
| `get_slab_bands_BdG` | 表面BdG能带 |
| `get_bulk_bands_BdG_pwave` | p波超导 |
| `get_bulk_BdG_pwave_wilsonloop` | p波Wilson loop |
| `get_slab_bands_BdG_at_q` | 特定q的BdG |
| `get_bulk_bands_BdG_sparse` | 稀疏BdG |

### 拓扑性质（+ham）
| 函数名 | 功能 |
|-------|------|
| `get_Berry_curvature` | Berry曲率（loop法） |
| `get_Berrycurvature_dip` | Berry曲率（偶极法） |
| `get_Berrycurvature_cop` | Berry曲率（协变法） |
| `get_wilsonloop` | Wilson loop（陈数） |
| `get_wilsonloop_xy` | 2D Wilson loop |
| `get_wilsonloop_mirror` | 镜面Wilson loop |
| `get_parity_singleK` | 宇称（Z2不变量） |
| `get_bcd` | Berry曲率偶极矩 |
| `get_bcd2` | Berry曲率偶极矩v2 |

### 量子几何（+ham）
| 函数名 | 功能 |
|-------|------|
| `quantum_geometry_general_k` | 单k点量子几何 |
| `quantum_geometry_general_plane` | 平面量子几何 |
| `get_quan_metric` | 量子度规 |
| `get_sigma_quantum_metric_dipole` | 量子度规偶极 |
| `get_D_metric_kresolved_from_gv` | 从群速度算D度规 |
| `get_Dk_qmd_plainD_core` | D度规核心 |

### 态密度和输运（+ham）
| 函数名 | 功能 |
|-------|------|
| `get_dos` | 态密度（高斯展宽） |
| `get_dos_FermiDirac` | 费米分布态密度 |
| `get_ahc` | 反常霍尔电导 |

### 表面和格林函数（+ham）
| 函数名 | 功能 |
|-------|------|
| `get_surfstates` | 表面态 |
| `get_surfgreen` | 表面格林函数 |
| `get_anc` | 迭代格林函数 |
| `get_slab_hk` | 表面哈密顿量 |
| `get_slab_h00` | 表面onsite矩阵 |

### 超胞和展开（+ham）
| 函数名 | 功能 |
|-------|------|
| `get_supercell` | 构建超胞 |
| `get_supercell_wannier` | Wannier超胞 |
| `get_supercell_wannier_3d` | 3D Wannier超胞 |
| `get_bulk_unfolding_bands` | 能带展开 |

---

## 物理背景知识

### 1. 紧束缚模型

紧束缚哈密顿量：
```
H(k) = ∑_R H(R) e^(ik·R)
```
其中 H(R) 是实空间跳跃矩阵。

### 2. Berry曲率和陈数

Berry曲率：
```
Ω_n(k) = ∇_k × A_n(k)
```
其中 Berry联络 A = i⟨u|∇_k|u⟩

陈数：
```
C = (1/2π) ∫_BZ Ω(k) d²k
```

### 3. Wilson Loop

Wilson loop算符：
```
W(ky) = exp[i ∮ A·dk]
```
其本征值相位即为Berry相位，积分可得陈数。

### 4. 量子度规

量子度规：
```
g_ij = Re[⟨∂_i u|P_⊥|∂_j u⟩]
```
其中 P_⊥ = 1 - |u⟩⟨u| 是投影算符。

### 5. BdG形式

超导准粒子满足Bogoliubov-de Gennes方程：
```
[ H-μ     Δ   ] [ u ] = E [ u ]
[ Δ†   -H*+μ  ] [ v ]     [ v ]
```

---

## 计算技巧和注意事项

### 1. k网格选择
- **能带计算：** nk = 50-100 通常足够
- **Berry曲率：** knum = 100-200（需要较密）
- **态密度：** knum = 200-500（需要很密）

### 2. 表面计算
- **层数选择：** nslab = 30-100，需要测试收敛性
- **表面态识别：** 观察波函数局域化

### 3. BdG计算
- **化学势μ：** 需要先计算正常态能带确定
- **配对Δ：** 通常 Δ << 带宽

### 4. 并行计算
代码中使用 `parfor`，自动利用多核：
```matlab
% 检查并行池
pool = gcp('nocreate');
if isempty(pool)
    parpool('local', 8);  % 8核
end
```

### 5. 大规模计算
对于大体系，使用稀疏矩阵版本：
- `get_bulk_bands_sparse`
- `get_bulk_bands_BdG_sparse`

---

## 扩展和自定义

### 添加新的物理效应

**示例：添加自旋轨道耦合（SOC）**

```matlab
% 1. 定义SOC哈密顿量
lambda_soc = 0.1;  % SOC强度
% 假设轨道顺序为 [px↑, py↑, px↓, py↓]
H_soc = zeros(4, 4);
H_soc(1, 4) = 1i * lambda_soc;   % ⟨px↑|H_soc|py↓⟩
H_soc(2, 3) = -1i * lambda_soc;  % ⟨py↑|H_soc|px↓⟩
H_soc = H_soc + H_soc';

% 2. 添加到在位能
g.add_zeeman(H_soc);
```

**示例：添加应变效应**

```matlab
% 修改跳跃参数
epsilon = 0.01;  % 应变
for i = 1:size(g.hopr, 1)
    R = g.hopr(i, :);
    if norm(R) > 0
        % 沿x方向应变修改跳跃
        strain_factor = 1 + epsilon * abs(R(1));
        g.ham(:, :, i) = g.ham(:, :, i) * strain_factor;
    end
end
```

---

## 引用和参考

如果使用本工具箱，建议引用以下相关文献：

### 紧束缚方法
- Marzari, N. & Vanderbilt, D. Maximally localized generalized Wannier functions. *Phys. Rev. B* **56**, 12847 (1997).
- Mostofi, A. A. et al. wannier90: A tool for obtaining maximally-localised Wannier functions. *Comput. Phys. Commun.* **178**, 685-699 (2008).

### 拓扑性质
- Fukui, T., Hatsugai, Y. & Suzuki, H. Chern numbers in discretized Brillouin zone. *J. Phys. Soc. Jpn.* **74**, 1674-1677 (2005).
- Soluyanov, A. A. & Vanderbilt, D. Computing topological invariants without inversion symmetry. *Phys. Rev. B* **83**, 235401 (2011).

### 量子几何
- Resta, R. The insulating state of matter: A geometrical theory. *Eur. Phys. J. B* **79**, 121-137 (2011).
- Gao, Y. & Xiao, D. Nonreciprocal directional dichroism induced by the quantum metric dipole. *Phys. Rev. Lett.* **122**, 227402 (2019).

---

## 版本历史

- **v0.1** (2024-07-26): 初始版本，包含基础能带和拓扑计算
- **当前版本** (2025):
  - 新增量子几何模块
  - 完善BdG形式计算
  - 优化大规模计算性能
  - 添加更多绘图功能

---

## 联系方式

如有问题或建议，请联系开发者或提交Issue。

**主要贡献者：** JXLI

**项目路径：** `/Volumes/T9/work/tb/matlab/+MTB/`

---

## 附录：常见问题

### Q1: 如何判断计算的拓扑不变量是否可靠？

**A:**
1. 增加k网格密度，检查陈数是否收敛到整数
2. 检查能隙是否足够大（Δ > 0.1 eV）
3. 对比不同方法（Berry曲率积分 vs Wilson loop）

### Q2: 表面态计算层数如何选择？

**A:**
1. 从小层数（nslab=30）开始
2. 逐渐增加，观察表面态能量是否收敛
3. 通常拓扑绝缘体需要 nslab > 50

### Q3: 如何处理能带交叉点？

**A:**
1. Berry曲率在交叉点奇异，需要使用自适应网格
2. 或者使用更高阶数值导数方法
3. 量子几何计算中设置合适的 `delta` 参数（如1e-6）

### Q4: 如何加速大规模计算？

**A:**
1. 使用稀疏矩阵函数
2. 开启MATLAB并行计算
3. 只计算费米面附近的能带
4. 使用GPU加速（需要自行修改代码）

### Q5: 如何导出数据到Python？

**A:**
```matlab
% 保存为.mat文件
save('data.mat', 'Energy', 'kpath', 'Omega_k');

% Python中读取
import scipy.io
data = scipy.io.loadmat('data.mat')
Energy = data['Energy']
```

---

**文档版本：** v1.0
**最后更新：** 2025-02-09
