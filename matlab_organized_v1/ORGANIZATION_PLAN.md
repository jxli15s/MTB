# 项目整理计划

## 📊 当前状态

### 已完成
✅ 创建新项目文件夹 `matlab_organized_v1/`
✅ 复制所有 .m 文件（148个）
✅ 复制 +MTB 工具箱
✅ 复制 +tbHFMF 模块
✅ 创建 data 符号链接

---

## 🎯 整理方案

### 方案A：按材料体系分类（推荐）

```
matlab_organized_v1/
├── +MTB/                          # 保持不变
├── +tbHFMF/                       # 保持不变
├── README.md
├── ORGANIZATION_PLAN.md
│
├── projects/                      # 新建：所有项目文件
│   ├── TaIrTe4/                  # TaIrTe4 相关（20+文件）
│   │   ├── basic/                # 基础计算
│   │   ├── quantum_geometry/     # 量子几何
│   │   ├── hfmf/                 # HFMF 计算
│   │   └── README.md
│   │
│   ├── Graphene/                 # 石墨烯系列
│   │   ├── bilayer/              # 双层
│   │   ├── rhombohedral/         # 菱形堆叠
│   │   ├── twisted/              # 扭转
│   │   └── README.md
│   │
│   ├── TMDs/                     # 过渡金属硫化物
│   │   ├── WS2/
│   │   ├── MoS2/
│   │   ├── NbSe2/
│   │   └── README.md
│   │
│   ├── Topological/              # 拓扑材料
│   │   ├── MnBiTe/
│   │   ├── SrSnO/
│   │   ├── TaRhTe4/
│   │   └── README.md
│   │
│   ├── Models/                   # 理论模型
│   │   ├── Haldane/
│   │   ├── Weyl/
│   │   ├── Kane_Mele/
│   │   └── README.md
│   │
│   └── HFMF_Applications/        # HFMF 应用项目
│       ├── BLG_valley/
│       ├── moire_lattice/
│       └── README.md
│
├── demos/                        # 示例和教程
│   ├── demo_basic_bands.m
│   ├── demo_berry_curvature.m
│   ├── demo_hfmf.m
│   └── README.md
│
├── utilities/                    # 辅助函数
│   ├── visualization/
│   ├── data_processing/
│   └── README.md
│
└── data/                         # 符号链接（保持）
```

### 方案B：按功能分类

```
matlab_organized_v1/
├── +MTB/
├── +tbHFMF/
│
├── band_structure/               # 能带结构计算
├── topology/                     # 拓扑性质
├── quantum_geometry/             # 量子几何
├── transport/                    # 输运性质
├── hfmf/                        # HFMF 计算
├── spectroscopy/                # 谱学
└── demos/
```

---

## 📝 文件分类清单

### TaIrTe4 系列（~25个文件）
- `TaIrTe4_all_data.m`
- `TaIrTe4_quantum_metric_dipole_*.m`
- `TaIrTe4_2d*.m`
- `TaIrTe4_8band*.m`
- `TaIrTe4_hfmf*.m`
- `TaIrTe4_dos.m`
- 等...

### 石墨烯系列（~15个文件）
- `BLG_MIC_derive_skew_flavors.m`
- `rhombohedral*.m`
- `Rhg_*.m`
- `Gra.m`
- `Hexagonal.m`
- `Moire.m`

### TMDs（~10个文件）
- `WS2_*.m`
- `NbSe2_*.m`
- `MoS2_*.m`（通过 fitExciton）

### 拓扑材料（~10个文件）
- `MnBiTe_*.m`
- `SrSnO_*.m`
- `TaRhTe4_*.m`

### 模型（~5个文件）
- `Haldane*.m`
- `Weyl.m`
- `TIT_TSC.m`

### HFMF 应用（~10个文件）
- `HFMF_*.m`
- `lattice_HFMF_*.m`
- `demo_tbhfmf_*.m`

### Demo/Test（~15个文件）
- `demo_*.m`
- `test*.m`
- `fig*.m`

### 辅助工具（~10个文件）
- `fitExciton*.m`
- `spectral*.m`
- `BZ_patch.m`
- `write_*.m`
- `read_*.m`

### 其他/杂项（~20个文件）
- 待分类的单独研究

---

## 🔧 整理步骤

### 第1步：清理废弃文件
- [ ] 删除所有 .asv 文件
- [ ] 识别并删除重复文件
- [ ] 删除测试文件（test22.m等）

### 第2步：创建目录结构
- [ ] 创建 projects/ 文件夹
- [ ] 创建各材料子文件夹
- [ ] 创建 demos/ 文件夹
- [ ] 创建 utilities/ 文件夹

### 第3步：移动文件
- [ ] 移动 TaIrTe4 系列
- [ ] 移动石墨烯系列
- [ ] 移动 TMDs 系列
- [ ] 移动拓扑材料
- [ ] 移动模型文件
- [ ] 移动 HFMF 应用
- [ ] 移动 demo 文件

### 第4步：文档化
- [ ] 为每个材料系统创建 README.md
- [ ] 添加使用说明
- [ ] 列出文件清单
- [ ] 说明物理背景

### 第5步：代码优化
- [ ] 统一命名规范
- [ ] 添加函数注释
- [ ] 检查路径依赖

---

## 🤔 需要决策的问题

### 1. 选择哪个整理方案？
- **方案A（按材料）**：优点是清晰，缺点是跨材料的方法难归类
- **方案B（按功能）**：优点是方法论清晰，缺点是材料分散

**建议：** 方案A，因为你的研究以材料为主线

### 2. 如何处理多材料共用的函数？
- 保留在顶层
- 创建 `shared/` 或 `utilities/` 文件夹

### 3. Demo 文件如何组织？
- 按材料分到各自文件夹
- 统一放在 `demos/`

**建议：** 统一放 `demos/`，方便新用户

### 4. 是否需要重命名文件？
- 保持原名（兼容性好）
- 统一命名（更清晰）

**建议：** 先保持原名，在各文件夹的 README 中说明

---

## 📅 实施时间表

### 立即执行（今天）
- 删除废弃文件
- 创建目录结构
- 移动 TaIrTe4 文件（主要工作）

### 短期（本周）
- 移动其他材料文件
- 创建各 README

### 中期（本月）
- 代码优化
- 添加注释
- 统一规范

---

## ✅ 执行确认

请确认以下问题后开始整理：

1. **整理方案**：使用方案A（按材料体系）？
2. **文件处理**：删除 .asv 和测试文件？
3. **命名规范**：保持原名还是重命名？
4. **立即执行**：现在开始整理？

---

**创建日期：** 2025-02-09
