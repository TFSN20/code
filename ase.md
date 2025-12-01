## 环境
### ase mace chgnet
```
pip install ase
pip install mace-torch
pip install chgnet
```
### wsl2 Ubuntu
- 临时python3.12.10
  ```
  python3 --version
  alias python3='python3.12'
  python3 --version
  # python3 -m venv ~/calc_env
  source ~/calc_env/bin/activate
  pip install torch --index-url https://download.pytorch.org/whl/cpu
  ```
- 
## 文档
### ase
- https://ase-lib.org/
### mace
- https://github.com/ACEsuit/mace
- 模型：https://github.com/ACEsuit/mace-foundations
### chgnet
- https://github.com/CederGroupHub/chgnet
## 模型选择
### 1. 无机晶体与材料科学 (Inorganic Materials)

如果你的研究对象是 **块体材料、晶体结构、缺陷或固态相变**：

*   **首选推荐 (Standard)**: `MACE-MP-0b3`
    *   **适用场景**: 大多数无机化合物（氧化物、硫化物、电池材料等）。
    *   **特点**: 基于 Materials Project 数据 (DFT PBE+U)。相比早期版本 (0a/0b)，它修复了高压下的稳定性问题，并显著改善了声子谱（Phonons）的预测精度。
*   **增强鲁棒性 (Robustness)**: `MACE-MPA-0`
    *   **适用场景**: 高温分子动力学 (MD)、高压相变或极端环境模拟。
    *   **特点**: 在 MP 训练集基础上增加了 sAlex 数据，对非平衡构型的预测更稳健。
*   **高精度需求 (High Accuracy)**: `MACE-MATPES-r2SCAN-0`
    *   **适用场景**: 对晶格常数、形成能有更高精度要求的体系。
    *   **特点**: 基于 r2SCAN (Meta-GGA) 泛函训练。通常比 PBE 泛函更准确，且推理速度相同。
*   **纯金属/合金 (Metals)**: `MACE-MATPES-PBE-0`
    *   **适用场景**: 不适合加 Hubbard U 的体系（如纯金属、金属间化合物）。
    *   **特点**: 纯 PBE 泛函，去除了 Materials Project 默认的 +U 修正。
*   **声子性质研究 (Phonons)**: `MACE-OMAT-0`
    *   **特点**: 专门标注了 "Excellent phonons"，在晶格动力学性质预测上表现优异。

---

### 2. 有机分子与药物化学 (Organic & Molecules)

如果你的研究对象是 **小分子、有机液体、生物大分子**：

*   **标准有机分子**: `MACE-OFF23` (Small/Medium/Large)
    *   **适用范围**: 仅支持 10 种元素 (H, C, N, O, F, P, S, Cl, Br, I)。
    *   **场景**: 药物分子筛选、有机化学反应（中性体系）。
    *   **理论水平**: wB97M+D3 (高精度混合泛函 + 色散修正)。
*   **复杂分子/过渡金属配合物**: `MACE-OMOL-0`
    *   **适用范围**: 覆盖 89 种元素。
    *   **场景**: 有机金属催化剂、带电荷的分子、阳离子体系。
    *   **特点**: 支持电荷和自旋嵌入 (Charge/Spin embedding)，分子精度极高。

---

### 3. 表面科学与多相催化 (Surfaces & Catalysis)

如果你的研究对象涉及 **吸附、表面反应、固-气/固-液界面**：

*   **跨域全能王**: `MACE-MH-0` 或 `MACE-MH-1`
    *   **适用场景**: **表面科学的最佳选择**。
    *   **特点**: 训练集混合了无机晶体 (OMAT/MATPES)、分子 (OMOL) 和 催化吸附 (OC20) 数据。
    *   **优势**: 解决了单一模型在“块体”和“分子”之间泛化能力差的问题，非常适合模拟分子吸附在表面的过程。

---

### 总结速查表

| 研究场景 | 推荐模型 | 核心理由 |
| :--- | :--- | :--- |
| **通用无机晶体** | `MACE-MP-0b3` | MP数据库标准，修复了主要Bug，稳定。 |
| **高精度晶体结构** | `MACE-MATPES-r2SCAN-0` | r2SCAN 泛函通常比 PBE 更准。 |
| **金属/合金** | `MACE-MATPES-PBE-0` | 无 +U 修正，物理图像更符合金属。 |
| **表面/吸附/催化** | `MACE-MH-1` | 同时懂“晶体”和“分子”，跨域性能最强。 |
| **有机小分子** | `MACE-OFF23` | 专精有机化学，精度高 (wB97M+D3)。 |
| **金属配合物** | `MACE-OMOL-0` | 覆盖元素广，支持带电体系。 |

### 各 Head 详细对比

| Head 名称 (关键字) | 训练数据集 | 理论水平 (泛函) | 擅长领域 | 备注 |
| :--- | :--- | :--- | :--- | :--- |
| **`matpes_r2scan`** | **MATPES** | **r2SCAN** (Meta-GGA) | **高精度无机材料** | **首选**。比 PBE 更准（特别是晶格参数和结合能），且计算成本相近。适合大多数晶体和表面。 |
| **`mp_pbe_refit_add`** | MP (Materials Project) | PBE / PBE+U | 通用无机材料 | 标准 MP 精度。如果你需要和 Materials Project 数据库的数据做对比，选它。 |
| **`spice_wB97M`** | SPICE | **wB97M-D3/V** | **有机小分子/药物** | 极高精度，包含色散修正（范德华力）。适合纯有机体系、生物分子。 |
| **`oc20_usemppbe`** | OC20 (Open Catalyst) | PBE | **表面催化/吸附** | 专门训练于“催化剂表面+吸附物”体系。适合金属表面的吸附能计算。 |
| **`omol`** | OMOL | wB97M-D3/V | **带电分子/金属配合物** | 相比 SPICE，它更能处理带电体系（如离子）和过渡金属配合物。 |
| **`omat_pbe`** | OMAT | PBE | 有机材料/晶体 | 介于分子和晶体之间，如分子晶体。 |

---

## 注意
- 注意不同体系使用不同模型
- 注意带电体系
- 注意周期性体系
- 注意周期性体系的坐标转换
- 注意设置真空层
- 注意绝对数值只有在同一模型下才有意义
- 注意尽量使用大模型，专用模型，而不是多头模型（'matpes_r2scan', 'mp_pbe_refit_add', 'spice_wB97M', 'oc20_usemppbe', 'omol', 'omat_pbe'）
