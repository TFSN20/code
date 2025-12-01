## 环境
### ase mace chgnet
```
pip install ase
pip install mace-torch
pip install chgnet
```
## 文档
### ase
- https://ase-lib.org/
### mace
- https://github.com/ACEsuit/mace
- 模型：https://github.com/ACEsuit/mace-foundations
### self
#### 各 Head 详细对比

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
