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
  pip install torch-ema
  pip install torchmetrics
  ```
  ```
  pip install ase
  pip install mace-torch
  pip install chgnet
  pip install gpaw
  ```
## 文档
### ase
- https://ase-lib.org/
### mace
- https://github.com/ACEsuit/mace
- 模型：https://github.com/ACEsuit/mace-foundations
### chgnet
- https://github.com/CederGroupHub/chgnet
### gpaw
- https://gpaw.readthedocs.io/documentation/basic.html
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
# 知识
## 材料计算的底层原理
材料计算（Computational Materials Science）的底层原理是一个跨越多个时间和空间尺度的物理与数学体系。它并非单一的理论，而是根据研究对象的尺度大小，由量子力学、统计力学、经典力学和连续介质力学共同构成的金字塔结构。

截至2025年，随着人工智能（AI）的深度介入，材料计算正在经历从“完全基于物理方程求解”向“物理+数据驱动”双引擎模式的转变。

以下是按照尺度从微观到宏观的底层原理深度解析：

### 1. 电子尺度（Ångstroms, $10^{-10}$ m）：量子力学与第一性原理
这是材料计算的最底层，也是最根本的原理，被称为“第一性原理”（First Principles / Ab initio），意味着不需要经验参数，仅需原子序数和晶体结构即可计算。

*   **核心方程：** **薛定谔方程 (Schrödinger Equation)**
    $$ \hat{H}\Psi = E\Psi $$
    其中 $\hat{H}$ 是哈密顿算符（包含动能和势能），$\Psi$ 是波函数，E是能量。
*   **计算困境与解决：** 多电子体系的薛定谔方程极其复杂，无法解析求解。
*   **主流方法：密度泛函理论 (DFT, Density Functional Theory)**
    *   **原理：** 霍恩贝格-科恩 (Hohenberg-Kohn) 定理指出，基态能量由电子密度唯一决定。科恩-沈 (Kohn-Sham) 方程将多体问题简化为单电子在有效势场中运动的问题。
    *   **作用：** 计算材料的能带结构、电子态密度、光电性质、磁性以及原子间的成键本质。
    *   **2025现状：** 虽然DFT仍是主流，但为了追求更高精度（如强关联体系），**GW近似**和**动力学平均场理论 (DMFT)** 也在广泛使用。同时，量子计算（Quantum Computing）开始在小分子模拟上展现出超越经典计算机求解薛定谔方程的潜力。

### 2. 原子尺度（Nanometers, $10^{-9}$ m）：统计力学与经典力学
在这个尺度下，我们不再关注电子的云状分布，而是将原子视为“球体”，关注原子核的位置和动量。

*   **核心方程：** **牛顿运动定律**
    $$ F = ma = -\nabla V(r) $$
    力 $F$ 是势能面 $V(r)$ 的梯度。
*   **核心难点：势函数 (Potential Energy Surface)**
    如何描述原子间的相互作用力是核心。
    *   **传统方法：** 经验势函数（如Lennard-Jones, EAM, Tersoff），基于实验数据拟合，计算快但精度有限，无法处理断键/成键复杂的化学反应。
    *   **现代方法（AI介入）：机器学习势函数 (Machine Learning Potentials, MLP)**。利用神经网络（如DeepMD, NequIP）学习DFT的高精度数据，生成既有DFT精度又有经典MD速度的势函数。这是近五年来的革命性突破。
*   **主流方法：分子动力学 (Molecular Dynamics, MD)**
    通过积分牛顿运动方程，模拟原子随时间的演化轨迹。用于研究扩散、相变、力学形变（位错运动）和热传导。
*   **辅助方法：蒙特卡洛模拟 (Monte Carlo, MC)**
    基于随机采样和概率统计（Metropolis准则），用于寻找能量最低构型或研究热力学平衡态（如合金的有序/无序转变）。

### 3. 介观尺度（Micrometers, $10^{-6}$ m）：热力学与相变动力学
这个尺度关注的是材料的**微观组织（Microstructure）**，如晶粒、晶界、析出相、枝晶生长等。

*   **核心原理：** **热力学（最小自由能原理）与 动力学（扩散）**
    系统总是趋向于吉布斯自由能（Gibbs Free Energy）最低的状态。
*   **主流方法：相场法 (Phase Field Method, PFM)**
    *   引入连续变化的“序参数”（Order Parameter）来区分不同的相（如液相=0，固相=1，界面=0~1之间变化）。
    *   **核心方程：**
        *   **Cahn-Hilliard方程**（守恒场，描述组分扩散/相分离）。
        *   **Allen-Cahn方程**（非守恒场，描述晶粒生长/畴壁运动）。
    *   **作用：** 模拟凝固过程中的枝晶生长、固态相变、裂纹扩展、铁电畴翻转等。

### 4. 宏观尺度（Meters, $10^0$ m）：连续介质力学
在宏观层面，材料被视为连续介质，忽略原子结构，关注应力、应变、温度场等宏观物理量。

*   **核心方程：** **本构方程 (Constitutive Equations)** + 守恒定律（质量、动量、能量）。
    例如：胡克定律 $\sigma = E\epsilon$ 及其塑性修正。
*   **主流方法：有限元分析 (Finite Element Method, FEM)**
    将复杂几何体划分为网格，求解偏微分方程组。
*   **跨尺度连接：** 宏观计算需要的参数（如弹性模量、屈服强度、导热系数）往往由原子尺度（DFT/MD）计算提供，或者通过晶体塑性有限元（CPFEM）将晶粒取向信息引入宏观计算。

### 5. 第四范式：数据驱动的材料科学 (AI & Informatics)
这是当前（2025年）最前沿的“底层”逻辑变化。

*   **原理：** 放弃直接求解物理方程，转而通过高维数据寻找 **"结构-效能" (Structure-Property)** 之间的映射关系。
*   **方法：**
    *   **材料基因组 (Materials Genome)：** 高通量计算生成海量数据库。
    *   **图神经网络 (GNNs)：** 直接输入晶体图结构，预测带隙、形成能等性质。
    *   **生成式AI (Generative AI)：** 类似于ChatGPT，但用于生成新的晶体结构或分子式（Inverse Design，逆向设计），即根据所需性能反推材料结构。

### 总结：多尺度耦合 (Multiscale Modeling)
真正的材料计算往往不是孤立的，现在的趋势是**多尺度耦合**：
1.  **DFT** 算出电子结构，训练 **MLP**。
2.  **MLP** 驱动 **MD** 模拟，算出界面能和迁移率。
3.  界面能参数输入 **相场法 (PFM)**，模拟显微组织演化。
4.  显微组织特征输入 **有限元 (FEM)**，预测宏观工件的服役寿命。

这就是材料计算的完整底层逻辑链条。

## DFT材料计算的历史发展

---

### 一、 核心概念与底层原理 (Underlying Principles)

#### 1. 核心思想：从“波函数”到“电子密度”
在DFT出现之前，量子力学计算的核心是求解**薛定谔方程（Schrödinger Equation）**。
*   **困难点**：对于一个包含 $N$ 个电子的系统，波函数 $\Psi$ 是 $3N$ 个空间坐标的函数。随着电子数增加，计算量呈指数级爆炸（“指数墙”），根本无法计算宏观材料。
*   **DFT的解决之道**：DFT证明了系统的基态性质完全由**电子密度 $\rho(\vec{r})$** 决定。电子密度只是 3 个空间坐标的函数，与电子数量 $N$ 无关。
    *   **降维打击**：将 $3N$ 维问题简化为 3 维问题。

#### 2. 关键方程：Kohn-Sham 方程
为了实际求解，Kohn 和 Sham 引入了一个虚构的**“无相互作用电子系统”**，使其电子密度与真实的**“相互作用电子系统”**相同。
*   **能量泛函形式**：
    $$E[\rho] = T_s[\rho] + E_{ext}[\rho] + E_{H}[\rho] + E_{xc}[\rho]$$
    *   $T_s$：无相互作用动能。
    *   $E_{ext}$：外部势能（原子核吸引）。
    *   $E_{H}$：哈特里能（电子-电子库伦排斥）。
    *   **$E_{xc}$（交换-关联泛函，Exchange-Correlation Functional）**：这是DFT的“黑匣子”或“垃圾桶”，包含了所有未知的量子力学复杂相互作用（交换作用和电子关联）。**DFT发展的历史，本质上就是寻找更精确 $E_{xc}$ 的历史。**

---

### 二、 历史发展时间轴与关键人物 (Timeline & Key Figures)

#### 1. 史前时代：雏形与尝试 (1920s)
*   **1927年 - Thomas-Fermi 模型**
    *   **人物**：Llewellyn Thomas, Enrico Fermi.
    *   **细节**：最早尝试用电子密度描述能量。
    *   **缺点**：完全忽略了电子交换和关联，甚至无法描述化学键（分子不能成键），粗糙且不准确。

#### 2. 奠基时代：理论确立 (1960s)
*   **1964年 - Hohenberg-Kohn (HK) 定理**
    *   **人物**：Pierre Hohenberg, **Walter Kohn** (1998年诺贝尔化学奖得主).
    *   **底层原理**：
        1.  **定理一**：基态电子密度唯一决定外部势（从而决定哈密顿量和波函数）。
        2.  **定理二**：变分原理，真实的基态密度使能量泛函最小化。
    *   **意义**：从数学上证明了DFT的可行性，但没给出具体怎么算。

*   **1965年 - Kohn-Sham (KS) 方程**
    *   **人物**：**Walter Kohn**, **Lu Jeu Sham (沈吕九)**.
    *   **底层原理**：提出了轨道概念（KS轨道），将多体问题转化为单电子问题。引入了**局部密度近似 (LDA)**。
    *   **意义**：让DFT真正变得可计算，是现代DFT软件的基石。

#### 3. 发展时代：泛函的进化 (Jacob's Ladder) (1980s - 1990s)
John Perdew 将泛函精度的提升比作通往天堂的“雅各布天梯”（Jacob's Ladder）：

*   **第一阶：LDA (Local Density Approximation)**
    *   **时间**：1965 (提出), 1980 (Ceperley-Alder通过QMC给出精确参数).
    *   **原理**：假设微区内的电子气是均匀的。
    *   **特点**：结构计算尚可，但结合能不准，带隙严重低估。

*   **第二阶：GGA (Generalized Gradient Approximation)**
    *   **时间**：1980s末 - 1996.
    *   **人物**：Axel Becke, John Perdew, Yue Wang.
    *   **关键节点**：
        *   **PW91** (Perdew-Wang 1991).
        *   **PBE** (Perdew-Burke-Ernzerhof, 1996)：物理和材料界最著名的泛函，至今仍是标准配置。
    *   **原理**：不仅考虑密度 $\rho$，还考虑密度的梯度 $\nabla \rho$（即电子分布的不均匀性）。

*   **第三阶：Hybrid Functionals (杂化泛函)**
    *   **时间**：1993 - 2003.
    *   **人物**：Axel Becke.
    *   **关键节点**：
        *   **B3LYP** (1993)：化学界的神级泛函，混合了部分Hartree-Fock交换能。
        *   **HSE06** (Heyd-Scuseria-Ernzerhof, 2003)：解决了固体带隙计算不准的问题，半导体计算标配。

#### 4. 现代：精度修正与AI融合 (2010s - 至今)
*   **范德华力修正 (vdW corrections)**：
    *   **人物**：Stefan Grimme.
    *   **关键**：DFT-D3, DFT-D4。弥补了DFT无法描述弱相互作用（如石墨层间作用、气体吸附）的缺陷。
*   **Meta-GGA (SCAN泛函)**：
    *   **时间**：2015.
    *   **细节**：满足了所有已知的物理约束条件，计算成本接近GGA但精度更高。
*   **ML-DFT (机器学习势函数)**：
    *   **时间**：2018 - 至今.
    *   **代表**：DeepMD, MACE, NequIP.
    *   **原理**：利用DFT的高精度数据训练神经网络，以经典力学的速度实现量子力学的精度。

---

### 三、 关键名词中英文对照 (Terminology)

| 简称 | 英文全称 | 中文名称 | 说明 |
| :--- | :--- | :--- | :--- |
| **DFT** | Density Functional Theory | 密度泛函理论 | 总称 |
| **HK Theorems** | Hohenberg-Kohn Theorems | 霍恩贝格-科恩定理 | 理论基石 |
| **KS Eq.** | Kohn-Sham Equations | 科恩-沈方程 | 实际计算方程 |
| **XC** | Exchange-Correlation | 交换-关联 | 未知相互作用的总和 |
| **LDA** | Local Density Approximation | 局域密度近似 | 最简单的泛函 |
| **GGA** | Generalized Gradient Approximation | 广义梯度近似 | 引入密度梯度，比LDA准 |
| **PBE** | Perdew-Burke-Ernzerhof | PBE泛函 | 固体物理最常用的GGA泛函 |
| **Hybrid** | Hybrid Functionals | 杂化泛函 | 混合了HF交换能，如B3LYP, HSE06 |
| **VASP** | Vienna Ab initio Simulation Package | 维也纳从头算模拟包 | 商业软件，材料界最常用 |
| **QE** | Quantum ESPRESSO | 量子咖啡（意译） | 开源软件，主要用于赝势平面波 |
| **HOMO/LUMO** | Highest Occupied/Lowest Unoccupied Molecular Orbital | 最高占据/最低未占据轨道 | 决定化学反应活性和带隙 |
| **DOS** | Density of States | 态密度 | 描述能量状态分布 |

---

### 四、 优缺点分析 (Pros & Cons)

#### 优点 (Pros)
1.  **性价比极高**：在保证量子力学精度的前提下，计算量随原子数 $N$ 呈 $O(N^3)$ 增长（相比波函数法的 $O(e^N)$ 或 $O(N^4+)$），使得处理几百个原子的超胞成为可能。
2.  **普适性**：不需要经验参数（First-principles/Ab initio），只要知道原子序数和位置，理论上可以计算宇宙中任何材料。
3.  **基态性质准确**：对晶格常数、键长、弹性模量等结构性质预测非常准确（误差通常<1-2%）。

#### 缺点 (Cons)
1.  **带隙低估 (Band Gap Underestimation)**：标准DFT（如PBE）计算的半导体带隙通常比实验值偏小30%-50%。需用HSE06或GW方法修正。
2.  **自相互作用误差 (Self-Interaction Error, SIE)**：电子会错误地与自己产生的电场发生相互作用，导致强关联体系（如过渡金属氧化物）计算失败。需用 DFT+U 方法修正。
3.  **弱相互作用缺失**：早期泛函无法描述范德华力（vdW），导致层状材料或气体吸附距离计算错误。需用 DFT-D 修正。
4.  **激发态困难**：标准DFT是基态理论。计算光吸收、发射等激发态需要 **TD-DFT (Time-Dependent DFT)**，且精度有限。

---

### 五、 关联与应用 (Connections & Applications)

DFT 已经成为连接微观物理与宏观工程的桥梁：

1.  **新能源电池**：
    *   计算锂/钠离子在电极材料中的**扩散能垒**（决定充电速度）。
    *   预测电压平台和体积膨胀率。
2.  **半导体与芯片**：
    *   计算**能带结构 (Band Structure)**，筛选合适的半导体材料。
    *   研究掺杂对导电性的影响。
3.  **催化 (Catalysis)**：
    *   **Nørskov 理论**：计算反应中间体在催化剂表面的**吸附能**，寻找“火山图”顶点的最佳催化剂（如析氢、CO2还原）。
4.  **高通量筛选 (High-Throughput Screening)**：
    *   结合 **Materials Project (MP)** 数据库，利用DFT数据训练AI模型，从成千上万种候选材料中筛选新材料。

### 总结
DFT 从1964年的一个数学定理，发展到今天成为物理、化学、材料学博士生的必备工具，其核心在于**Walter Kohn**和**沈吕九**将复杂的量子多体问题简化为了单电子问题。虽然它有带隙低估等缺陷，但随着**杂化泛函**和**机器学习**的加入，DFT依然是目前人类探索微观物质世界最强大的理论工具。
