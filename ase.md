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
- https://huggingface.co/mace-foundations/mace-mh-1
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
    $$\hat{H}\Psi = E\Psi$$  
    其中 $\hat{H}$ 是哈密顿算符（包含动能和势能）， $\Psi$ 是波函数，E是能量。
*   **计算困境与解决：** 多电子体系的薛定谔方程极其复杂，无法解析求解。
*   **主流方法：密度泛函理论 (DFT, Density Functional Theory)**
    *   **原理：** 霍恩贝格-科恩 (Hohenberg-Kohn) 定理指出，基态能量由电子密度唯一决定。科恩-沈 (Kohn-Sham) 方程将多体问题简化为单电子在有效势场中运动的问题。
    *   **作用：** 计算材料的能带结构、电子态密度、光电性质、磁性以及原子间的成键本质。
    *   **2025现状：** 虽然DFT仍是主流，但为了追求更高精度（如强关联体系），**GW近似**和**动力学平均场理论 (DMFT)** 也在广泛使用。同时，量子计算（Quantum Computing）开始在小分子模拟上展现出超越经典计算机求解薛定谔方程的潜力。

### 2. 原子尺度（Nanometers, $10^{-9}$ m）：统计力学与经典力学
在这个尺度下，我们不再关注电子的云状分布，而是将原子视为“球体”，关注原子核的位置和动量。

*   **核心方程：** **牛顿运动定律**  
    $$F = ma = -\nabla V(r)$$  
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

## 历史

| 时间 | 人物(中文, 英文) | 成果(中文, 英文全称, 英文缩写) | 底层原理, 底层细节 | 关联性/角色 | 应用(过去的主流应用, 现阶段主要的应用) | 其他补充 |
| --- | --- | --- | --- | --- | --- | --- |
| 约 18 世纪 (欧拉) / 20 世纪 (计算应用) | 欧拉等 (Leonhard Euler et al.) | 有限差分法 (Finite-difference Method, FDM) | 核心原理： 用差分（函数在网格点上的取值之差）来近似函数在某点的导数（微分）。底层细节： 将空间离散化为三维网格 (Grid)。网格上的波函数值是待求的变量。将薛定谔方程中的动能项（二阶导数）替换为有限差分近似，从而将微分方程转化为一个大型稀疏矩阵的本征值问题。例如，一维二阶导数近似为: $$\small \frac{\partial^2 \Psi}{\partial x^2} \approx \frac{\Psi(x+h) - 2\Psi(x) + \Psi(x-h)}{h^2}$$ 其中 $h$ 为网格间距。 | 它是求解微分方程的一种通用数值方法；在量子化学中，它提供了无需基组（Basis-set Free）的求解途径，直接在实空间（Real Space）中工作。 | 过去的主流应用： 主要用于求解简单的量子力学问题（如一维势阱）。现阶段主要应用： 用于求解原子和分子的定态/含时薛定谔方程，特别是应用于纳米结构、量子点、非均匀电场下的体系，以及含时密度泛函理论 (TD-DFT) 的实空间模拟。 | 优点： 无需基组，精度可由网格间距 $h$ 系统控制；计算力和应力较困难；计算成本与网格点数成正比，对大体系计算量可能很大。 与PW对比： 平面波在倒空间工作，FDM 在实空间工作。 |
| 1926 | 薛定谔 (Schrödinger) | 薛定谔方程 (Schrödinger Equation, SE) | 核心原理： 描述量子系统状态的波函数/量子系统中粒子的波函数 ($\Psi$) 随时间（含时）或不随时间（定态）的演化，是量子力学的基本方程之一。定态方程： $$\hat{H}\Psi = E\Psi$$ (其中 $\hat{H}$ 为哈密顿量算符，代表系统的总能量；$E$ 为能量本征值，即系统的特定能量； $\Psi$ 为定态波函数)。含时方程： $$i\hbar \frac{\partial \Psi}{\partial t} = \hat{H}\Psi$$ (描述 $\Psi$ 如何随时间 $t$ 变化)。物理诠释 (玻恩诠释)： $\Psi^2$ 代表在某一时刻某一位置找到粒子的概率密度。 | 它是所有非相对论量子化学计算的基石 | 量子化学计算的基础：所有现代量子化学计算方法（如 Hartree-Fock, DFT, Møller-Plesset Perturbation Theory, Coupled Cluster 等）都是通过求解或近似求解薛定谔方程来预测分子性质。 | 它将波粒二象性融入一个统一的数学框架 |
| 1927-1928 | 托马斯 (Llewellyn Hilleth Thomas)，费米 (Enrico Fermi) | 托马斯-费米模型 (Thomas-Fermi Model, TF) | 基于电子密度 $\rho(\mathbf{r})$ 的原子势场理论。它将非相互作用的自由电子气模型应用于原子中的电子，用局部电子密度近似动能。其总能量泛函只依赖于电子密度 $\rho(\mathbf{r})$。 | 它是密度泛函理论 (DFT) 的最早雏形和里程碑。它首次提出用电子密度取代复杂的波函数来描述多体电子系统。 | 过去主流应用： 早期用于粗略估计重原子的电子密度分布和结合能；原子核物理中用于描述核物质。 现阶段主要的应用： 作为理论物理和量子化学中分析和教学的概念性模型；它是许多更高级密度泛函模型的数学起点或极限情况（如大原子核）。 | 模型的缺陷在于没有包含电子交换能和关联能，且不能很好地描述分子结合；在核附近和远离核的区域，密度衰减不准确。 |
| 1926 (由薛定谔方程确定) |  | 非相对论近似 (Non-relativistic Approximation, NRA) | 核心原理： 假设系统中的所有粒子，尤其是电子，运动速度远小于光速 ($v \ll c$)。底层细节： 薛定谔方程基于经典的能量表达式 (动能 $T = \frac{p^2}{2m}$)，忽略了相对论效应，如电子运动质量随速度的增加、自旋-轨道耦合等。总哈密顿量 $\hat{H}$ 中不包含任何相对论修正项。这意味着，薛定谔方程本质上是经典力学能量守恒定律在量子世界的表达。 | 它是所有基于薛定谔方程的量子化学方法的基本限制和前提。它定义了“常规”量子化学的适用范围，即对不含重原子或高速运动粒子的体系。 | 过去的主流应用： 所有对轻元素体系 (如 C、H、N、O 等组成的有机分子) 的计算。现阶段主要应用： 依然是大多数计算 (如有机化学、生物化学体系) 的默认和有效近似。当体系涉及重元素（如第四周期及以后的元素）时，必须放弃此近似并使用相对论方法。 | 局限性： 无法准确描述原子序数 $Z$ 较大的重元素体系的性质，特别是其内层电子、能级分裂和磁性性质。 重要性： 在 $v \ll c$ 的情况下，该近似极大地简化了哈密顿量，使得计算成为可能。它是狄拉克方程 (相对论量子力学方程) 的低速极限。 |
| 1927 | 玻恩、奥本海默 (Max Born, J. Robert Oppenheimer) | 玻恩-奥本海默近似 (Born-Oppenheimer Approximation, BOA) | 核心原理： 由于原子核的质量远大于电子的质量 ($m_{nucleus} \gg m_{electron}$)，原子核的运动速度比电子慢得多，因此可以将原子核的运动与电子的运动解耦 (Decoupling)。底层细节： 假设电子的运动瞬时适应原子核的某个固定位置。在数学上，这允许将总波函数 $\Psi_{total}(\mathbf{R}, \mathbf{r})$ 分解为电子波函数 $\Psi_{electron}(\mathbf{R}, \mathbf{r})$ 和原子核波函数 $\Psi_{nuclear}(\mathbf{R})$ 的乘积: $$\Psi_{total} \approx \Psi_{electron} \cdot \Psi_{nuclear}$$ 电子波函数在核的势场中求解，其电子能量 $E_{electron}$ 依赖于核坐标 $\mathbf{R}$，形成势能面 (Potential Energy Surface, PES)。 | BO 近似是将多体薛定谔方程简化为可求解的电子薛定谔方程的关键步骤，是整个分子和材料电子结构计算的基石。没有它，几乎所有传统的量子化学计算都无法进行。波恩-奥本海默近似是连接薛定谔方程到实际电子结构计算（如HF、DFT等）的关键桥梁。它通过分离电子和核运动，使得电子结构计算成为可能。 | 过去的主流应用： 使化学家首次能够计算分子的平衡结构和振动频率。现阶段主要应用： 用于构建势能面 (PES)，是计算分子结构优化、反应过渡态、反应路径、振动光谱以及分子动力学模拟的核心前提。 | 局限性： 在处理涉及电子激发态的交叉或靠近点 (如圆锥交叉, Conical Intersections) 时会失效，此时原子核和电子的运动发生强耦合。 重要性： 它将复杂的量子力学问题转化为更容易处理的电子结构问题 (在固定核场中) 和核运动问题 (在势能面上)。 |
| 1928 | 狄拉克 (Paul Dirac) | 狄拉克方程 (Dirac Equation, DE) | 核心原理： 它是描述自旋为 $1/2$ 粒子（如电子）运动的相对论性波动方程。它满足狭义相对论的洛伦兹协变性。底层细节： 方程采用四分量旋量 ($\Psi$) 作为波函数，自然地包含了电子的自旋 ($1/2$)。 $$i\hbar \frac{\partial \Psi}{\partial t} = \left( c\hat{\boldsymbol{\alpha}} \cdot \hat{\mathbf{p}} + \hat{\boldsymbol{\beta}}mc^2 \right) \Psi$$ (其中 $\hat{\boldsymbol{\alpha}}$ 和 $\hat{\boldsymbol{\beta}}$ 是狄拉克矩阵， $\hat{\mathbf{p}}$ 是动量算符)。重要推论： 该方程预测了负能量态的存在，从而预言了反物质（如正电子）的存在。 | 它是相对论量子力学的基石；是描述电子运动的最基本方程；是相对论量子化学方法的理论基础。它将电子的自旋和相对论效应自然地引入量子理论。 | 过去的主流应用： 完美解释了氢原子光谱的精细结构；理论上预言了正电子 (Positron) 的存在 (1932 年被实验证实)。现阶段主要应用： 在相对论量子化学中用于精确计算包含重元素（如锕系和镧系元素）的体系，以获得精确的电子结构、核磁共振参数、光谱和重元素化学性质。 | 意义： 统一了量子力学、狭义相对论和电子自旋。 局限性： 狄拉克方程是单粒子方程，当应用于多电子体系时，仍需引入近似。 与SE的关系： 薛定谔方程是非相对论近似下（低速极限）的狄拉克方程 |
| 1928, 1930 | 道格拉斯·哈特里、弗拉基米尔·福克 (Douglas Hartree, Vladimir Fock) | 哈特里-福克方法 (Hartree-Fock Method, HF) | 核心原理： 采用平均场近似。假设多电子波函数可用一个斯莱特行列式 (Slater Determinant) 近似表示，即每个电子只在一个由所有其他电子运动的平均势场中运动，忽略电子间的瞬时相关性。底层细节 (HF 方程)： 通过解福克方程 (Fock Equation) 获得分子轨道 $\phi_i$ 和轨道能 $\epsilon_i$。 $$\hat{f} \phi_i = \epsilon_i \phi_i$$ (其中 $\hat{f}$ 是福克算符)。计算过程： 采用自洽场 (Self-Consistent Field, SCF) 迭代过程，直到计算得到的电子密度不再变化。 | 它是第一个定量求解多电子体系薛定谔方程的方法；是从头计算法 (Ab Initio) 的基石；是更精确方法（如耦合簇 CC、组态相互作用 CI）的起点 (参考波函数)。 | 过去的主流应用： 计算小型分子的基态电子结构、分子轨道和电离能。现阶段主要应用： 作为许多高级相关能方法的第一步；用于周期性体系的计算（如晶体）；为定性分析分子轨道、电子密度和电势面提供有效且计算成本相对较低的工具。 | 局限性： HF 理论完全忽略了电子的相关能 (Electron Correlation)，导致计算出的结合能和解离能误差较大。它低估了键长，高估了振动频率。 重要性： 奠定了量子化学计算的轨道理论基础，引入了轨道、轨道能等现代化学语言。 |
| 1929 | 鲍林、斯莱特等 (Linus Pauling, John C. Slater et al.) | 原子轨道线性组合 (Linear Combination of Atomic Orbitals, LCAO) | 核心原理： 假设分子轨道 (MO) 可以由分子中所有原子上的原子轨道 (Atomic Orbitals, AO) 以线性组合的方式构造出来。底层细节 (数学形式)： 分子轨道 $\psi_i$ 表示为: $$\psi_i = \sum_{\mu=1}^{N} c_{i\mu} \phi_{\mu}$$ 其中: $\psi_i$ 是第 $i$ 个分子轨道； $\phi_{\mu}$ 是原子轨道（基函数）； $c_{i\mu}$ 是变分系数（或线性组合系数），表示 $\phi_{\mu}$ 对 $\psi_i$ 的贡献程度。这些系数通过求解哈特里-福克 (HF) 方程或 科恩-沙姆 (KS) 方程（通常是它们的矩阵形式，即 Roothaan 方程）获得。 | 分子轨道理论 (MOT) 的基石；它是将薛定谔方程应用于多原子分子时，在空间表示上所做的第一个且最重要的近似。所有现代 ab initio 和 DFT 计算都基于此近似。 | 过去的主流应用： 解释化学键的形成、分子对称性和定性预测轨道能级，如 $\text{H}_2$ 分子的 $\sigma$ 键和 $\sigma^*$ 键。现阶段主要应用： 它是所有量子化学计算（包括 HF、MP2、CC、DFT 等）中分子轨道构建的核心方法。计算结果的精度直接取决于所选原子轨道集合（即基组）的大小和质量。 | 局限性： 严格来说，原子轨道只在孤立原子中是准确的；在分子中，轨道会发生形变（极化）。LCAO 只是一个近似。 重要性： 它将复杂的偏微分方程转化为易于处理的线性代数方程组，是实现量子化学计算的实用化关键。 |
| 1930s | 费利克斯·布洛赫 (Felix Bloch) (布洛赫定理) | 平面波基组 (Plane-wave Basis Set, PW) | 核心原理： 利用晶体中势场的周期性。根据布洛赫定理 (Bloch Theorem)，在周期性势场中的电子波函数可以表示为一个平面波与一个具有晶格周期性的函数的乘积。在实践中，这个周期性函数可以被傅里叶级数展开，其展开项就是一组平面波。 $$\small \Psi_k(\mathbf{r}) = e^{i\mathbf{k}\cdot\mathbf{r}} u_k(\mathbf{r})$$ 底层细节： 基组大小由动能截止值 $E_{\text{cut}}$ 控制。所有动能小于 $E_{\text{cut}}$ 的平面波都被包含在基组中。 | 周期性体系（如晶体、表面）计算的首选基组；它与赝势/PAW 方法结合，构成了现代固态 DFT 计算的主流框架。 | 过去的主流应用： 用于早期的能带结构计算和晶体结构预测。现阶段主要应用： 几乎所有使用 PAW/赝势的固态 DFT 软件包 (如 VASP, Quantum ESPRESSO) 都采用平面波基组，用于计算晶格常数、能带、态密度、缺陷和表面吸附等性质。 | 优点： 系统性（增加 $E_{\text{cut}}$ 可系统地提高精度）；完备性；计算力和应力非常方便；基组重叠误差 (BSSE) 为零。 局限性： 仅适用于周期性体系；需要与赝势/PAW 结合使用，不适合孤立分子。 |
| 1934 | 莫勒、普莱塞特 (Christian Møller, Milton S. Plesset) | 莫勒-普莱塞特二阶微扰理论 (Møller-Plesset Perturbation Theory, MP2) | 核心原理： 采用微扰理论 (Perturbation Theory)，将精确哈密顿量 $\hat{H}$ 分解为零阶项  $\hat{H}_ 0$  (通常是福克算符 $\hat{F}$ 的和) 和微扰项 $\hat{V}$ (即电子相关性引起的误差)。底层细节： MP2 计算的是二阶能量修正 $E^{(2)}$ 。这个修正项包含了单电子激发和双电子激发（即电子间的瞬时相互作用），并将其添加到 HF 的总能量 $E_{\text{HF}}$ 上: $$E_{\text{MP2}} = E_{\text{HF}} + E^{(2)}$$ 它是从头计算法 (Ab Initio) 中计算相关能的最简单方法。 | 它是波函数理论中第一个系统性地引入动态相关性的计算方法；是比 HF 更精确的分子能量和结构计算方法，是耦合簇 (CC) 等更高级方法的铺路石。 | 过去的主流应用： 用于对分子结构和热化学性质的初步相关性修正，特别是对 HF 结果不佳的弱相互作用（如范德华力）体系。现阶段主要应用： 用于计算中小分子体系的相关能修正、几何结构优化和频率分析；常用于作为 DFT 方法精度的校验基准之一；是后-HF 方法中计算成本最低的选择。 | 计算成本： $\mathcal{O}(N^5)$，高于 HF ($\mathcal{O}(N^4)$) 和 DFT ($\mathcal{O}(N^3)$)。 局限性： 并非变分法 (Non-Variational)，即计算的能量不一定高于真实能量；如果体系的参考波函数 (即 HF 解) 质量差，则 MP2 结果可能恶化。 |
| 1930s (晶体)/ 1950s-1980s (发展) | 费米、赫尔曼、海涅等 (E. Fermi, F. Herman, V. Heine et al.) | 赝势 (Pseudopotential, PP) / 有效核电荷近似 (Effective Core Potential, ECP) | 核心原理： 利用化学反应主要涉及价电子的原理，将原子核和内层核心电子作为一个整体，用一个有效势场（即赝势 $V_{\text{pseudo}}$）来代替其对价电子的作用。底层细节： 这样做的好处是消除了价电子波函数在核附近的高频振荡，使得可以用更少的基函数（如平面波或高斯基）来精确描述价电子，同时计算的电子数量大幅减少。 | 它是实现高效固态物理和重元素量子化学计算的关键近似；它将多电子问题的复杂性集中到价电子身上，是平面波计算的核心伙伴。 | 过去的主流应用： 早期应用于晶体能带结构计算，用于处理周期性体系。现阶段主要应用： 几乎所有使用平面波基组的固态 DFT 软件包（如 VASP、Quantum ESPRESSO）都依赖赝势（或 PAW，它是赝势的高级形式）来处理原子核和核心电子。也常用于包含重元素（考虑相对论效应）的分子计算。 | 优点： 大幅减少计算量；可以内置相对论效应（相对论赝势），从而解决薛定谔方程的非相对论局限性。 局限性： 并非全电子方法，会丢失核心电子信息；赝势的选择会对计算结果精度产生影响。 |
| 1960s | 约瑟夫·西尼斯特拉、约瑟夫·皮普尔斯 (J. Čížek, J. Paldus, R. J. Bartlett et al.) | 耦合簇理论 (Coupled Cluster Theory, CC) | 核心原理： 采用指数化激发算符 $\hat{T}$ 作用于 Hartree-Fock (HF) 参考波函数 $\Psi_0\rangle$ 上，从而系统性地包含电子相关性。 $$\small\Psi_{\text{CC}}\rangle = e^{\hat{T}}\Psi_0\rangle$$ 底层细节： 激发算符 $\hat{T} = \hat{T}_1 + \hat{T}_2 + \hat{T}_3 + \dots$ 包含了单激发 ($\hat{T}_1$)、双激发 ($\hat{T}_2$) 等。常用的近似： $\text{CCSD}$ (包含单和双激发)， $\text{CCSD}(\text{T})$ (在 $\text{CCSD}$ 基础上通过微扰理论加入三激发修正)， $\text{CCSD}(\text{T})$ 是目前公认的高精度计算标准。 | --- | --- | --- |
| 1964 | 瓦尔特·霍恩伯格、皮埃尔·科恩 (Walter Kohn, Pierre Hohenberg) | 霍恩伯格-科恩定理 (Hohenberg-Kohn Theorems, HK) | 核心原理 (第一定理)： 系统的基态电子密度 $\rho(\mathbf{r})$ 唯一地决定了外部势能 $V_{ext}(\mathbf{r})$，从而唯一地决定了哈密顿量和所有物理性质。核心原理 (第二定理)： 存在一个普适泛函 $E[\rho]$，对于正确的基态密度，该泛函给出最小能量值 (变分原理)。 $$\frac{\delta E[\rho]}{\delta \rho(\mathbf{r})} = 0$$ 意义： 将基态问题从复杂的 $N$ 电子波函数 $\Psi(\mathbf{r}_1, \ldots, \mathbf{r}_N)$ 降维至仅依赖于三维空间坐标的电子密度 $\rho(\mathbf{r})$。 | 它是密度泛函理论 (DFT) 的数学基石和理论基础；它证明了用电子密度取代波函数来描述量子体系是可行且等价的。 | 过去的主流应用： 作为 DFT 的理论证明，开启了基于电子密度而非波函数的新计算方法研究。现阶段主要应用： 它是所有 DFT 计算 (无论是基态几何优化、频率、还是光谱性质计算) 的理论前提。DFT 是目前计算成本与准确性之间平衡得最好的方法，被广泛应用于化学、物理和材料科学。 | 局限性： HK 定理只证明了普适泛函 $E[\rho]$ 的存在性，但没有给出它的具体解析形式。这导致实际应用中必须对泛函的交换-相关部分 ($E_{xc}[\rho]$) 进行近似。 重要性： 瓦尔特·科恩因此获得了 1998 年诺贝尔化学奖。 |
| 1965 | 瓦尔特·科恩、陆·沙姆 (Walter Kohn, Lu Sham) | 科恩-沙姆方程 (Kohn-Sham Equation, KS) | 核心原理： 引入一个虚拟的、无相互作用的参考体系，该体系的电子密度与真实体系的电子密度相同。通过解这个虚拟体系的单电子方程 (KS 方程)，来获得真实体系的电子密度和总能量。底层细节 (KS 方程)： 形式上与 Hartree-Fock (HF) 方程相似: $\hat{h}_ {KS} \phi_ i = \epsilon_ i \phi_ i$ ,其中  $\hat{h}_ {KS}$  是科恩-沙姆哈密顿量，包含动能、核吸引、库仑排斥和交换-相关势 $V_ {xc}$ 。总能量 $E$ 是动能、核吸引能、库仑能以及交换-相关泛函 $E_ {xc}[\rho]$ 的和。 | 它是将 Hohenberg-Kohn 定理 转化为可操作的计算方法的桥梁；是所有现代密度泛函理论 (DFT) 计算的核心方程。 | 过去的主流应用： 比 HF 方法更准确地处理电子相关性，计算分子性质，特别是过渡金属配合物。现阶段主要应用： 工业界和学术界中最常用的量子化学方法。广泛应用于计算分子结构、反应机理、光谱、材料性质（如带隙、费米面）等，是计算化学领域的主力工具。 | 局限性： KS 轨道 $\phi_i$ 及其能量 $\epsilon_i$ 严格来说不是真正的物理意义上的轨道和轨道能（除非是最高占据轨道 $HOMO$）。计算的准确性完全依赖于使用的交换-相关泛函 $E_{xc}$ 的精度。 |
| 1965 (科恩-沙姆论文中提出) | 瓦尔特·科恩、陆·沙姆等 (Walter Kohn, Lu Sham et al.) | 局域密度近似 (Local Density Approximation, LDA) | 核心原理： 假设体系的交换-相关能 $E_{xc}$ 仅局域地取决于该点的电子密度 $\rho(\mathbf{r})$。底层细节： 体系在每一点 $\mathbf{r}$ 的性质，被近似为具有相同密度的均匀电子气 (Homogeneous Electron Gas, HEG) 的性质。 $$E_{xc}^{LDA}[\rho] = \int \rho(\mathbf{r}) \epsilon_{xc}^{HEG}(\rho(\mathbf{r})) d\mathbf{r}$$ (其中 $\epsilon_{xc}^{HEG}$ 是均匀电子气的单位体积交换-相关能，可以通过精确的蒙特卡洛模拟获得)。 | 它是历史上第一个可用于实际量子化学计算的交换-相关泛函；是所有后续更复杂泛函（如 GGA, meta-GGA, Hybrid 泛函）的理论起点。 | 过去的主流应用： 主要用于固体物理和材料科学，特别是计算金属和半导体的晶格常数、带隙等周期性体系性质。现阶段主要应用： 依然是固体物理中某些体系（如金属）计算的首选；作为更高级泛函（如 GGA）的构成模块之一；为后续泛函的发展提供了基准。 | 局限性： LDA 假设密度变化缓慢，但对原子和分子中快速变化的电子密度区域 (如原子核附近) 精度较差。它系统性地高估结合能 (Overbinding) 和低估能带隙。 重要性： 尽管存在缺陷，它证明了基于密度泛函进行实际计算的可行性和有效性。 |
| 1965 | 伦德奎斯特、赫丁 (Lars Hedin) | GW 近似 (GW Approximation) | 核心原理： 求解戴森方程 (Dyson Equation) 以获得单粒子激发能（即能隙和准粒子能）。GW 近似用自能 $\hat{\Sigma}$ 来代替 Kohn-Sham (KS) 方程中的交换-相关势 $V_{xc}$。 $$\small \hat{\Sigma}(\mathbf{r}, \mathbf{r}', \omega) \approx i G(\mathbf{r}, \mathbf{r}', \omega) W(\mathbf{r}, \mathbf{r}', \omega)$$ 底层细节： $G$ 是单粒子格林函数，描述电子的传播；$W$ 是动态屏蔽库仑相互作用（Screened Coulomb Interaction），描述电子在其他电子形成的极化云中运动时受到的有效库仑作用。它克服了 DFT 缺乏正确的长程相关性描述的问题。 | 它是计算准粒子激发能和能带隙的黄金标准；是对标准 DFT 结果进行高精度修正的关键方法，是多体微扰理论在凝聚态物理中的主力。 | 过去的主流应用： 理论上用于验证和修正 DFT 的能带结构计算。现阶段主要应用： 精确计算半导体和绝缘体的能带隙（Band Gap）、电离能和电子亲和能。在光伏材料、半导体器件和激发态动力学研究中至关重要。 | 计算成本： 通常为 $\mathcal{O}(N^4)$ 到 $\mathcal{O}(N^6)$，远高于标准 DFT。 与DFT关系： 通常采用 $G_0W_0$ 近似，即使用 DFT (通常是 PBE) 的结果作为零阶近似来计算 $G$ 和 $W$，然后进行一次修正。 |
| 1980s - 1990s | 帕尔、王、派德、贝克等 (R. G. Parr, Q. S. Wang, J. P. Perdew, A. D. Becke et al.) | 广义梯度近似 (Generalized Gradient Approximation, GGA) | 核心原理： 假设体系的交换-相关能 $E_{xc}$ 不仅依赖于局域电子密度 $\rho(\mathbf{r})$，还依赖于电子密度的梯度 $\nabla\rho(\mathbf{r})$，从而允许泛函响应密度变化。底层细节： 泛函形式为: $$E_{xc}^{GGA}[\rho] = \int f(\rho(\mathbf{r}),\nabla\rho(\mathbf{r})$$ | --- | --- | --- |
| 1990 | 大卫·范德比尔特 (David Vanderbilt) | 超软赝势 (Ultrasoft Pseudopotential, USPP) | 核心原理： 传统赝势要求赝波函数在核心区（原子核附近）与全电子波函数完全一致，这导致赝波函数在核心区仍然有一定的振荡。USPP 放宽了这一限制，允许赝波函数更平滑。底层细节： 为了弥补在核心区放弃的规范化条件，USPP 引入了增广项（或非局域电荷修正），确保电荷总数和库仑作用的精确性。 | 它是实现低动能截止值 ($E_{\text{cut}}$) 的关键，大幅提高了计算效率；是早年固态计算（如 Quantum ESPRESSO）中常用的高效赝势类型之一。 | 过去的主流应用： 用于处理需要相对较低计算资源的体系，特别是包含 $d$ 轨道或 $p$ 轨道元素（如过渡金属、氧、氮等）的体系。现阶段主要应用： 尽管 PAW 已成为主流，但 USPP 仍然在某些计算中作为 PAW 的高效替代或基准测试的选择之一。 | 优点： 动能截止值 (Cutoff Energy) 比传统范式赝势 (Norm-Conserving PP) 低得多，计算速度快。 局限性： 数学处理上比范式赝势复杂；在某些高度非局域的体系中，可能会引入比 PAW 更大的误差。 |
| 1991 | 佩德、王 (J. P. Perdew, Y. Wang) | PW91 泛函 (Perdew-Wang 1991) | 核心原理： 属于 GGA 家族。旨在通过引入电子密度梯度来修正 LDA，以更好地满足已知的理论约束条件。底层细节： PW91 在设计上大量使用了非经验性的方法，通过插值等数学技巧，构建了交换和相关泛函。它是基于 Perdew-Wang 1991 均匀电子气相关能构建的。 $$\small E_{XC}^{\text{PW91}} = E_X^{\text{PW91}} + E_C^{\text{PW91}}$$ 它在当时满足了大部分已知的理论限制。 | 它是 GGA 泛函的先驱之一；在 PBE 出现之前，PW91 是许多固态和分子计算的首选 GGA 泛函，具有极大的影响力和应用历史。 | 过去的主流应用： 早期固态物理和材料科学计算中，作为 LDA 的高精度非经验替代品，用于计算晶格常数、能带结构和表面性质。现阶段主要应用： 尽管在许多领域被 PBE 取代，但由于其历史一致性，在某些分子动力学模拟和特定固态基准测试中仍有应用。 | 与 PBE 的关系： PBE (1996) 是对 PW91 的简化版本。PBE 的数学形式比 PW91 更简单、更平滑，因此在计算中数值稳定性更好。在精度上，两者在大多数体系中非常接近。 重要性： PW91 的成功直接铺平了 PBE 等后续非经验 GGA 泛函的发展道路。 |
| 1990s 至今 | 斯坎达里斯、佩德等 (N. C. Handy, J. P. Perdew, J. Tao, A. Savin et al.) | 元广义梯度近似 (Meta-Generalized Gradient Approximation, meta-GGA) | 核心原理： 在 GGA 的基础上，泛函除了依赖于电子密度 $\rho(\mathbf{r})$ 和梯度 $\nabla\rho(\mathbf{r})$ 外，还引入了轨道动能密度 $\tau(\mathbf{r})$。 $$\tau(\mathbf{r}) = \sum_i^{occ} \frac{1}{2}\nabla\phi_i(\mathbf{r})$$ | --- | --- | --- |
| 1996 | 佩德、伯克、恩泽霍夫 (J. P. Perdew, K. Burke, M. Ernzerhof) | PBE 泛函 (Perdew-Burke-Ernzerhof) | 核心原理： 属于 GGA 家族，交换和相关部分都依赖于电子密度 $\rho$ 及其梯度 $\nabla\rho$。底层细节： PBE 的关键在于，它在设计时严格满足了 DFT 普适交换-相关泛函的已知理论约束条件（如均匀电子气极限、密度缩放性质等），且不含任何经验参数。 $$\small E_{XC}^{\text{PBE}} = E_X^{\text{PBE}} + E_C^{\text{PBE}}$$ (其中 $E_X^{\text{PBE}}$ 和 $E_C^{\text{PBE}}$ 是经过梯度修正的交换和相关部分)。 | --- | --- | 它是 GGA 泛函的黄金标准；是目前最流行的非经验纯 DFT 泛函；是其混合版本 PBE0 的基础。 |
| 1995 (Anisimov et al.) | 阿尼西莫夫等 (V. I. Anisimov et al.) | 局域密度近似加 U (LDA+U) / PBE+U | 核心原理： 针对 DFT 对强关联电子 (如过渡金属、稀土元素的 $d$ 和 $f$ 轨道电子) 描述不准确的缺陷进行修正。标准 DFT 错误地将 $d/f$ 电子视为离域的，从而低估了轨道局域化和库仑排斥。底层细节： 在标准的 DFT 能量表达式中，额外添加一个经验性的哈伯德项 ($U$) 来惩罚 $d$ 或 $f$ 轨道上的电子过度离域化。 $$E^{\text{DFT}+U} = E^{\text{DFT}} + E^U(\rho_{d/f}, U, J)$$ $U$ 代表有效库仑排斥能（On-site Coulomb Repulsion），它通常是一个根据实验或更精确计算确定的经验参数。 | 修正强关联体系的关键技术；允许 DFT 方法处理那些标准方法无法准确描述的体系，如 Mott 绝缘体、反铁磁体等。 | 过去的主流应用： 主要用于修正 LDA 和 PBE 在计算过渡金属氧化物（如 $\text{NiO}$、$\text{CoO}$）和稀土化合物的能带结构和磁性结构时的严重错误。现阶段主要应用： 固态物理和材料科学中，用于准确计算掺杂半导体、电极材料、催化剂表面以及所有含有未满 $d$ 或 $f$ 壳层元素的体系的能带隙、磁矩和缺陷能。 | 局限性： 修正参数 $U$ 是经验性或半经验性的，选择 $U$ 的值对结果影响巨大，缺乏普适的确定方法。 重要性： 在计算成本相对较低的情况下，显著改善了 DFT 对许多功能性材料（特别是电池材料）电子结构和磁性结构的描述。 |
| 1993 至今 | 贝克、斯托特、史考特等 (A. D. Becke, G. Scuseria, M. J. Frisch et al.) | 混合泛函 (Hybrid Functional) | 核心原理： 将纯 DFT 的交换泛函 $E_X^{DFT}$ 替换为一部分精确的 Hartree-Fock 交换能 $E_X^{HF}$。这种混合是基于近似线性化交换-相关能的电子-电子相互作用。底层细节： 泛函形式通常为: $$E_{XC}^{Hybrid} = A \cdot E_X^{HF} + (1-A) \cdot E_X^{DFT} + E_C^{DFT}$$ (其中 $A$ 是混合系数，通常为 $\approx 20\%$ 到 $25\%$)。代表泛函： B3LYP (最常用), PBE0, HSE (屏幕混合泛函)。 | 它是 DFT 泛函阶梯上的第五层 (最常用层)；极大地提升了 DFT 在计算热化学、反应势垒和激发态方面的精度，使其达到了“化学精度”的要求。 | 过去的主流应用： 计算有机分子和生物分子的精确热力学性质和反应势垒，取代 HF 成为主流计算方法。现阶段主要应用： 工业界和学术界计算反应过渡态、NMR/EPR 等光谱参数、准确的 HOMOs/LUMOs 能量和电子激发能 (通过 TD-DFT) 的默认首选方法。 | 缺点： 计算成本比纯 DFT (GGA/Meta-GGA) 高出约 $2-10$ 倍，因为它涉及计算 $E_X^{HF}$。 重要性： B3LYP (1993) 的出现被认为是 DFT 在化学领域取得巨大成功的关键因素之一。 |
| 1993 | 贝克、李、杨、帕尔 (A. D. Becke, C. Lee, W. Yang, R. G. Parr) | B3LYP 泛函 (Becke, 3-parameter, Lee-Yang-Parr) | 核心原理： 是一种混合 (Hybrid) 泛函，它将精确的 Hartree-Fock 交换能与三种不同的局域和梯度密度泛函组件以经验拟合的系数混合。底层细节 (公式)： $$E_{XC}^{\text{B3LYP}} = E_X^{\text{LDA}} + a_0 (E_X^{\text{HF}} - E_X^{\text{LDA}}) + a_X (E_X^{\text{B88}} - E_X^{\text{LDA}}) + a_C E_C^{\text{LYP}} + (1-a_C) E_C^{\text{VWN}}$$ 其中: $E_X^{\text{HF}}$ 为精确交换；$E_X^{\text{LDA}}$ 和 $E_C^{\text{VWN}}$ 为局域密度近似部分； $E_X^{\text{B88}}$ 为 Becke 交换梯度修正； $E_C^{\text{LYP}}$ 为 Lee-Yang-Parr 相关梯度修正。经验拟合系数: $a_0 = 0.20$, $a_X = 0.72$, $a_C = 0.81$。 | 应用最广泛的混合泛函；在准确性和计算成本之间取得了极好的平衡，是化学家群体中最受欢迎的 DFT 方法。 | 过去的主流应用： 首次实现了对大量有机分子的基态几何结构、相对能量和振动频率的精确计算，解决了 GGA 的过度离域问题。现阶段主要应用： 几乎所有中小型分子体系的基态结构优化、热化学、过渡态搜索以及含时密度泛函理论 (TD-DFT) 计算激发能的首选方法之一。 | 局限性： 并非普适泛函，对范德华力 (vdW) 描述较差；对某些电荷转移激发态和键解离过程的描述不准确。 重要性： 它的高可靠性和广泛可用性，使其成为现代计算化学的黄金标准之一。 |
| 1994 | 皮特·布洛赫尔 (Peter E. Blöchl) | 投影缀加波方法 (Projector Augmented Wave, PAW) | 核心原理： 建立了一个将真实的全电子波函数 $\Psi\rangle$ 映射到平滑的赝波函数 $\tilde{\Psi}\rangle$ 的可逆变换 $\hat{\mathcal{T}}$。 $$\small\Psi\rangle = \hat{\mathcal{T}}$$ | --- | --- | --- |
| 1998 | 扎克瑞亚斯、伯克 (M. R. Zacharias, K. Burke) | 修正 PBE 泛函 (Revised PBE, revPBE) | 核心原理： 属于 GGA 家族。保留 PBE 的相关泛函 $E_C^{\text{PBE}}$，但对 PBE 的交换泛函 $E_X^{\text{PBE}}$ 进行修改，使其在密度的梯度变化较大时，交换能的衰减更快。底层细节： 修正目的是为了更好地满足交换孔的渐进行为 (Asymptotic Behavior)，从而改善分子键解离和弱相互作用的描述。 | 它是 GGA 泛函中用于改善热化学和弱相互作用描述的尝试之一；是 PBE 泛函的一种重要修正和替代。 | 过去的主流应用： 用于计算分子的热化学性质和过渡态势垒，试图解决 PBE 对分子键结合能的高估倾向。现阶段主要应用： 在涉及表面吸附、分子间弱相互作用以及需要更精确分子键解离能的计算中，被用作比标准 PBE 更精确的选择。 | 与PBE比较： revPBE 在计算结合能时通常不如 PBE 准确（通常低估结合能），但在计算原子化能方面表现可能更佳。 后续发展： RPBE (Revised revPBE) 在 revPBE 的基础上进一步修正了交换能，以寻求更普适的精度。 |
| 1999 | 汉肯堡 (B. Hanke) | 修订 PBE 泛函 (Revised PBE, RPBE) | 核心原理： 属于 GGA 家族。它沿用了 PBE 的相关泛函 $E_C^{\text{PBE}}$，但对 revPBE 的交换泛函进行了再次修正。底层细节： RPBE 的交换泛函旨在降低 GGA 对过渡金属与表面之间结合能 (吸附能)的高估倾向。它通过修改 $E_X^{\text{revPBE}}$ 中的参数，使其更好地描述电子密度随空间的变化。 | 它是 PBE 家族中专门优化用于描述表面现象和化学吸附的重要泛函，是许多涉及催化计算的基准泛函。 | 过去的主流应用： 用于计算过渡金属表面对分子（如 $\text{CO}$、$\text{H}_2\text{O}$）的吸附能，修正标准 PBE 在此问题上的系统误差。现阶段主要应用： 在表面化学、异相催化、电化学等领域广泛使用，用于计算吸附质的几何结构、吸附能和表面反应的势垒。 | 与PBE比较： RPBE 对吸附能的低估倾向小于 revPBE，同时优于 PBE 的高估倾向，因此在描述表面化学时精度更高。 特点： 与 PBE 一样，RPBE 也是非经验性的，参数均由理论推导而非拟合而来。 |
| 1999 | 佩德、恩泽、艾克哈特 (J. P. Perdew, M. Ernzerhof, K. Burke) | PBE0 泛函 (Perdew, Burke, Ernzerhof 0) | 核心原理： 是一种混合 (Hybrid) 泛函，它基于PBE (广义梯度近似，GGA) 纯泛函，并以理论推导出的固定比例混合精确的 Hartree-Fock (HF) 交换能。底层细节 (公式)： $$E_{XC}^{\text{PBE0}} = \frac{1}{4} E_X^{\text{HF}} + \frac{3}{4} E_X^{\text{PBE}} + E_C^{\text{PBE}}$$ 混合系数： 精确 HF 交换项的混合系数 $A=1/4$ (即 25%) 是根据微扰理论的零阶近似，由理论推导得出的，而非经验拟合。 | 理论上最纯粹的混合泛函之一；是基于PBE 纯泛函的改进，继承了 PBE 对固体和周期性体系的良好描述，同时加入了 HF 交换来修正 GGA 的局限性。 | 过去的主流应用： 作为 B3LYP 的替代品，用于需要较高精度的分子计算，特别是在处理动力学相关性时表现优异。现阶段主要应用： 广泛应用于计算分子的热化学、光谱性质、过渡金属配合物和固态体系的性质。在计算能隙 (Band Gaps) 方面通常优于 B3LYP。 | 优势： 参数通过理论推导而非经验拟合获得，因此在理论上更具普适性，适用于更广泛的体系。 与B3LYP比较： PBE0 在计算键长、振动频率和反应势垒方面通常与 B3LYP 相当，但在描述电荷转移和能隙方面可能更准确。 |
| 2003 (Heyd, Scuseria, Ernzerhof) | 海德、史考特、恩泽霍夫 (J. Heyd, G. E. Scuseria, M. Ernzerhof) | 屏幕混合泛函 (Screened Hybrid Functional, HSE) | 核心原理： 承认电子相互作用在短程是强烈的，应用精确 HF 交换；而在长程是屏蔽 (Screened) 的，故使用纯 DFT 交换。这通过引入一个分离参数 ($\omega$) 来实现。底层细节 (公式)： 交换能被分为短程 ($SR$) 和长程 ($LR$): $$E_X = \frac{1}{4} E_X^{\text{HF}, SR}(\omega) + \frac{3}{4} E_X^{\text{PBE}, SR}(\omega) + E_X^{\text{PBE}, LR}(\omega)$$ 对于 HSE06，精确 HF 交换项的混合比例为 $25\%$，与 PBE0 相同，但仅在短程有效。分离参数 $\omega$ 通常取 $0.11\ a_0^{-1}$。 | 它是 Hybrid 泛函在固体和周期性体系中的优化版本；显著降低了计算精确 HF 交换项的成本，提高了计算效率。 | 过去的主流应用： 作为第一个成功将混合泛函的精度应用于计算固态能带隙的方法。现阶段主要应用： 在材料科学中，用于精确计算半导体和绝缘体的能带隙；计算固体缺陷和表面吸附；处理需要较高精度，但又不能承受完整 PBE0 计算成本的大分子或周期性体系。 | 优势： 相较于 PBE0，在周期性体系中收敛更快，计算成本更低，并且在计算能带隙方面具有更高的精度。 局限性： 依然是混合泛函，计算成本高于纯 GGA 泛函。对 $\omega$ 参数的选择会影响结果。 |
| 2010 | 斯蒂芬·格里姆 (Stefan Grimme) | DFT-D3 色散修正 (DFT-D3) | 核心原理： 假设色散力可以近似为原子对之间的吸引力之和，并根据距离的高次幂（如 $R^{-6}$ 和 $R^{-8}$）衰减。 $$E_{\text{disp}}(D3) = - \frac{1}{2} \sum_{A \neq B} \sum_{n=6,8} s_n \frac{C_n^{AB}}{R_{AB}^n} f_{\text{damp}}(R_{AB})$$ 底层细节： $C_n^{AB}$ 是原子对 $A, B$ 之间的色散系数，这些系数是非经验性计算得出的，它们取决于原子的局域环境（配位数）。$s_n$ 是经验性缩放因子，用于拟合不同的主体 DFT 泛函（如 PBE-D3、B3LYP-D3）。$f_{\text{damp}}$ 是阻尼函数，用于防止在近程距离处吸引力发散。 | 它是目前最流行的外加色散修正方法之一；它能够以极低的计算成本，显著提高标准 DFT 泛函在处理非共价相互作用时的精度。 | 过去的主流应用： 首次将 DFT 的应用范围扩展到能够精确处理大型非共价体系（如生物分子、分子晶体）。现阶段主要应用： 几乎所有 DFT 计算中，只要涉及到分子间弱相互作用、吸附现象、构象异构或分子晶体，都会默认使用 DFT-D3 或其变体（如 D3(BJ)） | 优点： 计算成本几乎为零；可用于任何DFT 泛函。 缺点： 是一种经验性修正；因为它依赖于原子配位数，对某些过渡金属体系的描述可能不够准确。 |
| 2015 | 佩德、苏丹、孙等 (J. P. Perdew, V. N. Staroverov, J. Sun et al.) | 强约束与适度规范泛函 (Strongly Constrained and Appropriately Normed, SCAN) | 核心原理： 属于 Meta-GGA 家族，依赖于电子密度 ($\rho$)、梯度 ($\nabla\rho$) 和**动能密度** ($\tau$)。底层细节： SCAN 泛函被设计为同时满足已知的 17 个理论约束条件（例如：对均匀电子气、慢变势场、非均匀密度缩放等极限的精确描述），这是在它之前的任何 GGA 或 Meta-GGA 泛函都未能完全做到的。 | --- | --- | 纯 DFT 泛函（不含 HF 交换）的最高理论代表；在理论上是泛函阶梯第四层的终极形式；其精度在许多基准测试中可以媲美甚至超越一些混合泛函。 |
| 2016 | 马丁·赫德-戈登组 (M. Head-Gordon et al.) | $\omega$ B97M-V 泛函 (Range-Separated Hybrid, Meta-GGA with VV10 correlation) | 核心原理： 采用范围分离 (Range-Separated) 策略，在短程使用 HF 交换，在长程使用纯 DFT 交换。同时，它是一个元广义梯度近似 (Meta-GGA) 泛函，包含了动能密度 ($\tau$)。底层细节： 它是通过组合优化 (Combinatorial Optimization) 的方法，在一个巨大的候选泛函空间中，基于大量基准数据集（约 5000 个数据点）训练和选择的最佳 12 参数泛函。它包含了 VV10 非定域相关性 ($\text{VV10 Nonlocal Correlation}$)，旨在精确描述范德华力。 | 它是 混合元 GGA (Hybrid Meta-GGA) 泛函中的最高精度代表之一；在精度测试中，它通常超越了许多主流的高精度泛函（如 M06-2X、$\omega$ B97X-D） | 过去/现阶段应用： 主要用于需要高精度热化学、非共价相互作用（特别是范德华力）、异构化能和反应势垒计算的学术研究和高精度基准测试。也被用于 TD-DFT 激发态几何优化。 | 计算成本： 高于 GGA 和一般的混合泛函，但低于耦合簇 ($\text{CCSD(T)}$)。 变体： 该泛函有其他变体，如 $\omega$ B97M-D3(BJ) 或 $\omega$ B97M-(2)（双重混合泛函）。 |
| 2013-2016 | 马丁·赫德-戈登组、施里姆 (M. Head-Gordon et al., S. Grimme et al.) | $\omega$ B97M+D3 泛函 (Range-Separated Hybrid, Meta-GGA with D3 Dispersion Correction) | 核心原理： 采用 $\omega$ B97M 的高精度 范围分离混合元 GGA (Hybrid Meta-GGA) 框架作为主体，并额外添加一个经验性/半经验性的 DFT-D3 色散修正。 $$\small E_{\text{DFT-D3}} = E_{\omega\text{B}97\text{M}} + E_{\text{disp}}(D3)$$ 底层细节： DFT-D3 修正项是基于原子间 $C_6/C_8$ 色散系数的对原子对势模型，以经验方式修正了主体泛函在长程相关性上的缺陷。这种方法比 $\omega$ B97M-V 中使用的 VV10 非定域相关项在计算上更高效。 | 它是 $\omega$ B97M 家族中最常用的范德华校正版本之一；被设计用于在保持 $\omega$ B97M 对热化学高精度的同时，可靠地计算所有类型的非共价相互作用。 | 过去/现阶段应用： 广泛应用于计算生物分子体系、分子晶体、吸附现象、以及涉及大量非共价相互作用的复杂体系的几何结构和相对能量。它提供了一种在精度上接近 $\text{CCSD(T)}$，但计算成本远低于 $\text{CCSD(T)}$ 的方法。 | 与 $\omega$ B97M-V 比较： $\omega$ B97M-V 将色散修正内置在相关泛函中；而 $\omega$ B97M+D3 是主体泛函加上外加校正。在许多情况下，两者都能提供相似的高精度，但 $\omega$ B97M+D3 的计算通常略微更快。 |
| 2016 | 马丁·赫德-戈登组 (M. Head-Gordon et al.) | $\omega$ B97M-VV10 泛函 (Range-Separated Meta-GGA with VV10 Correlation) | 核心原理： 结合了 $\omega$ B97 家族的范围分离、Meta-GGA 框架 (包含动能密度 $\tau$) 和 VV10 非定域相关修正。底层细节： 交换-相关泛函 $E_{\text{XC}}$ 经过组合优化，精确拟合大量化学数据集。VV10 项（由 Vydrov 和 Van Voorhis 提出）是一种非定域相关 (NLC) 泛函，它内嵌在总能量表达式中，以物理方式描述长程色散力（范德华力）。 | 它是现代 DFT 中追求化学精度的最高级别泛函之一；尤其擅长处理需要精确平衡共价键、热化学和非共价相互作用的复杂体系。 | 过去/现阶段应用： 用于对非共价作用主导的复合体（如 $\pi-\pi$ 堆叠、氢键网络）、异构化能、反应势垒和高精度热力学数据进行基准计算。它被认为是接近 $\text{CCSD(T)}$ 精度的最优 DFT 泛函之一。 | 与 $\omega$ B97M+D3 比较： $\omega$ B97M-VV10 将 vdW 修正项内嵌在相关泛函中，比外加的 D3 修正在理论上更具普适性和优雅性。 计算成本： 由于引入了 Meta-GGA 和非定域相关项，计算成本高于常规的 GGA 或 Hybrid 泛函（如 B3LYP）。 |
| 2020 | Bartók 等 (Á. T. Bartók et al.) | 修订的 SCAN 泛函 (Restored and Revised SCAN, r2SCAN) | 核心原理： 属于 Meta-GGA 家族，依赖于电子密度 ($\rho$)、梯度 ($\nabla\rho$) 和**动能密度** ($\tau$)。底层细节： r2SCAN 的设计是为了修正原始 SCAN 泛函在计算过程中对数值积分的敏感性和收敛性问题，同时保留 SCAN 满足 17 个已知理论约束的优势。 | 它在数学上对 SCAN 的交换和相关部分进行了更平滑的函数形式修订。它是 Meta-GGA 泛函中最先进、最受关注的版本之一；是纯 DFT 泛函（不含 HF 交换）中追求普适高精度的重要尝试。 | --- | --- |

## GGA大全
DFT雅各布天梯（Jacob's Ladder）的**第二阶梯（Rung 2）**被称为**广义梯度近似（GGA, Generalized Gradient Approximation）**泛函。

根据文献和工业生产中（如VASP, CP2K, Quantum ESPRESSO, ORCA等软件）最常用的GGA泛函，以下是Libxc中约定俗成的搭配全称、分类及特性总结：

### 常用GGA泛函列表

| 通用名称 | Libxc 约定搭配 (Exchange + Correlation) | 全称 (主要作者/年份) | 经验性 | 局域性 | 修正性 |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **PBE** | `GGA_X_PBE` + `GGA_C_PBE` | Perdew-Burke-Ernzerhof (1996) | 非经验 | 半局域 | 无 |
| **BLYP** | `GGA_X_B88` + `GGA_C_LYP` | Becke (1988) Exchange + Lee-Yang-Parr (1988) Correlation | 经验性 | 半局域 | 无 |
| **BP86** | `GGA_X_B88` + `GGA_C_P86` | Becke (1988) Exchange + Perdew (1986) Correlation | 经验性 | 半局域 | 无 |
| **PW91** | `GGA_X_PW91` + `GGA_C_PW91` | Perdew-Wang (1991) | 非经验 | 半局域 | 无 |
| **PBEsol** | `GGA_X_PBE_SOL` + `GGA_C_PBE_SOL` | PBE for Solids (Perdew et al. 2008) | 非经验 | 半局域 | 无 |
| **RPBE** | `GGA_X_RPBE` + `GGA_C_PBE` | Revised PBE (Hammer-Hansen-Nørskov 1999) | 非经验 | 半局域 | 无 |
| **revPBE** | `GGA_X_PBE_R` + `GGA_C_PBE` | Revised PBE (Zhang-Yang 1998) | 经验性* | 半局域 | 无 |
| **B97-D** | `GGA_XC_B97_D` | Becke 97 (Grimme parametrization 2006) | 经验性 | 半局域 | 含色散修正 (D2) |
| **AM05** | `GGA_X_AM05` + `GGA_C_AM05` | Armiento-Mattsson (2005) | 非经验 | 半局域 | 无 |

### 说明
1.  **Libxc 命名规则**：Libxc 通常将泛函拆分为交换（Exchange, `_X_`）和关联（Correlation, `_C_`）两部分。对于无法拆分或作为一个整体参数化的泛函（如 B97-D），使用 `_XC_` 前缀。
2.  **经验性/非经验性**：
    *   **非经验 (Non-empirical)**：基于物理约束条件推导，不依赖实验数据拟合（如 PBE, PW91, PBEsol, RPBE）。
    *   **经验性 (Empirical)**：包含拟合实验数据（如原子化能、稀有气体数据）的参数（如 BLYP, BP86 中的 B88 交换部分，revPBE 拟合了 $\kappa$ 参数，B97-D 高度参数化）。
3.  **局域性**：所有 Rung 2 GGA 泛函本质上都是**半局域（Semi-local）**的，因为它们仅依赖于电子密度 $\rho(\mathbf{r})$ 及其梯度 $\nabla \rho(\mathbf{r})$。
    *   *注*：B97-D 虽然基础形式是 GGA，但其包含的经验色散校正（Dispersion Correction）引入了长程相互作用描述，不过在泛函分类中通常仍归类为带修正的 GGA。
4.  **修正性**：大部分标准 GGA（如 PBE, BLYP）本身不包含长程范德华力（色散）描述，常需外挂 DFT-D3, DFT-D4 或 VV10 等修正（例如写作 PBE-D3）。B97-D 是特例，其定义中已内含 D2 色散修正。

DFT雅各布天梯的**第三阶梯（Rung 3）**被称为**Meta-GGA（Meta-Generalized Gradient Approximation，元广义梯度近似）**泛函。这一阶梯引入了动能密度 $\tau$（或拉普拉斯算符 $\nabla^2 \rho$），属于**半局域（Semi-local）**泛函。

以下是文献和工业生产中（如VASP, Gaussian, CP2K等）最常用的Meta-GGA泛函及其Libxc搭配：

### 常用Meta-GGA泛函列表

| 通用名称 | Libxc 约定搭配 (Exchange + Correlation) | 全称 (主要作者/年份) | 经验性 | 局域性 | 修正性 |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **SCAN** | `MGGA_X_SCAN` + `MGGA_C_SCAN` | Strongly Constrained and Appropriately Normed (Sun et al. 2015) | 非经验 | 半局域 | 无 |
| **r2SCAN** | `MGGA_X_R2SCAN` + `MGGA_C_R2SCAN` | Regularized-Restored SCAN (Furness et al. 2020) | 非经验 | 半局域 | 无 (数值稳定性优于SCAN) |
| **TPSS** | `MGGA_X_TPSS` + `MGGA_C_TPSS` | Tao-Perdew-Staroverov-Scuseria (2003) | 非经验 | 半局域 | 无 |
| **revTPSS** | `MGGA_X_REVTPSS` + `MGGA_C_REVTPSS` | Revised TPSS (Perdew et al. 2009) | 非经验 | 半局域 | 无 |
| **M06-L** | `MGGA_X_M06_L` + `MGGA_C_M06_L` | Minnesota 06 Local (Zhao-Truhlar 2006) | 经验性 | 半局域 | 隐式含中程色散 (拟合) |
| **PKZB** | `MGGA_X_PKZB` + `MGGA_C_PKZB` | Perdew-Kurth-Zupan-Becke (1999) | 非经验 | 半局域 | 无 (TPSS的前身) |
| **B97M-V** | `MGGA_XC_B97M_V` | Becke 97 Meta with VV10 (Mardirossian-Head-Gordon 2015) | 经验性 | **非局域*** | 含非局域 VV10 色散修正 |

### 说明
1.  **Libxc 命名规则**：前缀为 `MGGA_`。大多数常用的 Meta-GGA 泛函在 Libxc 中也是分为交换（`_X_`）和关联（`_C_`）两部分调用，但像 B97M-V 这种高度参数化且设计为整体的泛函常以 `_XC_` 形式存在。
2.  **经验性/非经验性**：
    *   **非经验 (Non-empirical)**：SCAN 和 TPSS 系列基于满足所有已知的精确物理约束条件构建，不依赖分子数据库拟合。
    *   **经验性 (Empirical)**：M06-L 是典型的明尼苏达系列泛函，通过拟合大量化学数据库及其参数化，虽然在过渡金属化学中表现优异，但依赖参数拟合。
3.  **局域性**：
    *   标准 Rung 3 泛函（如 SCAN, TPSS, M06-L）依赖 $\tau(\mathbf{r})$，仍属于**半局域**近似。
    *   **特例**：B97M-V 虽然基础形式是 Meta-GGA，但设计中强制包含了非局域关联（VV10），严格来说跨越了 Rung 3 和 Rung 4 之间的界限，但在 Libxc 分类及部分文献讨论中常作为现代高精度 Meta-GGA 的代表提及。
4.  **修正性**：
    *   **r2SCAN**：是 SCAN 的正则化版本，主要修正了 SCAN 在数值积分网格上的不稳定性（Numerical Noise），目前在 VASP 等平面波软件中逐渐取代 SCAN。
    *   **色散修正**：SCAN 和 TPSS 自身不包含长程范德华力，常需搭配 D3 或 D4 修正（如 r2SCAN-D4）。M06-L 通过拟合包含了部分中程相互作用，但长程仍不准确。

DFT雅各布天梯的**第四阶梯（Rung 4）**被称为**杂化泛函（Hybrid Functionals）**。这一阶梯在DFT交换关联项中引入了部分**精确交换（Exact Exchange / Hartree-Fock Exchange）**，属于**非局域（Non-local）**泛函。

杂化泛函在Libxc中通常作为**整体（Monolithic）**定义，因为HF交换的混合比例是泛函定义的固有部分，一般不需要用户分别指定交换和关联部分（但在某些代码中也可以手动混合）。

以下是文献和工业生产中（如Gaussian, VASP, ORCA, Q-Chem等）最常用的杂化泛函及其Libxc约定名称：

### 常用杂化泛函列表 (Hybrid Functionals)

| 通用名称 | Libxc 约定名称 (整体) | 全称 (主要作者/年份) | 经验性 | 局域性分类 | 修正性/特性 |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **B3LYP** | `HYB_GGA_XC_B3LYP` | Becke 3-parameter, Lee-Yang-Parr (1993/1994) | 经验性 | 全局杂化 | 有机化学标准，含3个拟合参数 |
| **PBE0** | `HYB_GGA_XC_PBEH` | PBE Hybrid (Adamo-Barone 1999) | 非经验 | 全局杂化 | 25% HF 混合，物理约束构建 |
| **HSE06** | `HYB_GGA_XC_HSE06` | Heyd-Scuseria-Ernzerhof (2006) | 非经验* | **屏蔽杂化** | 固态带隙标准，仅短程含HF交换 (Screened) |
| **M06-2X** | `HYB_MGGA_XC_M06_2X` | Minnesota 06 2x HF (Zhao-Truhlar 2008) | 经验性 | 全局**Meta**杂化 | 54% HF，优于主族热化学/动力学/弱相互作用 |
| **TPSSh** | `HYB_MGGA_XC_TPSSH` | TPSS Hybrid (Staroverov et al. 2003) | 非经验 | 全局**Meta**杂化 | 10% HF，TPSS的杂化版本 |
| **ωB97X-D**| `HYB_GGA_XC_WB97X_D` | ωB97X with Dispersion (Chai-Head-Gordon 2008) | 经验性 | **区间分离**杂化 | 含长程修正(LC) + 经验色散(D2) |
| **CAM-B3LYP**| `HYB_GGA_XC_CAM_B3LYP`| Coulomb-Attenuating Method B3LYP (Yanai 2004) | 经验性 | **区间分离**杂化 | 修正了B3LYP的长程电荷转移激发问题 |
| **SCAN0** | `HYB_MGGA_XC_SCAN0` | SCAN Hybrid (Hui-Chai 2016) | 非经验 | 全局**Meta**杂化 | 25% HF，SCAN的杂化版本 (注: r2SCAN0 也渐流行) |

### 说明
1.  **Libxc 命名规则**：
    *   杂化泛函的 ID 前缀通常为 `HYB_GGA_XC_` 或 `HYB_MGGA_XC_`。
    *   **重要提示**：在 Libxc 中，杂化泛函是一个整体 ID（例如 ID 402 对应 B3LYP）。调用时，代码会自动处理 HF 交换积分的混合比例（如 B3LYP 的 20%，PBE0 的 25%）。这与 Rung 2/3 常常拆分 Exchange/Correlation 的做法不同。
    *   **PBE0 的名称**：在 Libxc 中约定俗成的名称是 `PBEH` (PBE Hybrid)，对应通用名称 PBE0。
    *   **B3LYP 的定义**：Libxc 的 `HYB_GGA_XC_B3LYP` 对应最通用的 Gaussian 软件定义（使用 VWN 1-RPA 版本的局域关联）。

2.  **经验性/非经验性**：
    *   **经验性 (Empirical)**：如 B3LYP, M06-2X, ωB97X-D。参数是通过拟合特定的热化学或动力学数据库得到的。
    *   **非经验 (Non-empirical)**：如 PBE0, TPSSh, SCAN0。混合参数（如 0.25）是基于微扰理论或物理常数推导的，而非数据拟合。HSE06 属于物理模型驱动，虽然屏蔽参数 $\omega$ 有选定过程，但通常归类为非经验或物理导向。

3.  **局域性与分类**：
    *   **非局域 (Non-local)**：所有 Rung 4 泛函因包含 Hartree-Fock 交换积分，在数学上都是非局域的。
    *   **全局杂化 (Global Hybrid)**：HF 交换比例在全空间是常数（如 B3LYP, PBE0）。
    *   **区间分离杂化 (Range-Separated Hybrid, RSH)**：将电子相互作用分为短程和长程，HF 比例随距离变化。
        *   **长程修正 (Long-range Corrected, LC)**：如 CAM-B3LYP, ωB97X-D，主要用于解决电荷转移（CT）和里德堡态问题。
        *   **屏蔽杂化 (Screened Hybrid)**：如 HSE06，仅在短程保留 HF 交换，长程使用纯 DFT，计算固体时极大降低计算量。

4.  **修正性**：
    *   **色散修正**：ωB97X-D 自身集成了经验色散（-D）。其他如 B3LYP, PBE0 通常需要外挂 D3(BJ) 或 D4 修正（例如写作 B3LYP-D3(BJ)）。
    *   **M06-2X**：通过高度参数化，内部隐式包含了一些中程色散作用，但对极长程作用仍有限。

DFT雅各布天梯的**第五阶梯（Rung 5）**被称为**双杂化泛函（Double Hybrid Functionals）**或包含**虚轨道（Virtual Orbitals）**信息的泛函（如RPA）。这一阶梯不仅引入了精确交换（HF Exchange），还引入了基于微扰理论（如MP2）或随机相位近似（RPA）的**非局域关联**。

这是目前DFT应用的最高阶梯，属于**完全非局域（Fully Non-local）**泛函。

以下是文献和工业生产中（如ORCA, Gaussian, Q-Chem等）最常用的双杂化泛函及其Libxc搭配：

### 常用双杂化泛函列表 (Double Hybrids)

| 通用名称 | Libxc 约定名称 (整体) | 全称 (主要作者/年份) | 经验性 | 局域性 | 修正性/特性 |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **B2PLYP** | `HYB_GGA_XC_B2PLYP` | Becke 2-parameter, Lee-Yang-Parr with MP2 (Grimme 2006) | 经验性 | 完全非局域 | 首个广泛使用的双杂化，含MP2关联 |
| **mPW2-PLYP** | `HYB_GGA_XC_MPW2PLYP` | Modified PW91 2-parameter (Schwabe-Grimme 2008) | 经验性 | 完全非局域 | 改进了B2PLYP的弱相互作用描述 |
| **B2GP-PLYP** | `HYB_GGA_XC_B2GP_PLYP` | B2 General Purpose (Karton-Martin-Radom 2008) | 经验性 | 完全非局域 | 针对热化学数据重新参数化 (General Purpose) |
| **DSD-PBEP86**| `HYB_GGA_XC_DSD_PBEP86`| Dispersion Spin-component-scaled Double-hybrid (Kozuch-Martin 2011) | 经验性 | 完全非局域 | **DSD系列**：含自旋分量缩放(SCS)及色散修正 |
| **ωB97X-2** | `HYB_GGA_XC_WB97X_2` | ωB97X Double Hybrid (Chai-Head-Gordon 2009) | 经验性 | **区间分离**非局域 | 结合了长程修正(LC)与双杂化，减少自相互作用误差 |
| **XYG3** | `HYB_GGA_XC_XYG3` | X3LYP-based Doubly Hybrid (Zhang-Xu-Goddard 2009) | 经验性 | 完全非局域 | 基于绝热连接形式构建，拟合G3热化学数据 |
| **PBE-QIDH** | `HYB_GGA_XC_PBE_QIDH` | PBE Quadratic Integrand Double Hybrid (Bremond et al. 2014) | 非经验* | 完全非局域 | 尝试减少经验参数，基于物理模型构建的双杂化 |

### 说明
1.  **Libxc 命名规则与实现机制**：
    *   **特殊性**：Rung 5 泛函依赖虚轨道（未占据轨道），计算量通常为 $O(N^5)$。Libxc 本身**不计算**微扰部分（MP2/PT2）的能量，它只提供泛函定义的**混合系数**（Scaling factors for HF, DFT Exchange, DFT Correlation, MP2 Correlation）。
    *   宿主软件（如 ORCA, CP2K）会调用 Libxc 获取这些系数，然后自行计算 MP2 积分并进行混合。
    *   前缀通常仍为 `HYB_GGA_XC_`，因为其基础仍是杂化泛函架构。

2.  **经验性/非经验性**：
    *   **经验性 (Empirical)**：绝大多数双杂化泛函（B2PLYP, DSD系列, XYG3）都是高度经验性的。它们通常包含两个或更多个经验参数（$c_{HF}$, $c_{MP2}$ 等），通过拟合大型热化学数据库（如 G3/99, GMTKN）得到。
    *   **非经验 (Non-empirical)**：严格来说，Rung 5 中最著名的非经验方法是 **RPA (Random Phase Approximation)**。但在 Libxc 中，RPA 通常不作为一个单一的泛函 ID 存在，而是通过组合 `PBE` 轨道 + `EXX` + `RPA Correlation` 的方法实现。PBE-QIDH 是少数尝试走非经验路线的双杂化泛函之一。

3.  **局域性**：
    *   **完全非局域 (Fully Non-local)**：不仅交换项是非局域的（如 Rung 4），关联项也是非局域的（依赖占据和虚轨道的所有双电子积分）。

4.  **修正性与高级特性**：
    *   **DSD (Dispersion-corrected, Spin-component-scaled)**：这是现代双杂化泛函的黄金标准（如 DSD-PBEP86）。它引入了两个主要修正：
        1.  **D (Dispersion)**：明确包含色散校正（如 D3 或 D4）。
        2.  **S (Spin-component-scaled)**：对 MP2 关联能中的“自旋平行”和“自旋反平行”成分赋予不同的权重，以提升精度。
    *   **应用场景**：Rung 5 泛函通常用于追求化学精度（Chemical Accuracy, ~1 kcal/mol）的计算，精度通常优于所有低阶梯泛函，但计算成本极高。

除了Rung 2-5的常见泛函外，Libxc中还包含**Rung 1（LDA）**基础泛函，以及**特定领域的专用泛函**（如固体物理中的vdW-DF系列交换项、现代高性能组合泛函）。

以下是补充的Libxc中常被使用的搭配泛函：

### 1. Rung 1：局域密度近似 (LDA)
这是所有DFT泛函的基础，常用于固体物理能带计算或作为高阶泛函的组成部分（如PBE的关联部分基于PW92）。

| 通用名称 | Libxc 约定搭配 (Exchange + Correlation) | 全称 (主要作者/年份) | 经验性 | 局域性 | 修正性 |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **SVWN5** | `LDA_X` + `LDA_C_VWN` | Slater Exchange + Vosko-Wilk-Nusair (1980) fit 5 | 非经验 | **完全局域** | 无 (均匀电子气模型) |
| **PW92** | `LDA_X` + `LDA_C_PW` | Slater Exchange + Perdew-Wang (1992) | 非经验 | **完全局域** | 无 (比VWN更准确的HEG拟合) |
| **SVWN3** | `LDA_X` + `LDA_C_VWN_3` | Slater + VWN fit 3 | 非经验 | **完全局域** | 无 (Gaussian软件中B3LYP定义的各分量) |

### 2. 固体物理专用：vdW-DF 系列交换泛函
在VASP、Quantum ESPRESSO等软件中，计算层状材料或吸附时，常使用特定的**GGA交换泛函**搭配非局域关联核（Non-local Kernel）。Libxc提供了这些特定的交换部分。

| 通用名称 | Libxc 约定搭配 (Exchange Only) | 全称 (主要作者/年份) | 经验性 | 局域性 | 搭配说明 |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **optB88-vdW**| `GGA_X_OPTB88_VDW` | Klimeš-Bowler-Michaelides (2010) | 经验性 | 半局域 | 需搭配 `LDA_C_VWN` + vdW-DF1 Kernel |
| **optB86b-vdW**| `GGA_X_OPTB86B_VDW`| Klimeš-Bowler-Michaelides (2011) | 经验性 | 半局域 | 需搭配 `LDA_C_VWN` + vdW-DF1 Kernel |
| **cx-vdW** | `GGA_X_LV_RPW86` | Berland-Hyldgaard (2014) | 非经验 | 半局域 | 需搭配 vdW-DF1 Kernel (一致性交换) |
| **C09-vdW** | `GGA_X_C09_X` | Cooper (2009) | 经验性 | 半局域 | 需搭配 vdW-DF1 Kernel |

### 3. 现代高性能杂化泛函 (Advanced Hybrids)
近年来Head-Gordon组和Truhlar组开发的顶级精度泛函，常用于高精度计算基准。

| 通用名称 | Libxc 约定名称 (整体) | 全称 (主要作者/年份) | 经验性 | 局域性 | 修正性/特性 |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **ωB97X-V** | `HYB_GGA_XC_WB97X_V` | ωB97X with VV10 (Mardirossian-Head-Gordon 2014) | 经验性 | **区间分离**杂化 | 内含非局域 VV10 色散 (比 -D3 更严谨) |
| **ωB97M-V** | `HYB_MGGA_XC_WB97M_V`| ωB97 Meta with VV10 (Mardirossian-Head-Gordon 2016) | 经验性 | **区间分离**Meta杂化 | 现有泛函中综合精度最高的之一 |
| **MN15** | `HYB_MGGA_XC_MN15` | Minnesota 15 (Yu-He-Truhlar 2016) | 经验性 | 全局Meta杂化 | 广谱精度，改善了多参考态问题的描述 |
| **LRC-ωPBEh**| `HYB_GGA_XC_LRC_WPBEH`| Long-Range Corrected ωPBE (Rohrdanz et al. 2009) | 非经验*| **区间分离**杂化 | 常用作调优 $\omega$ 参数的模板 (Tuned-LC) |

### 总结说明
*   **LDA (Rung 1)**：是所有 DFT 计算的基石，虽然精度不如 GGA/Hybrid，但在 Libxc 中是 `LDA_C_PW` 等关联泛函被广泛调用作为高阶泛函的一部分。
*   **vdW-DF 交换项**：在 Libxc 中，这些泛函通常只提供**交换能**部分（如 `GGA_X_OPTB88_VDW`），因为**非局域关联**（Kernel）通常涉及双重空间积分，由宿主代码（如 VASP, CP2K）直接处理，而不是由 Libxc 的单点核函数处理。
*   **ωB97X-V / ωB97M-V**：代表了当前 DFT 发展的最高水平（组合了 RSH、Meta 和 VV10 非局域关联），在 Libxc 中作为整体封装，是追求极高精度的首选。

## Pretrained Foundation Models

We provide a series of pretrained foundation models for various applications.

### Latest Recommended Foundation Models

| Model Name | Elements Covered | Training Dataset | Level of Theory | Target System | Model Size | GitHub Release | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- |
| MACE-MP-0a           | 89               | MPTrj            | DFT (PBE+U)         | Materials         | small: mace_mp_0/2023-12-10-mace-128-L0_energy_epoch-249.model; medium: mace_mp_0/2023-12-03-mace-128-L1_epoch-199.model; large: mace_mp_0/2024-01-07-mace-128-L2_epoch-199.model | >=v0.3.6       | Initial release of foundation model.                               |
| MACE-MP-0b        | 89               | MPTrj             | DFT (PBE+U)           | Materials            | models: mace_mp_0b/mace_agnesi_medium.model | >=v0.3.10      | Improve pair repulsion and correct isolated atoms. |
| MACE-MP-0b2        | 89               | MPTrj             | DFT (PBE+U)           | Materials            | models: mace_mp_0b2/mace-medium-density-agnesi-stress.model | >=v0.3.9      | Improve stability at high pressure. |
| MACE-MP-0b3        | 89               | MPTrj             | DFT (PBE+U)           | Materials            | models: mace_mp_0b3/mace-mp-0b3-medium.model | >=v0.3.9      | Fixed some phonons issues compared to b2. Improved high pressure stability and reference energies. |
| MACE-MPA-0           | 89               | MPTrj + sAlex    | DFT (PBE+U)         | Materials         | medium-mpa-0: mace_mpa_0/mace-mpa-0-medium.model | >=v0.3.10      | Improved accuracy for materials, improved high pressure stability. |
| MACE-OMAT-0          | 89               | OMAT             | DFT (PBE+U) VASP 54 | Materials         | medium-omat-0: mace_omat_0/mace-omat-0-medium.model | >=v0.3.10      |                                                                    |
| MACE-OFF23           | 10               | SPICE v1         | DFT (wB97M+D3)      | Organic Chemistry | small: mace_off23/MACE-OFF23_small.model; medium: mace_off23/MACE-OFF23_medium.model; large: mace_off23/MACE-OFF23_large.model; medium: mace_off23/MACE-OFF24_medium.model | >=v0.3.6       | Initial release covering neutral organic chemistry.                |
| MACE-MATPES-PBE-0    | 89               | MATPES-PBE       | DFT (PBE)           | Materials         | medium: mace_matpes_0/MACE-matpes-pbe-omat-ft.model | >=v0.3.10      | No +U correction.                                                  |
| MACE-MATPES-r2SCAN-0 | 89               | MATPES-r2SCAN    | DFT (r2SCAN)        | Materials         | medium: mace_matpes_0/MACE-matpes-r2scan-omat-ft.model | >=v0.3.10      | Better functional for materials.                                   |
| MACE-OMOL-0 | 89               | OMOL    | DFT (wB97M-VV10)        | Molecules/Transition metals/Cations         | large: mace_omol_0/MACE-omol-0-extra-large-1024.model | >=v0.3.14      | Charge/Spin embedding, very good molecular accuracy.                                   |
| MACE-MH-0/1 | 89               | OMAT/OMOL/OC20/MATPES    | DFT (PBE/R2SCAN/wB97M-VV10)        | Inorganic crystals, molecules and surfaces. | mh-0: mace_mh_1/mace-mh-0.model; mh-1: mace_mh_1/mace-mh-1.model | >=v0.3.14      | Very good cross domain performance on surfaces/bulk/molecules.   |

The first generation of models are available in the MACE-MP-0.

We subsequently released a second generation of models in the MACE-MP-0b, MACE-MP-0b2 and MACE-MP-0b3 releases.
This release includes models with improved stability during MD simulations, using core repulsion, a new repulsion regularization for high pressure, and a few extra high pressure training examples. Moreover,
we recommend using the second generation models for fine-tuning.

We have also released a model trained on an enlarged dataset containing 3.5M new crystals, obtained by combining the MPtraj and sAlex datasets. This model, released as MACE-MPA-0, achieves state-of-the-art accuracy on the Matbench benchmarks and significantly improves accuracy compared to the MACE-MP-0 models on material systems.

We do not guarantee that the second generation models are better than the first generation models in all cases, but they are expected to be more stable during MD simulations.

### MACE-MP: Materials Project Force Fields

We have collaborated with the Materials Project (MP) to train a universal MACE potential covering 89 elements on 1.6 M bulk crystals in the [MPTrj dataset](https://figshare.com/articles/dataset/23713842) selected from MP relaxation trajectories.

> **!CAUTION**
>
> The MACE-MP models are trained on MPTrj raw DFT energies from VASP outputs, and are not directly comparable to the MP's DFT energies or CHGNet's energies, which have been applied MP2020Compatibility corrections for some transition metal oxides, fluorides (GGA/GGA+U mixing corrections), and 14 anions species (anion corrections). For more details, please refer to the MP Documentation and MP2020Compatibility.yaml.

### MACE-OFF: Transferable Organic Force Fields

There is a series (small, medium, large) transferable organic force fields. These can be used for the simulation of organic molecules, crystals and molecular liquids, or as a starting point for fine-tuning on a new dataset.

### MACE-MH-1

#### Available Model Heads

MACE-MH-1 contains multiple task-specific heads trained on different levels of theory:

| Head Name | Level of Theory | Best For | Access |
| --- | --- | --- | --- |
| omat_pbe (default) | PBE/PBE+U | General materials, balanced performance across tasks | Specify in model |
| omol | ωB97M-VV10 | 1% of OMOL data: Molecular systems, organic chemistry, Organometallic | Specify in model |
| spice_wB97M | ωB97M-D3(BJ) | Molecular systems and organic chemistry | Specify in model |
| rgd1_b3lyp | B3LYP | Reaction chemistry | Specify in model |
| oc20_usemppbe | PBE | Surface catalysis, adsorbates | Specify in model |
| matpes_r2scan | r2SCAN meta-GGA | High-accuracy materials | Specify in model |

By default, the OMAT head (PBE) is used, which provides the best cross-domain performance.

#### Best Practices
- For fine-tuning: Use OMAT head first. Test other heads if needed.
- For materials: Use OMAT head. Use D3 corrections for systems with dispersions. Test matpes_r2scan head if r2scan better reference.
- For molecules: Consider using OMOL head (ωB97M-VV10) for improved intramolecular interactions. OMAT head good for condensed phase molecular systems, test it too.
- For surfaces: OMAT head provides excellent performance; OC20 head available for specialized applications
### other
for some models(mh-1), need provide a head keyword('matpes_r2scan', 'mp_pbe_refit_add', 'spice_wB97M', 'oc20_usemppbe', 'omol', 'omat_pbe'......) to specify the head you want to use.
