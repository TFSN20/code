# 引言（Introduction）

GPAW 的计算**通过 Python 脚本来驱动和控制**。在使用上，你通常会：

1. 用 **ASE（Atomic Simulation Environment）** 创建/读取一个 **`ase.Atoms`** 对象来描述原子结构（元素、坐标、晶胞、周期性等）。
2. 创建一个 **`gpaw.GPAW`** 计算器（calculator），设置计算模式、XC 泛函、k 点、网格等参数。
3. 将计算器绑定到 `Atoms`：`atoms.calc = calc`
4. 通过调用能量/力/应力等性质触发计算，例如 `atoms.get_potential_energy()`、`atoms.get_forces()`。

GPAW 强依赖 ASE。ASE 不仅负责“描述原子”，还提供了分子动力学（MD）、几何优化、轨迹与结果分析、可视化等大量实用功能。**如果你对 ASE 不熟悉，建议至少先阅读 ASE 的 About/基本概念部分**，否则后续在建模（晶胞、真空层、周期性、k 点路径、结构优化等）上会很吃力。

> 文档中有些 Python 代码示例会以 `>>>` 开头（`...` 表示续行），这是交互式解释器（REPL）的风格。建议你打开 Python 解释器或 Jupyter，把示例逐段运行并观察输出。

---

## 单位约定（ASE conventions）

GPAW 计算器使用的单位遵循 **ASE conventions**，最重要的是：

* 能量：**electron volt（eV）**
* 长度：**ångström（Å）**
* 力：通常是 **eV/Å**
* 电子温度/展宽：常用单位 **eV**

---

# 进行一次 PAW 计算（Doing a PAW calculation）

要用 GPAW 做一次 PAW（Projector Augmented-Wave）计算，你至少需要：

* 一个 **ASE `Atoms`** 对象（结构）
* 一个 **GPAW calculator**（计算设置）

下面是一个最小示例：计算单个 **H₂ 分子** 的力（forces）。

```python
from ase import Atoms
from gpaw import GPAW

d = 0.74
a = 6.0

atoms = Atoms('H2',
              positions=[(0, 0, 0),
                         (0, 0, d)],
              cell=(a, a, a))
atoms.center()

calc = GPAW(mode='fd', nbands=2, txt='h2.txt')
atoms.calc = calc
print(atoms.get_forces())
```

如果执行上述代码，会启动一个单分子 H₂ 的 DFT/PAW 计算。该计算使用一个大小为 `6.0 × 6.0 × 6.0 Å` 的**超胞（supercell）**，并采用**簇/分子边界条件（cluster boundary conditions）**的思路（即：在大真空盒子里放一个分子，避免与周期镜像相互作用）。

在**不额外指定**更多参数时，该计算等效于使用了一组默认且“合理”的选择（这些默认值会随 GPAW 版本演进，建议以输出文件为准）。在该例中，典型默认设置包括：

* `nbands=2`：2 条电子能带（bands）
* 交换-相关（exchange-correlation, **XC**）泛函：默认 **LDA**
* 自旋：默认**非自旋极化**（spin-paired / spin-unpolarized）
* 实空间网格（real-space grid）：默认会选择对应的网格点数（示例中常见为 `32 × 32 × 32`，但以实际输出为准）

这些参数都可以在文本输出文件 `h2.txt` 中找到。

---

## 显式指定参数（推荐做法）

GPAW 会尽量为你未指定的参数做“合理默认”，但在实际材料计算中，**明确写出关键参数**能显著提升可重复性与可维护性。比如：

```python
calc = GPAW(mode='fd',
            nbands=1,
            xc='PBE',
            gpts=(24, 24, 24))
```

这里表示：

* `nbands=1`：每个自旋通道只算 1 条 band（对 H₂ 这种体系足够）
* `xc='PBE'`：使用 GGA 类的 **PBE** 泛函
* `gpts=(24,24,24)`：每个方向使用 24 个网格点

> 实用建议：在“材料体系（周期晶体、表面、界面）”里，通常你至少会显式指定 `mode / xc / kpts / (h 或 gpts) / occupations / convergence / txt` 这几类参数。

---

# 参数速查表（Parameters）

下表是常用关键字（keyword）的总览。**keyword 名称保持英文**（因为脚本里就是这些关键字），说明部分为中文并保留必要英文术语。

| keyword         | type                       | default value | 说明（description）                               |
| --------------- | -------------------------- | ------------- | --------------------------------------------- |
| `basis`         | `str` or `dict`            | `{}`          | 原子基组（Atomic basis set），主要用于 LCAO 模式           |
| `charge`        | `float`                    | `0`           | 体系总电荷（以电子电荷为单位，`charge=-1` 表示多 1 个电子）         |
| `communicator`  | Object                     |               | MPI 通信器（communicator），用于只用部分进程跑某个计算           |
| `convergence`   | `dict`                     |               | SCF（自洽场）收敛标准                                  |
| `eigensolver`   | `str`                      | `'dav'`       | 本征求解器（eigensolver），如 Davidson / RMM-DIIS / CG |
| `external`      | Object                     |               | 外加势（external potential），如外电场等                 |
| `gpts`          | seq                        |               | 网格点数（grid points），如 `(nx, ny, nz)`            |
| `h`             | `float`                    | `0.2`         | 网格间距（grid spacing，单位 Å），主要对 FD 实空间网格相关        |
| `hund`          | `bool`                     | `False`       | 是否用 Hund’s rule 猜初始磁矩/占据                      |
| `kpts`          | seq / dict                 | Γ-point       | k 点采样（Brillouin-zone sampling）                |
| `maxiter`       | `int`                      | `333`         | 最大 SCF 迭代步数                                   |
| `mixer`         | Object                     |               | 密度混合（Pulay mixing）方案                          |
| `mode`          | `str` or `dict`            |               | 计算模式：`pw` / `lcao` / `fd`（平面波/原子轨道/有限差分）      |
| `nbands`        | `int/str`                  |               | 电子能带数（每自旋通道）                                  |
| `occupations`   | occ. obj. / dict           |               | 占据数与展宽（smearing）设置                            |
| `parallel`      | `dict`                     |               | 并行策略选项                                        |
| `poissonsolver` | Object / dict              |               | Poisson 求解器设置、偶极修正（dipole correction）等        |
| `random`        | `bool`                     | `False`       | 用随机数初始化波函数（Wave function initialization）      |
| `setups`        | `str` or `dict`            | `'paw'`       | PAW 数据集（setup）或赝势（pseudopotential）选择          |
| `spinpol`       | `bool`                     |               | 是否自旋极化（spin-polarized）                        |
| `symmetry`      | `dict/str`                 | `{}`          | 对称性使用策略（point-group / time-reversal）          |
| `txt`           | `str` / `None` / file obj. | `'-'`（stdout） | 文本输出位置：标准输出/文件/关闭输出                           |
| `xc`            | `str`                      | `'LDA'`       | 交换-相关泛函（XC functional）                        |
| `extensions`    | `list[ExtensionInput]`     | `[]`          | 扩展功能（Extensions）                              |

---

# 选择计算模式：PW、LCAO 还是 FD？

GPAW 支持三类“波函数表示（representation）”模式：

* **PW（plane-wave）平面波模式**
* **LCAO（linear combination of atomic orbitals）原子轨道线性组合**
* **FD（finite-difference）有限差分实空间网格**

不同模式在内存、速度、可实现功能、收敛行为上差异很大，选对模式会显著影响效率和结果可靠性。

---

## Plane-waves（PW）平面波模式

用平面波展开波函数。使用 `mode='pw'` 会采用默认平面波截断能 `E_cut = 340 eV`，并使用满足：

`|G + k|^2 / 2 < E_cut`

的平面波基。

若要设置不同截断能（cutoff），可用：

```python
from gpaw import GPAW, PW
calc = GPAW(mode=PW(200))
```

> 实用建议：PW 模式对周期体系（晶体、表面 slab）非常常用；当你需要应力张量（stress tensor）或某些响应函数功能时，PW 往往是唯一选择（见后表）。

---

## LCAO 模式

用“原子型轨道”构成的基组展开波函数，即 **LCAO**（linear combination of atomic orbitals）。设置：

```python
mode='lcao'
```

> 实用建议：LCAO 往往在大体系（大超胞、缺陷、界面）中更省内存更快，适合做初筛、预优化、或需要大量构型扫描的任务。缺点是**基组完备性（complete basis set limit）**更难系统逼近。

---

## Finite-difference（FD）有限差分模式

用实空间网格（real-space grid）表示波函数，通过有限差分离散化。设置：

```python
mode='fd'
```

> **Warning（重要提醒）**
> 未来版本中，“不显式指定 `mode` 参数”可能会变成错误（error）。目前如果你不写 `mode`，GPAW 可能会**隐式选择 FD 模式并给出 warning**。因此建议你**总是显式写出 `mode=...`**，让脚本更稳定可复现。

---

## PW / LCAO / FD 的对比（实用角度）

### 1) 内存消耗（Memory）

* **LCAO**：自由度最少，**内存最低**
* **PW**：较高
* **FD**：通常**最高**（实空间网格自由度大）

### 2) 速度（Speed）

* **小体系 + 很多 k 点**：PW 往往最快
* **大体系**：LCAO 常最有效率
* PW 与 FD：小体系 PW 往往快；大体系 FD 可能更容易并行扩展（parallelize better）而反超

### 3) 绝对收敛（Absolute convergence）

* **FD**：通过减小 `h`（更密网格）可系统提高精度
* **PW**：通过增大 `E_cut` 可系统提高精度
* **LCAO**：基组完备极限更难达到，能量的“绝对收敛”更依赖基组质量与体系

### 4) Eggbox 误差

* **LCAO 与 FD**：会出现小的 “eggbox error”（原子平移时能量出现与网格间距相关的周期性微小起伏）
* **PW**：GPAW 的 PW 实现没有这个问题

### 5) 功能支持（Features）

* FD 历史最久、功能最多
* **应力张量（stress tensor）**与某些**响应函数（response functions）**只在 PW 模式可用

---

## 三种模式的功能支持表（Features）

| 功能/关键词                                         | FD                       | LCAO | PW  |
| ---------------------------------------------- | ------------------------ | ---- | --- |
| GPU ground state calculations                  | experimental             |      | Yes |
| Time-propagation TDDFT                         | YES                      | YES  |     |
| Dielectric function（介电函数）                      |                          |      | YES |
| Casida equation                                | YES                      |      |     |
| Hybrid functionals                             | （no forces, no k-points） |      | YES |
| Stress tensor（应力张量）                            |                          |      | YES |
| GW                                             |                          |      | YES |
| BSE                                            |                          |      | YES |
| Direct orbital optimization（广义 mode-following） |                          | YES  |     |
| Non-collinear spin（非共线自旋）                      |                          |      | YES |
| Solvent models（溶剂模型）                           | YES                      | YES  |     |
| MGGA                                           | YES                      |      | YES |
| Constrained DFT                                | YES                      | YES  |     |
| Ehrenfest                                      | YES                      |      |     |
| Spin-spirals                                   |                          |      | YES |

> 注：表中 YES/experimental 等表述保持原意。不同 GPAW 版本对功能覆盖可能会变化，遇到功能不可用或行为不同请以当前版本文档/报错信息为准。

---

# 关键参数详解与实用建议

下面按材料计算中最常“踩坑/最常需要手动控制”的参数来解释。

---

## `nbands`：电子能带数（Number of electronic bands）

`nbands` 决定每个自旋通道（spin channel）计算多少条 band。

* 对于**非自旋极化**体系，若价电子数为 10，则占据态为 5 条（每条 band 2 个电子），所以 `nbands=5` 可包含全部占据态。
* 对于**自旋极化**体系，spin-up 与 spin-down 占据数不同：例如 10 个价电子、总磁矩为 2 的体系，至少需要 `nbands=6`（示例解释：up 需要 6 条占据，down 4 条占据 + 2 条空带）。

默认规则（原文）：
默认 `nbands` = `4 + 1.2 ×（占据带数）`。

实用理解：

* **绝缘体/分子**：占据态与非占据态间隙大时，`nbands` 取“刚够放下占据态”通常就能跑。
* **金属**：需要更多空带来描述费米面附近部分占据，`nbands` 往往要更大；增加空带也可能改善收敛。

> **Tip**
>
> * `nbands=0`：不加空带（zero empty bands）
> * `nbands=-n`：给 n 条空带
> * `nbands='n%'`：给占据带数的 `n/100` 倍
> * `nbands='nao'`：band 数等于原子轨道数（atomic orbitals），对应 LCAO 可用的最大 `nbands`

---

## `xc`：交换-相关泛函（Exchange-Correlation functional）

常用 XC 泛函（保留原表与关键字）：

| `xc`       | full libxc keyword        | 说明                          | reference |
| ---------- | ------------------------- | --------------------------- | --------- |
| `'LDA'`    | `'LDA_X+LDA_C_PW'`        | Local density approximation | 1         |
| `'PBE'`    | `'GGA_X_PBE+GGA_C_PBE'`   | Perdew, Burke, Ernzerhof    | 2         |
| `'revPBE'` | `'GGA_X_PBE_R+GGA_C_PBE'` | revised PBE                 | 3         |
| `'RPBE'`   | `'GGA_X_RPBE+GGA_C_PBE'`  | revised revPBE              | 4         |
| `'PBE0'`   | `'HYB_GGA_XC_PBEH'`       | Known as PBE0               | 5         |
| `'B3LYP'`  | `'HYB_GGA_XC_B3LYP'`      | B3LYP (as in Gaussian Inc.) | 6         |

要点：

* 默认 `xc='LDA'`
* `PBE/revPBE/RPBE` 属于 **GGA**
* `PBE0/B3LYP` 属于 **hybrid functional（杂化泛函）**

GPAW 默认使用 **libxc** 提供的泛函（但 LDA/PBE/revPBE/RPBE/PW91 等可能用 GPAW 自己的实现）。libxc 的关键字来自 `xc_funcs.h`，通常去掉前缀 `XC_`，并用 `exchange + correlation` 用 `+` 拼接。

例如化学里常见的 LDA 近似可以写成：`'LDA_X+LDA_C_VWN'`。

### 用字典形式指定 XC（带参数的泛函）

某些泛函需要额外参数，可以用 dict：

```python
xc = {'name': 'PBE', 'stencil': 2}
```

这里 `stencil` 控制密度梯度的差分模板（适用于 GGA/MGGA）。

### Hybrid 泛函与 `--exact-exchange` setups

Hybrid（精确交换，Exact exchange / EXX）需要 setup 中包含 exx 信息。可以检查 setup 文件是否包含 `<exact_exchange_`：

```bash
$ zcat $GPAW_SETUP_PATH/O.PBE.gz | grep "<exact_exchange_"
```

对缺少 exx 信息的元素生成 setup：

```bash
$ gpaw-setup --exact-exchange -f PBE H C
```

目前所有 hybrid 泛函以 PBE setup 为基础。

---

## `kpts`：布里渊区采样（Brillouin-zone sampling）

默认只用 **Γ-point（Gamma-point）**。仅用 Γ 点的好处是波函数可以选为实数（real）。

### Monkhorst-Pack 网格

最常见写法：

```python
kpts = (N1, N2, N3)
```

即 `N1 × N2 × N3` 的规则 k 点网格。

更灵活的 dict 形式（推荐）：

```python
kpts={'size': (4, 4, 4)}  # 4x4x4 Monkhorst-pack
kpts={'size': (4, 4, 4), 'gamma': True}  # shifted 4x4x4 Monkhorst-pack
```

### 用 k 点密度指定（points per Å^-1）

```python
kpts={'density': 2.5}  # 最小密度 2.5 points/Ang^-1
kpts={'density': 2.5, 'even': True}  # 向上取最近的偶数
kpts={'density': 2.5, 'gamma': True}  # 包含 gamma 点
```

密度计算关系（原文）：
`N * (a / 2pi)`
其中 `N` 是该方向的 k 点数，`a` 是晶胞在对应倒格矢方向的长度。

### 显式给出任意 k 点集合

```python
kpts=[(0, 0, -0.25), (0, 0, 0), (0, 0, 0.25), (0, 0, 0.5)]
```

这些是相对于倒格基矢的**缩放坐标（scaled coordinates）**。

---

## `spinpol`：自旋极化（Spin-polarized calculation）

默认行为：

* 如果原子带磁矩（magnetic moments），则会进行**自旋极化**计算；
* 否则进行**自旋成对**（spin-paired）计算。

可用 `spinpol=True` 强制开启自旋极化。

> 实用建议：对于可能有磁性的材料（过渡金属、含缺陷的氧化物、表面吸附等），建议你**显式写出 `spinpol=True/False`** 并在结构中合理设置初始磁矩（ASE 的 `atoms.set_initial_magnetic_moments([...])`），避免“默认行为”导致结果不一致。

---

## `gpts` 与 `h`：网格点数与网格间距（Grid）

网格质量会显著影响总能量与力的收敛。

* 更密的网格（更小的 `h`，或更大的 `gpts`）通常更准确。
* 经验上，多数元素 `h=0.2 Å` 可获得合理的能量收敛（原文给出的经验值）。

### 直接指定网格点数 `gpts`

```python
gpts=(n1, n2, n3)
```

其中 `n1/n2/n3` 是正整数，并且**都要能被 4 整除**（原文要求）。

### 用 `h` 让程序自动选取 `gpts`

```python
h=0.25
```

程序会尽量选择使网格点密度约为 `1/h^3` 的 `gpts`（注意：要匹配晶胞尺寸，`h` 往往是近似值）。

### 网格间距与平面波 cutoff 的换算（了解即可）

Briggs 等给出过从平面波截断能到实空间网格间距的换算关系（PRB 54, 14362 (1996)）。在 GPAW 中可这样转换：

```python
>>> from gpaw.utilities.tools import cutoff2gridspacing, gridspacing2cutoff
>>> from ase.units import Rydberg
>>> h = cutoff2gridspacing(50 * Rydberg)
```

---

## 想要精确控制 `h`：`adjust_cell`

由于 `h` 需要与晶胞“整除匹配”，你指定的 `h` 常常是近似值。若你确实希望**精确的 `h`**，你需要相应调整晶胞。示例：

```python
from gpaw.utilities.adjust_cell import adjust_cell
from ase import Atoms

d = 0.74
a = 6.0
atoms = Atoms('H2', positions=[(0, 0, 0), (0, 0, d)])
# set the amount of vacuum at least to 4 Å
# and ensure a grid spacing of h=0.2
adjust_cell(atoms, 4.0, h=0.2)
```

这里 `adjust_cell(atoms, 4.0, h=0.2)` 的含义是：至少留 4 Å 真空，并调整晶胞使网格间距满足 `h=0.2 Å`。

---

## `symmetry`：对称性使用（Use of symmetry）

默认：使用所有点群对称（point-group symmetry）与时间反演对称（time-reversal symmetry），将 k 点约化到不可约布里渊区（irreducible BZ）。

注意：如果你移动原子导致某些对称性被破坏，可能会报错。可通过关闭点群对称避免：

```python
symmetry={'point_group': False}
```

这会只保留时间反演对称（哈密顿量在 `k -> -k` 下不变）。

若你希望完全不做 k 点对称性约化（比如调试、能带计算等）：

```python
symmetry={'point_group': False, 'time_reversal': False}
```

或简写为：

```python
symmetry='off'
```

`symmetry` dict 的完整键如下：

| key             | default | 说明                        |
| --------------- | ------- | ------------------------- |
| `point_group`   | `True`  | 使用点群对称                    |
| `time_reversal` | `True`  | 使用时间反演对称                  |
| `symmorphic`    | `True`  | 只使用 symmorphic symmetries |
| `tolerance`     | `1e-7`  | 相对容差                      |

---

## 波函数初始化（Wave function initialization）

默认会用 **原子轨道线性组合（LCAO guess）** 作为初猜。如果你要求的 `nbands` 多于预计算原子轨道数量，剩余 band 会用随机数初始化。

> 实用建议：当你发现 SCF 很难收敛或不同运行不一致时，初始化方式与 `nbands`/`occupations`/`mixer` 的组合往往是关键。必要时可在可重复性要求较高的场景避免随机初始化（例如减少多余空带、或保证初始轨道数足够）。

---

## `occupations`：占据数与展宽（Smearing）

用 dict 控制占据数展宽，例如费米-狄拉克（Fermi-Dirac）：

```python
from gpaw import GPAW
calc = GPAW(...,
            occupations={'name': 'fermi-dirac',
                         'width': 0.05},
            ...)
```

对应分布（width = `k_B T`）：

$$f(E) = 1 / (1 + e^{E/(k_B T)})$$

默认值（原文）：

* **周期体系（periodic boundary conditions）**默认 `width=0.1 eV`，并把总能量外推到 `T=0 K`。
* **分子（无周期边界）**默认 `width=0`，得到整数占据。

其它分布函数（原文列举）：

* `{'name': 'marzari-vanderbilt', 'width': ...}`
* `{'name': 'methfessel-paxton', 'width': ..., 'order': ...}`

### 固定总磁矩（spin-polarized）

自旋极化计算中，可以把总磁矩固定在初始值：

```python
occupations={'name': ..., 'width': ..., 'fixmagmom': True}
```

### 固定占据数（Fixed occupations）

用 `gpaw.occupations.FixedOccupationNumbers`：

```python
from gpaw.occupations import FixedOccupationNumbers
calc = GPAW(...,
            occupations=FixedOccupationNumbers([[1, 1, ..., 0, 0],
                                                [1, 1, ..., 0, 0]]))
```

> 实用建议：
>
> * 金属体系：通常需要非零 `width`（如 0.05–0.2 eV，具体取值需收敛测试）。
> * 绝缘体/分子：可用更小的 `width` 或 0（但对收敛困难体系仍可用小展宽辅助收敛）。
> * 做能带/态密度后处理时，建议记录你使用的 `occupations` 参数，方便结果复现。

---

## Compensation charges（补偿电荷）

补偿电荷会展开到正确的多极矩（multipoles）到 `l = lmax`，默认 `lmax=2`。

---

## `charge`：体系总电荷

默认电中性。单位是“负电子电荷”为 1 的计数方式（原文表述），例如：

* `charge = -1` 表示比中性多 1 个电子（整体带 -1 电荷）

---

## `convergence`：SCF 收敛标准（Accuracy of the self-consistency cycle）

`convergence` 用 dict 指定 SCF 收敛阈值。默认值（原文）：

```python
{'energy': 0.0005,  # eV / electron
 'density': 1.0e-4,  # electrons / electron
 'eigenstates': 4.0e-8,  # eV^2 / electron
 'bands': 'occupied'}
```

含义解释：

* `energy`：最近 3 次迭代的能量变化 < 0.5 meV / 价电子
* `density`：电荷密度变化（绝对值积分） < 0.0001 / 价电子
* `eigenstates`：Kohn–Sham 方程残差平方积分 < `4.0e-8 eV^2` / 价电子（**不影响 LCAO**）
* `bands='occupied'`：只要求占据态收敛

若你只给部分 dict，其余使用默认值。例如：

```python
convergence={'energy': 0.0001}
```

会把能量阈值改为 0.1 meV，其它保持默认。

还可以设置额外关键字，如 `'eigenvalues'`, `'forces'`, `'work function'`, `'minimum iterations'` 等（其中 `'eigenvalues'` 只在 “New GPAW” 实现）。你也可以自定义收敛判据（custom criteria）。

关于 `bands` 的高级用法（原文要点）：

* `{'bands': 'all'}`：连空带也强制收敛
* `{'bands': 200}`：收敛最低 200 条 band
* `{'bands': -10}`：除最后 10 条以外全部 band 收敛（最后几条空带往往最难收敛）
* `{'bands': 'CBM+5.0'}`：收敛到导带底（CBM）以上 5.0 eV（对金属 CBM 取费米能级）

---

## `maxiter`：最大 SCF 迭代步数

如果在 `maxiter` 次自洽迭代内仍未收敛，计算会报错停止。若你需要“至少迭代多少步”，可以通过自定义收敛准则实现。

---

## `txt`：文本输出位置（Where to send text output）

* 默认 `txt='-'`：输出到标准输出（stdout）
* `txt='filename.txt'`：输出到文件
* `txt=None`：关闭所有文本输出
* 也可以传入一个 file-like 对象（只要有 `write()` 方法）

> 实用建议：做可复现的材料计算，建议总是写 `txt='xxx.txt'`，并和 `*.gpw` 一起保存，便于追溯计算设置与收敛过程。

---

## `mixer`：密度混合（Density mixing / Pulay mixing）

GPAW 使用 Pulay mixing 混合电荷密度，核心参数：

* `beta`：线性混合系数
* `nmaxold`：参与混合的历史密度数量
* `weight`：衡量密度变化时，对长波长分量加权（long wavelength 更重要）

原文给出的经验选择：

* 小分子（零边界）：`Mixer(beta=0.25, nmaxold=3, weight=1.0)`（GPAW 在零边界时也会倾向这么选）
* 大分子/团簇/周期体系：`Mixer(beta=0.05, nmaxold=5, weight=50.0)`（周期边界时也会倾向这么选）

自旋极化计算会用 `MixerDif` 替代 `Mixer`。

> 实用建议：收敛困难（尤其金属、强极化、真空 slab）时，`mixer` 往往是最有效的“调参入口”之一。一般做法是减小 `beta`、增加 `nmaxold`、提高 `weight`。

---

## 固定密度计算：`fixed_density()`

在做能带结构（band structure）或增加空带并希望它们收敛时，常常希望**复用已收敛密度**，避免再做密度更新。可用：

* `gpaw.calculator.GPAW.fixed_density()`

这会在整个 SCF 过程中使用固定密度（类似 Harris 计算），密度可以来自 `.gpw` 文件或前一次计算。

---

## `setups`：PAW 数据集 / 赝势（PAW datasets or pseudopotentials）

`setups` 用于指定每个元素使用哪个 setup 文件。

对元素 `E`、setup 名 `NAME`、XC 名 `XC`，GPAW 会在 setup 路径里按顺序寻找：

* `E.NAME.XC` 或 `E.NAME.XC.gz`
* 若 `NAME='paw'`（默认），则寻找 `E.XC` 或 `E.XC.gz`

`setups` 可以是：

* 单个字符串：对所有原子使用同一 setup 名
* dict：可对不同元素或不同原子编号指定不同 setup

特殊 key：

* `'default'`：指定默认 setup 名
  `setups={'default': 'paw'}` 等价于 `setups='paw'`

例：Na 最新 PAW setup 可能把 6 个 semicore p 态也放进价电子。若要用只含 1 个 s 价电子的非默认 setup（`Na.1.XC.gz`），可：

```python
setups={'Na': '1'}
```

### 三个/四个特殊名字（不对应文件名）

* `'ae'`：全电子（all-electron）模式，不用 PAW/赝势
* `'sg15'`：SG15 优化的守恒范数 Vanderbilt 赝势（PBE），需要额外安装

  * 下载：`gpaw install-data --sg15 {<dir>}`
  * 查看/绘图：`gpaw-upfplot <pseudopotential>`
  * 注意：目前在 GPAW 中仍建议视为 experimental（原文提示）
* `'hgh'`：Hartwigsen-Goedecker-Hutter 守恒范数赝势（无需安装）

  * 某些元素有更好的 semicore 版本：用 `'hgh.sc'`
* `'ghost'`：LCAO 模式中的 ghost atom（用于 BSSE 等，见 Ghost atoms）

若 dict 同时对元素符号与原子编号都指定，**原子编号优先**。例如：

```python
setups={'default': 'soft', 'Li': 'hard', 5: 'ghost', 'H': 'ae'}
```

含义：

* Li 用 `hard`
* H 用全电子 `ae`
* 原子序号 5 是 ghost（无论它是什么元素）
* 其它原子用 `soft`

---

## `basis`：原子基组（Atomic basis set）

`basis` 主要用于 **LCAO 模式**，也会影响 FD/PW 模式下的 LCAO 初始猜测（先在 LCAO 基里解 Kohn–Sham 得到初始波函数）。

* `basis='basisname'` 时，GPAW 会在 setup 路径里找 `symbol.basisname.basis` 文件

  * 若元素用了非默认 setup 名，则会找 `symbol.setupname.basisname.basis`

dict 形式可对不同元素/原子指定不同基组：

```python
basis={'H': 'sz', 'C': 'dz', 7: 'dzp'}
```

表示：H 用 `sz`，C 用 `dz`，Atoms 对象中编号 7 的原子用 `dzp`。

> **Note**
> 如果你想从 `dzp` 基组里只取 `sz` 子集，可以写：`basis='sz(dzp)'`。
> 如果基组文件有自定义名：`'szp(mybasis.dzp)'`。

`basis=None`（默认）表示使用 setup 中的 pseudo partial waves 作为基组（总是可用）。其它基组要求对应 `.basis` 文件存在于 setup 路径。

---

## `eigensolver`：本征求解器（Eigensolver）

默认：`eigensolver='dav'`（Davidson）通常表现良好。

其它选项：

* `eigensolver='rmm-diis'`：RMM-DIIS（少量空带时往往效率高）
* `eigensolver='cg'`：共轭梯度（稳定但更慢）

若需要按 band 并行，必须用 Davidson 或 RMM-DIIS。

也可以直接传 eigensolver 对象以获得更细控制：

```python
from gpaw.eigensolvers import CG
calc = GPAW(..., eigensolver=CG(niter=5, rtol=0.20), ...)
```

* `niter`：每条 band 在一次 SCF step 内最多 CG 迭代次数
* `rtol`：残差相对变化小于该值则提前停止该 band 的 CG

LCAO 有自己的 eigensolver。`DirectLCAO` 会直接对哈密顿矩阵对角化。也可以用 ETDM（Exponential Transformation Direct Minimization），但**不推荐用于金属**（因为 ETDM 中占据数不是变分地确定的）。

---

## `poissonsolver`：Poisson 求解器与偶极修正（Dipole correction）

`poissonsolver` 可指定 Poisson solver 类，或启用 dipole-layer 修正。

默认行为（原文）：

* FD/LCAO 模式默认 `FastPoissonSolver`：结合 Fourier/Fourier-sine 变换 + 并行转置
* PW 模式：在倒空间直接除以 `|G|^2`

### 使用旧的 multigrid Jacobian Poisson solver（示例）

```python
from gpaw import GPAW, PoissonSolver
calc = GPAW(...,
            poissonsolver=PoissonSolver(
                name='fd', nn=3, relax='J', eps=2e-10),
            ...)
```

* `nn`：差分模板范围（stencil）
* `relax`：`'J'`（Jacobian）或 `'GS'`（Gauss-Seidel）
* `eps`：收敛阈值

> **Note**
> Poisson solver 通常不是性能瓶颈，但在某些网格布局下可能表现差（LCAO 更常见）。了解这一点有助于排查异常慢或异常不收敛的情况。

### 偶极层修正（dipole-layer correction）

当体系在某一方向**非周期**，而另外两个方向周期（典型：表面 slab，真空方向非周期）时，可用 dipole 修正：

```python
from gpaw import GPAW

correction = {'dipolelayer': 'xy'}
calc = GPAW(..., poissonsolver=correction, ...)
```

这里 `'xy'` 表示偶极层平行于 xy 平面（即沿 z 方向修正）。
无修正时，势在所有非周期边界趋向 0；有修正时，可在体系两侧引入与偶极矩相关的势差。

还可以同时传给 Poisson solver 其它参数：

```python
GPAW(...,
     poissonsolver={'dipolelayer': 'xy', 'name': 'fd', 'relax': 'GS'},
     ...)
```

### FFT Poisson solver（完全周期体系）

```python
GPAW(..., poissonsolver={'name': 'fft'}, ...)
```

FFT Poisson solver 一般对网格更“不挑剔”，能减少对网格间距的依赖。也可用于非周期体系，但要求你把体系显式设为周期并在非周期方向放足够真空，避免跨边界的非物理相互作用。

---

## 有限差分模板与插值（Finite-difference stencils）

GPAW 可为 Laplacian 使用不同范围的有限差分模板（stencil）。

### Poisson 方程的 stencil

```python
from gpaw import GPAW, PoissonSolver
calc = GPAW(..., poissonsolver=PoissonSolver(nn=n), ...)
```

误差阶：$O(h^{2n})$，`n` 取 1–6，默认 `n=3`。

### Kohn–Sham 方程的 stencil（FD 模式）

```python
from gpaw import GPAW, FD
calc = GPAW(mode=FD(nn=n))
```

默认也 `n=3`。

### 密度插值（FD/LCAO）

PW 模式中粗网格到细网格用 FFT；FD/LCAO 用 tri-quintic 插值（5 次多项式）。可设置插值阶数：

```python
from gpaw import GPAW, FD
calc = GPAW(mode=FD(interpolation=n))
# or
from gpaw import GPAW, LCAO
calc = GPAW(mode=LCAO(interpolation=n))
```

多项式阶数为 `2n-1`，默认 `n=3`，允许 1–4（线性、三次、五次、七次）。

---

## `hund=True`：用 Hund’s rule 猜初始磁矩

设置 `hund=True` 会：

* 强制自旋极化（spinpolarized）
* 初始占据与各原子磁矩按 Hund’s rule 设定
* 用户显式指定的磁矩会被忽略（原文提示）

如希望固定总磁矩，可结合：

```python
occupations={'name': ..., 'fixmagmom': True}
```

默认 `hund=False`。

---

## `external`：外加势（External potential）

例如施加恒定电场：

```python
from gpaw.external import ConstantElectricField
calc = GPAW(..., external=ConstantElectricField(2.0, [1, 0, 0]), ...)
```

更多见 `gpaw.external`。

---

## 输出详细程度（Output verbosity）

默认每个 SCF step 只输出有限信息。若要更详细信息（尤其用于诊断收敛问题），可用 `verbose=1`。

---

## `communicator`：通信器（Communicator object）

可指定 communicator，让 calculator 只使用一部分 MPI 进程，用于并行跑多个计算（例如不同原子构型 image）。见“Running different calculations in parallel”。

---

# 总能量的含义（Total Energies）

GPAW 计算的能量是相对于“分离参考原子”的能量：参考原子是**自旋成对、带电中性、球对称**的状态——也就是生成 setup 时使用的原子状态。

因此：

* 分子：计算得到的能量 ≈ 负的原子化能（atomization energy）
* 固体：计算得到的能量 ≈ 负的内聚能（cohesive energy）

如果你得到**正的总能量**，意味着你的体系相对于参考原子是**不稳定态**（原文强调）。

> **Note**
> 得到的并不是“真实的”原子化/内聚能。真实值通常更低，因为许多原子的真实基态是**自旋极化且非球对称**，能量低于该参考原子。

---

# 保存与重启计算（Restarting a calculation）

你可以把计算状态写入 `.gpw` 文件：

```python
calc.write('H2.gpw')
```

`H2.gpw` 是二进制文件，包含波函数、密度、原子位置以及计算器参数等。

之后可在另一个 Python 会话中重启：

```python
>>> from gpaw import *
>>> atoms, calc = restart('H2.gpw')
>>> print(atoms.get_potential_energy())
```

重启后你可以小幅修改参数。例如改变网格点数：

```python
atoms, calc = restart('H2.gpw', gpts=(20, 20, 20))
print(atoms.get_potential_energy())
```

> **Tip（另一种常用写法）**
>
> ```python
> atoms, calc = restart('H2.gpw')
> calc.set(gpts=(20, 20, 20))
> print(atoms.get_potential_energy())
> ```

更多细节见 Restart files。

---

# 通过 observer 自定义计算行为（attach / observers）

你可以给 calculator 绑定一个 observer，使其每 N 次迭代执行一次。示例：每 5 次迭代输出一个不同名字的 restart 文件。

```python
calc = GPAW(...)

occasionally = 5

class OccasionalWriter:
    def __init__(self):
        self.iter = 0

    def write(self):
        calc.write('filename.%03d.gpw' % self.iter)
        self.iter += occasionally

calc.attach(OccasionalWriter().write, occasionally)
```

另见 `attach()`。

> 实用建议：
>
> * 长时间计算（大体系、难收敛）中，定期写 `.gpw` 能显著降低意外中断造成的损失。
> * 也可以 attach 自定义的收敛监控、打印关键量、输出轨迹等。

---

# 调试模式（Debug mode）

## `GPAW_DEBUG`

以 debug-mode 运行 GPAW（例如检查传入 C 扩展的数组一致性等）。

---

# 命令行选项（Command-line options）

如果你通过 `gpaw python` 运行脚本，可以使用 dry-run 模式：

```bash
$ gpaw python --dry-run=N script.py
```

dry-run 会：

* 打印计算参数
* 估算内存
* **不执行实际计算**
* 还会打印在 `N` 核下将采用的并行设置

> Tip（开发/调试用）
> 如果你需要从命令行额外传入一些参数给 Python（用于开发调试）：
>
> ```bash
> $ python3 -X a=1 -X b
> >>> import sys
> >>> sys._xoptions
> {'a': '1', 'b': True}
> ```
>
> 见 Python 的 `-X` 选项。

---

# 关于持续重构与 “New GPAW”

某些功能可能只在 “new” 或 “old” GPAW 中可用。这源于 GPAW 正在进行的计算后端重构（refactoring）。更多信息见 `“New GPAW”` 相关页面/说明。

---

#（附）一个更偏“材料计算”的常用脚本骨架（建议收藏）

下面给一个更符合周期材料计算习惯的“骨架”，便于你快速起步（参数需要按体系做收敛测试）：

```python
from ase.io import read
from gpaw import GPAW, PW

atoms = read('POSCAR')          # 或 CIF/XYZ 等
atoms.pbc = True                # 周期体系
# atoms.center(vacuum=10.0, axis=2)  # slab 时可加真空并居中（示意）

calc = GPAW(
    mode=PW(340),               # 或 PW(400/500...)，需收敛测试
    xc='PBE',
    kpts={'density': 2.5, 'gamma': True},
    occupations={'name': 'fermi-dirac', 'width': 0.1},  # 金属常用；绝缘体可更小
    convergence={'energy': 1e-4}, # 仅示例，按需求调整
    txt='gpaw.out'
)

atoms.calc = calc
energy = atoms.get_potential_energy()
print('E =', energy, 'eV')
calc.write('state.gpw')
```

> 关键提醒：
>
> * **收敛测试（convergence test）**仍然是材料计算最关键的一步：k 点密度、PW cutoff 或 FD 网格 `h`、展宽 `width`、`nbands` 都可能影响结果。
> * 对表面 slab 的真空方向与偶极修正（dipole correction）也常常很重要。

# Convergence Issues（收敛问题排查与处理指南）

> 核心原则：**先尽量使用 GPAW 的默认参数**（default parameters）。默认往往“简单且够用”，也更不容易把体系带进奇怪的数值状态。只有在出现收敛问题（SCF 不收敛、能量/力振荡、密度发散等）时，再逐项调整。

下面给出一份**按优先级排序**的实用建议清单。通常从第 1～4 条就能解决大多数收敛问题。

---

## 1. 先确认结构与自旋态是否合理（Geometry & spin-state）

* **单位检查：ASE 默认长度单位是 Ångström（Å）**，不是 Bohr，也不是 nm。
  很多“怎么都不收敛”的根源其实是结构尺度错了（例如把 nm 当 Å）。
* 自旋极化（spin-polarized）体系：给出合理的初始磁矩（initial magnetic moments）。
  例如过渡金属、含缺陷体系、开壳层分子等。
* **奇数电子数的分子不要做 spin-paired（非自旋极化）**计算，否则物理上不一致，SCF 往往会挣扎或给出错误结果。
* 做“孤立原子”相关计算前，建议先阅读 *Calculation of atomization energies*（原子化能计算）相关说明，避免参考原子态与设置不匹配。

> **实用检查清单（建议你每次卡收敛都先看一眼）**
>
> * 结构是否有严重重叠/过近键长（<0.6 Å 之类）？
> * 真空层是否足够（分子/表面 slab）？
> * `atoms.pbc` 是否与你的体系一致？
> * 自旋是否合理：奇数电子？磁性元素？初始磁矩是否全 0？

---

## 2. 降低密度混合强度（Use less aggressive density mixing）

收敛振荡/发散时，**最常用也最有效**的手段之一是让密度混合（density mixing / Pulay mixing）更“保守”。

可以尝试：

* 非自旋体系：`Mixer` 或 `MixerSum`
* 自旋极化体系：`MixerDif`

示例（把混合系数调小、增大长波权重）：

```python
mixer=Mixer(0.02, 5, 100)
# 或
mixer=MixerSum(0.02, 5, 100)
# 自旋极化：
mixer=MixerDif(0.02, 5, 100)
```

记得导入 mixer 类：

```python
from gpaw import Mixer, MixerSum, MixerDif
```

> **额外经验**
>
> * 某些体系（例如过渡金属原子、强局域磁性体系）可能对“历史密度”很敏感，这时把 mixer 的历史步数从 `5` 降到 `1` 往往有帮助（更像线性混合，减少“过度纠正”）。

---

## 3. 每步 SCF 更精确地解本征问题（Solve eigenproblem more accurately）

当 SCF “能量在来回跳”或“密度几乎收敛但最后几步不稳定”时，可能是每个 SCF step 里本征值问题（eigenvalue problem）解得不够充分。

你可以：

1. 用 Davidson eigensolver，并增加每步迭代次数：

```python
from gpaw import Davidson
```

例如：

```python
eigensolver=Davidson(3)
```

2. 或者试试共轭梯度（CG）：

```python
eigensolver='cg'
```

> 原文提示：CG 往往对**空带（unoccupied bands）**收敛更快。
> **实用理解**：如果你加了很多空带（例如做光学/激发态/介电函数前准备），CG 有时能减少“空带拖累 SCF”的问题。

---

## 4. 使用更平滑的占据数分布（Use smoother occupations）

* 对**非周期体系（分子）**，GPAW 默认 `Fermi temperature = 0`（等价于整数占据），这在某些体系会让 SCF 很“硬”，容易卡。
* 可以设置一个有限的展宽（smearing），如 `fermi-dirac`，并**检查结果对温度/展宽的收敛性**（很重要！）

（占据数设置详见 *Occupation numbers* 部分。）

> **实用建议**
>
> * 金属：通常必须用有限展宽（例如 0.05–0.2 eV，需做收敛测试）。
> * 绝缘体/分子：可以用更小展宽（例如 0.01–0.05 eV）辅助收敛，然后再逐步降低验证能量/结构是否稳定。

---

## 5. 增加空带数量（Try adding more empty states）

如果你手动设置了 `nbands`，而 SCF 不稳定，可以尝试：

* 增加空带（empty states）
* 或者干脆让 GPAW 使用默认的 `nbands`（一般足够大）

> **为什么有用？**
> 对金属、半金属、或费米能附近态密集的体系，空带太少会让占据/密度更新不稳定。

---

## 6. 用足够的 k 点（Use enough k-points）

周期体系收敛困难时，k 点采样不足也可能导致异常振荡或不稳定（尤其金属）。

可尝试用 k 点密度：

```python
kpts={'density': 3.5, 'even': True}
```

（见 *Brillouin-zone sampling*）

---

## 7. 结构优化不要一步走太大（Optimization step size）

如果你在做几何优化（geometry optimization），**优化器步长过大**可能会把结构带到数值上很差的位置，导致下一步 SCF 崩溃。

> **实用建议**
>
> * 减小 optimizer 的最大步长（max step）
> * 或先用更“粗”的设置把结构拉回合理区域，再逐步提高精度

---

## 8. 改善波函数初猜（Better initial guess for wave functions）

GPAW 的波函数初猜默认总是走 LCAO（linear combination of atomic orbitals）方案，并使用默认的 **single-zeta（`sz`）** 基组（每个价电子一个轨道）。

你可以尝试用更好的初猜基组，例如从 `dzp`（double-zeta polarization）中提取 `szp`：

```python
basis='szp(dzp)'
```

也可以通过增大 `Atomic basis set`（原子基组）来改善初猜。但要注意：

* 你需要先生成对应的 basis 文件（见 *LCAO mode* 的说明）
* **Warning：某些情况下这反而会让收敛变差**
  原文强调：通常只有当你同时显著增加空带数（empty states）时，这才更可能改善收敛。

---

# Restart files（重启文件 .gpw）

`.gpw` 文件是 GPAW 的重启文件（restart file），用于保存一个计算的电子结构状态，便于中断后恢复、或在已有密度基础上继续做后续计算（加空带、算能带、算响应等）。

---

## 写入重启文件（Writing restart files）

使用 `write()` 方法：

```python
calc.write('xyz.gpw')
```

这会保存：密度（density）、势（potential）、本征值（eigenvalues）以及其它必要信息。

### 同时保存波函数（wavefunctions）

波函数可能非常大。若你确实需要保存波函数（例如后处理需要），用：

```python
calc.write('xyz.gpw', mode='all')
```

---

## “New GPAW” 下更小的 .gpw 文件（节省空间技巧）

在 **“New GPAW”** 中还提供了额外方式减小 `.gpw` 文件：

1. 不保存（可能很大）的 PAW projections：

```python
calc.write('xyz.gpw', include_projection=False)
```

2. 进一步用单精度保存（single-precision）：

```python
calc.write('xyz.gpw', precision='single')
```

> **Note**
>
> * 重新加载时会以双精度（double-precision）形式加载，因此后续处理仍是双精度。
> * 单精度重启文件适合“磁盘/IO 很紧张”的场景，但如果你在做非常敏感的能量差（meV 级别）比较，仍建议评估它是否影响可重复性。

---

## Tip：自动定期写重启文件（attach）

可以把 `write()` 注册成 observer，让它每隔 n 次 SCF 迭代自动写一次：

```python
calc.attach(calc.write, n, 'xyz.gpw')
```

或保存全部（含波函数）：

```python
calc.attach(calc.write, n, 'xyz.gpw', mode='all')
```

> 这对昂贵计算非常有用：
> 即使 SCF 还没跑完、作业被队列中断，你也能从“中间态”的电子结构继续跑下去。

---

## 读取重启文件（Reading restart files）

两种常用方式：

方式 1：直接用文件名构造 calculator：

```python
calc = GPAW('xyz.gpw')
```

方式 2：同时得到 `atoms` 与 `calc`：

```python
atoms, calc = restart('xyz.gpw')
```

> **实用建议**
>
> * 用 `.gpw` 继续计算前，养成习惯：在输出里记录你重启时修改了哪些参数（kpts、nbands、mode、xc 等），便于复现实验链路。
> * 做能带/态密度常见流程：
>   先做一次收敛的自洽（SCF）→ 写 `.gpw` → 固定密度（fixed density）加密 k 点或加空带做后处理。

---

# “New GPAW”（新后端说明）

GPAW 的计算后端正在进行较大规模的重构（refactoring）。因此文档中有时会区分：

* **“new” GPAW**（新后端）
* **“old” GPAW**（旧后端）

## 如何显式选择新/旧后端

使用新后端创建 calculator：

```python
from gpaw.new.ase_interface import GPAW as NewGPAW
```

显式使用旧后端：

```python
from gpaw.calculator import GPAW as OldGPAW
```

默认行为：除非设置环境变量 `GPAW_NEW`，否则默认使用 **old GPAW**。

> **实用提示**
>
> * 在脚本/项目里建议明确写清你使用的是 `NewGPAW` 还是 `OldGPAW`，这样别人复现时不会因为环境变量不同而踩坑。
> * 如果你在集群上运行，环境变量可能由模块系统或作业脚本设置，出问题时优先检查 `echo $GPAW_NEW`。

---

# D3 Correction（Grimme DFT-D3 色散修正）

`D3` 扩展（extension）使用 Stefan Grimme 的 **DFT-D3** 方法，为常见 XC 泛函补充**长程色散/范德华（Van der Waals, vdW）**相互作用。

为什么重要？
大多数常规泛函（LDA/GGA 等）对长程 vdW 描述不足，因此在以下体系中 D3 往往会带来显著改进：

* 层状材料（2D）与其堆垛（bilayer / multilayer）
* 分子吸附在表面（physisorption）
* 分子晶体、弱相互作用复合物
* 多孔材料中客体-框架弱相互作用

---

## 用法（Usage）

D3 扩展有一个必需关键字：`xc`（用于指明与之匹配的 exchange-correlation functional）。

下面示例展示了一个 2D 双层体系（bilayer）中 D3 的用法（原文示例保留，加入少量注释）：

```python
from gpaw.new.ase_interface import GPAW
from gpaw.new.extensions import D3
from gpaw import PW
from ase.build import mx2

MoS2 = mx2('MoS2', a=3.2)
WSe2 = mx2('WSe2', a=3.2)

# 6.6Å of distance between layers
MoS2.positions[:, 2] += 3.3
WSe2.positions[:, 2] -= 3.3

bilayer = WSe2 + MoS2
bilayer.center(vacuum=6.0, axis=2)
bilayer.pbc = True  # suppress the D3 warning

calc = GPAW(mode=PW(400),
            xc='PBE',
            extensions=[D3(xc='PBE')],
            txt='out.txt')

bilayer.calc = calc
energy = bilayer.get_potential_energy()

# Acces the D3 correction energy
d3_correction = bilayer.calc.dft.d3.get_energy()
```

> **Note（周期性限制非常重要）**
> D3 要求体系要么：
>
> * **完全非周期**（not periodic），要么
> * **三维全周期**（periodic in all three directions）
>
> 如果体系只在 1 个或 2 个方向周期（典型 2D slab：x,y 周期，z 非周期），会出现 warning，并且 D3 会**假设体系三维周期**。

### 对 2D / slab 的实用建议（扩展说明）

由于 D3 在 2D slab 情况下会当作 3D 周期处理，你需要特别注意：

* 沿真空方向（通常 z）要留足够大 vacuum，以减小周期镜像之间的 **D3 虚假相互作用**。
* 做 **vacuum 收敛测试**：逐步增大真空层，看 D3 修正能量是否变化到可忽略（例如 < 1 meV/原子 或 < 1 meV/Å²，按你研究目标选择标准）。

---

## 输出与读取 D3 能量（Output）

D3 修正能量可通过 d3 对象的 `get_energy()` 读取（如上例最后一行）。

同时它也会写入输出文件（在迭代列表之后，extensions 区域）。例如上例的 `out.txt` 中相关片段：

```text
...
Energy contributions relative to reference atoms: (reference = -702828.048343)

kinetic:          71.457107
coulomb:         -50.449930
zero:             -0.629602
external:          0.000000
xc:              -44.310857
entropy:          -0.001037
spinorbit:         0.000000
--------extensions:---------
D3 (xc=PBE):      -1.446148
----------------------------
Free energy:     -25.380467
Extrapolated:    -25.379948
...
```

> **实用建议**
>
> * 发表或汇报时，建议同时记录：`xc`、是否启用 D3、D3 修正能量的量级、以及真空/k 点/cutoff 的收敛设置。
> * 如果你在比较不同体系的相对稳定性（能量差），请确保所有体系使用一致的 D3 设置与周期/真空处理方式。

---

# Atomic PAW Setups（原子 PAW 数据集 / Setup）

**Setup（PAW dataset）** 在 PAW 方法中的角色，类似于赝势法（pseudopotential method）中的 pseudopotential：
它包含了某元素在特定 XC 泛函下的原子参考信息、投影子（projectors）、部分波（partial waves）等，是 GPAW 计算的基础数据之一。

* 所有可用 setups 都打包在 **tar 文件（tar-files）**中
* 提供多种泛函版本：LDA、PBE、revPBE、RPBE、GLLBSC 等
* setups 文件以压缩的 **XML specification for atomic PAW datasets files** 存储（即压缩 XML 规范文件）

安装方法见下面 *Installation of PAW datasets*。

---

## Setup releases（版本发布说明）

> **注意（很重要）**
> 下面版本信息来自你提供的原文片段（截至 2024 年）。GPAW setups 可能会在之后继续更新。
> **做严肃材料计算时，建议你记录自己使用的 setups 版本号，并在需要时去官方页面确认是否有更新。**

### 发布摘要表（便于快速查看）

| 日期（Date）   | 版本（Tarfile） | 主要变化（简述）                                         |
| ---------- | ----------- | ------------------------------------------------ |
| 2024-11-27 | `24.11.0`   | 新增 lanthanides（镧系）PAW potentials（含 PBE 与 LDA 版本） |
| 2024-02-22 | `24.1.0`    | 新增 Cr 的 14-electron PAW potential（高精度推荐）         |
| 2016-03-22 | `0.9.20000` | 旧版本（无额外说明）                                       |

### 24.11.0（2024-11-27）：镧系 PAW potentials

新增 lanthanides（La–Lu）PAW potentials。示例生成命令（La）：

```bash
$ gpaw dataset La -sw -r2.2 -P5s,6s,5p,6p,5d,d,4f,f,G -fPBE -b
```

其它元素（Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu）类似生成。也提供 LDA 版本 potentials。

### 24.1.0（2024-02-22）：Cr 的 14-electron PAW potential

新增 Cr 的 **14-electron** PAW potential。原文建议：若追求高精度，它比旧的 6-electron 版本更推荐（旧 6-electron 仍是默认）。

使用方式：

```python
setups={'Cr': '14'}
```

生成命令（原文）：

```bash
$ gpaw dataset Cr -sw -r2.0 -P3s,4s,3p,4p,3d,d,F -fPBE -t 14 -b
```

也提供 LDA 版本。

---

## Installation of PAW datasets（安装 PAW 数据集）

默认安装通常包含一个基础 PAW dataset；其他 datasets 可以**自动安装**或**手动安装**。

---

### 方式 A：自动安装（推荐）

命令形式：

```bash
gpaw install-data --{<dataset>} {<dir>}
```

作用：

* 下载并解包最新包到：`<dir>/<name>-<version>`
* 过程中会提示你是否把路径写入 GPAW 配置文件（建议输入 `y`）

> **实用建议**
>
> * 在多用户服务器上，常把 `<dir>` 放到共享路径（例如 `/usr/share/...` 或某个公共软件目录），然后统一配置 `GPAW_SETUP_PATH`。

---

### 方式 B：手动安装（更可控，适合离线/集群环境）

1. 从 *Atomic PAW Setups* 页面获取 tar 文件：`gpaw-setups-<version>.tar.gz`
   解压到某处（建议 `$HOME`，或全局共享路径如 `/usr/share/gpaw-setups/`）：

```bash
cd
tar -xf gpaw-setups-<version>.tar.gz
```

解压后会得到目录：

```text
gpaw-setups-<version>/
```

其中包含常用泛函的所有原子数据。

2. 设置环境变量 `GPAW_SETUP_PATH` 指向该目录。以 bash 为例，在 `~/.bashrc` 中加入：

```bash
export GPAW_SETUP_PATH=~/gpaw-setups-<version>
```

> **Note**
> 如果你配置了多个 PAW dataset 路径（多个位置），GPAW 会使用**第一个找到的 setup 文件**。
> 因此路径顺序非常关键：你可能以为在用新版本，实际却被旧版本“抢先匹配”了。

---

### `GPAW_SETUP_PATH` 的含义（原文保留并补充说明）

**GPAW_SETUP_PATH：**
用冒号 `:` 分隔的一组目录路径（colon-separated paths），每个目录都包含 PAW datasets（setup 文件）。

例如：

```bash
export GPAW_SETUP_PATH=/path/to/gpaw-setups-24.11.0:/path/to/other-setups
```

> **实用建议（可复现性）**
>
> * 论文/报告/项目记录中建议写明：`GPAW` 版本 + `gpaw-setups` 版本 + `xc` + `mode` + `kpts` + `cutoff/h`。
> * 尤其当你使用非默认 setups（如 `setups={'Cr': '14'}`）时，更应明确记录，否则结果很难复现或对比。

# 高级主题（Advanced topics）

本节翻译并整理 GPAW *Advanced topics* 中的前 4 个主题：

1. **Advanced Poisson solvers（高级 Poisson 求解器）**
2. **Delta Self-Consistent Field（ΔSCF / dSCF）**
3. **External potential（外加势/外场）**
4. **Grids（网格）**

在保留必要英文术语（类名、关键字、模块名、方程符号）的同时，对结构做了适度优化，并补充了面向材料计算的实用说明与注意事项。

---

## 1) Advanced Poisson solvers（高级 Poisson 求解器）

### 1.1 背景：为什么默认 PoissonSolver 会不够？

GPAW 的 `PoissonSolver` 在默认参数下，通常在晶胞边界施加 **zero boundary conditions（零边界条件）**。对于存在**大偶极矩（large dipole moment）**的体系，这会成为问题。例如：

* 金属纳米颗粒的等离激元（plasmonic）电荷振荡导致强偶极（或更高多极）分布；
* 强极化的纳米结构、团簇（cluster）、非对称大分子等。

由于偶极势是 **long-ranged（长程）** 的，想得到收敛的电势（尤其 Hartree potential），往往需要非常大的真空层（vacuum），否则边界条件会“影响”电势。

但在很多 **LCAO** 计算里，大真空并不是结构/电子本身所必需（只是 Poisson 边界导致你被迫加大真空）。因此，为了避免“为了电势而加很大真空”，GPAW 提供两类思路（也可组合使用）：

1. **Multipole moment corrections（多极矩修正）**
2. **在更大 Poisson 网格上解 Poisson 方程（extended grid / extra vacuum grid）**

它们分别由：

* `MomentCorrectionPoissonSolver`
* `ExtraVacuumPoissonSolver`

实现。

> **实用提示（来自原文的强提醒）**
> 在纳米颗粒等离激元（nano-particle plasmonics）计算中，**使用多极矩修正通常是必要的**。
> 若不修正，往往需要 **> 10 Å** 真空才能得到收敛结果（代价很高）。

---

### 1.2 方法 A：Multipole moment corrections（多极矩修正）

#### 核心思想

通过向密度加入一组“修正密度”，使得修正后的密度在边界上对应的**多极矩（multipole moments）**为零（或被抵消）。再把这些修正密度对应的势加回到总势中，从而改善边界条件导致的误差。

#### 基本用法

```python
from gpaw.poisson import PoissonSolver
from gpaw.poisson_moment import MomentCorrectionPoissonSolver

poissonsolver = MomentCorrectionPoissonSolver(
    poissonsolver=PoissonSolver(),
    moment_corrections=4
)
```

* `moment_corrections=4`：修正前 4 个多极矩：

  * **s**（单极/总电荷项）
  * **p_x, p_y, p_z**（偶极项的三个分量）

修正多极矩的数量可通过 `moment_corrections` 调整。例如：

* `moment_corrections=9`：在上述基础上额外包含 **d** 型多极矩：
  $d_{xx}, d_{xy}, d_{yy}, d_{yz}, d_{zz}$ 等（原文列举）

> **实用建议**
>
> * 类球形的金属纳米颗粒：`moment_corrections=4` 或 `9` 往往就够。
> * 更复杂几何（长棒、二聚体、强非对称）：可能需要更高阶多极矩，或使用**多中心（multicenter）**多极矩修正（见下）。

#### 高级语法：多中心多极矩（multicenter multipole）

原文指出，上述简写等价于：

```python
from gpaw.poisson import PoissonSolver
from gpaw.poisson_moment import MomentCorrectionPoissonSolver

poissonsolver = MomentCorrectionPoissonSolver(
    poissonsolver=PoissonSolver(),
    moment_corrections=[{'moms': range(4), 'center': None}]
)
```

这里：

* `moment_corrections` 是一个由字典组成的列表
* 字典关键字：

  * `moms`：要修正的多极矩集合

    * `range(4)` 对应 $s, p_x, p_y, p_z$
  * `center`：修正项的中心位置（原文写 **atomic units**）

    * `None` 表示晶胞中心

> **单位提醒（很关键）**
> 原文写 `center` 使用 **atomic units**（通常指 Bohr）。但示例里 `center` 直接用 `(x, y, z)` Å 写入。
> **不同 GPAW 版本接口可能不同**：建议你在实际使用时确认 `center` 的单位要求。
> 如确实需要 Å → Bohr，可参考（ASE 中 `Bohr ≈ 0.529177 Å`）：
>
> ```python
> from ase.units import Bohr
> center_bohr = center_angstrom / Bohr
> ```

#### 示例：纳米颗粒二聚体（dimer）的双中心修正

若两个纳米颗粒中心在 `(x1, y1, z1)` Å 和 `(x2, y2, z2)` Å：

```python
import numpy as np
from gpaw.poisson import PoissonSolver
from gpaw.poisson_moment import MomentCorrectionPoissonSolver

moms = range(4)
center1 = np.array([x1, y1, z1])
center2 = np.array([x2, y2, z2])

poissonsolver = MomentCorrectionPoissonSolver(
    poissonsolver=PoissonSolver(),
    moment_corrections=[
        {'moms': moms, 'center': center1},
        {'moms': moms, 'center': center2}
    ]
)
```

当使用多个 center 时，GPAW 会把计算晶胞划分为**不重叠区域**：空间中每个点会被分配给距离最近的 center（类似 **Voronoi diagrams（沃罗诺伊划分）** 的思想）。

#### 如何在 GPAW 计算器中启用？

将上面定义的 `poissonsolver` 传给计算器：

```python
from gpaw import GPAW
calc = GPAW(..., poissonsolver=poissonsolver, ...)
```

---

### 1.3 方法 B：在更大 Poisson 网格上求解（ExtraVacuumPoissonSolver）

#### 适用场景

多极矩修正对**复杂几何**不一定总是成功（例如非常不规则的形状、强多中心极化等）。此时可以采用另一条路：
**只为 Hartree 势（Poisson 方程）使用更大的网格**，但电子结构/波函数仍在原来的“小网格”上计算，从而减少“必须加超大真空”的压力。

#### 基本用法：指定一个“大网格”解 Poisson

```python
from gpaw.poisson import PoissonSolver
from gpaw.poisson_extravacuum import ExtraVacuumPoissonSolver

poissonsolver = ExtraVacuumPoissonSolver(
    gpts=(256, 256, 256),
    poissonsolver_large=PoissonSolver()
)
```

说明：

* `gpts`：大网格的网格点数（**单位是 Poisson 网格的单位**，通常与 fine grid 相同）
* `poissonsolver_large`：用于在大网格上解 Poisson 的求解器实例

> **性能提示（原文）**
> 如果你使用 `FDPoissonSolver`（多重网格 multigrid），网格尺寸最好能被较高次幂的 2 整除（例如 128, 256, 512 …），这样 multigrid 会更快。

#### 用 coarses 加速：先粗化大网格再求解

```python
from gpaw.poisson import PoissonSolver
from gpaw.poisson_extravacuum import ExtraVacuumPoissonSolver

poissonsolver = ExtraVacuumPoissonSolver(
    gpts=(256, 256, 256),
    poissonsolver_large=PoissonSolver(),
    coarses=1,
    poissonsolver_small=PoissonSolver()
)
```

* `coarses`：对给定的大网格粗化的次数

  * `coarses=1` 意味着大网格 `(256,256,256)` 会先粗化一次，实际求解网格变成 `(128,128,128)`
  * 网格间距变为原来的 2 倍
* `poissonsolver_small`：在原“小而细”的网格上求解势，并用“大网格（粗化）势”来修正边界条件

#### 组合/嵌套：逐级更大的网格

由于 `ExtraVacuumPoissonSolver` 是 wrapper（包装器），可以叠加：

```python
from gpaw.poisson import PoissonSolver
from gpaw.poisson_extravacuum import ExtraVacuumPoissonSolver

poissonsolver0 = ExtraVacuumPoissonSolver(
    gpts=(256, 256, 256),
    poissonsolver_large=PoissonSolver(),
    coarses=1,
    poissonsolver_small=PoissonSolver()
)

poissonsolver = ExtraVacuumPoissonSolver(
    gpts=(256, 256, 256),
    poissonsolver_large=poissonsolver0,
    coarses=1,
    poissonsolver_small=PoissonSolver()
)
```

可通过 `poissonsolver.get_description()` 或 `txt` 输出查看最终使用的网格层级与设置。

---

## 2) Delta Self-Consistent Field（ΔSCF / dSCF）

ΔSCF（Delta Self-Consistent Field）是一类用于计算激发态/电离态/电子附着态能量差的方法。传统 ΔSCF 的思路常是“修改占据数（occupation numbers）”，然后重新做一次自洽计算，取总能量差得到激发能。

GPAW 这里介绍的是更一般的 **Linear expansion ΔSCF（线性展开 ΔSCF）**，适合处理**强杂化（strongly hybridized）**的轨道（例如分子化学吸附在过渡金属表面后，分子轨道与金属态强混合）。

---

### 2.1 Linear expansion ΔSCF（线性展开 ΔSCF）

#### 方法定义（保留原公式）

线性展开 ΔSCF 方法（原文标注 “1”）在每次 SCF 迭代中，把指定轨道 $\varphi_a(r)$ 的密度加到总密度中。额外电荷通常从费米能级（Fermi level）取走，以保持体系电中性：

$$n(r) = \sum_nf_{N-1}(T,\varepsilon_n)|\varphi_n(r)|^2 + |\varphi_a(r)|^2$$

其中：

* $N$：总电子数
* $f_{N-1}(T,\varepsilon_n)$：**N−1 电子体系**的 Fermi-Dirac 分布

为了让 band energy 处理正确，需要将目标轨道 $\varphi_a(r)$ 在 Kohn–Sham 轨道基中展开：

$$|\varphi_a\rangle = \sum_nc_{na}|\varphi_n\rangle, \qquad c_{na} = \langle\varphi_n|\varphi_a\rangle$$

目标轨道的 band energy 为：

$$\varepsilon_a = \sum_n|c_{na}|^2\varepsilon_n$$

#### 与传统 ΔSCF 的关系

这是对传统 ΔSCF（只修改 occupation numbers）的推广。如果 $\varphi_a(r)$ 的展开只包含一个归一化项，则会退化到传统 ΔSCF 的情况。

---

### 2.2 简单分子示例：CO 的 $5\sigma \rightarrow 2\pi$ 激发（co.py）

下面示例计算 CO 分子 $5\sigma\rightarrow2\pi$ 的激发能。我们只需要指定 “占据 LUMO（这里是 $2\pi$）”，方法会自动从最高占据轨道（此例为 $5\sigma$）取走一个电子以保持电中性。

关键点：

* `lumo` 是 `AEOrbital` 类的实例：每次迭代计算保存的 $2\pi$ 态在 Kohn–Sham 态中的展开
* 为了得到全电子重叠 $\langle\varphi_n|2\pi\rangle$，除了 pseudowavefunction 外，还需要提供 **projector overlaps（投影子重叠）**

> 代码中的 `[[1.0, lumo, 1]]` 表示：在 `spin=1` 的自旋通道中给 `lumo` 加 1.0 个电子。

```python
from ase.build import molecule
from gpaw import GPAW
from gpaw import dscf

# -------- Ground state calculation（基态）--------
calc = GPAW(mode='fd',
            nbands=8,
            h=0.2,
            xc='PBE',
            spinpol=True,
            convergence={'energy': 100,
                         'density': 100,
                         'eigenstates': 1.0e-9,
                         'bands': -1})

CO = molecule('CO')
CO.center(vacuum=3)
CO.calc = calc

E_gs = CO.get_potential_energy()

# Obtain the pseudowavefunctions and projector overlaps
# of the state which is to be occupied.
# n=5,6 is the 2pix and 2piy orbitals
n = 5
molecule = [0, 1]
wf_u = [kpt.psit_nG[n] for kpt in calc.wfs.kpt_u]
p_uai = [dict([(molecule[a], P_ni[n]) for a, P_ni in kpt.P_ani.items()])
         for kpt in calc.wfs.kpt_u]

# -------- Excited state calculation（激发态）--------
calc_es = GPAW(mode='fd',
               nbands=8,
               h=0.2,
               xc='PBE',
               spinpol=True,
               convergence={'energy': 100,
                            'density': 100,
                            'eigenstates': 1.0e-9,
                            'bands': -1})

CO.calc = calc_es
lumo = dscf.AEOrbital(calc_es, wf_u, p_uai)
# lumo = dscf.MolecularOrbital(calc, weights={0: [0, 0, 0,  1],
#                                             1: [0, 0, 0, -1]})
dscf.dscf_calculation(calc_es, [[1.0, lumo, 1]], CO)

E_es = CO.get_potential_energy()

print('Excitation energy: ', E_es - E_gs)
```

#### 关于 `MolecularOrbital`（代码中注释部分）

被注释的 `lumo = dscf.MolecularOrbital(...)` 提供另一种指定 CO 的 $2\pi$ 轨道的方法，它不需要先做一次分子基态计算。
在这个简单例子里两者结果相同，但对复杂体系通常推荐 `AEOrbital`。

> **重要注意（原文强调）**
> 使用 `AEOrbital` 时，**需要为 dSCF 构造一个新的 calculator 对象**（如示例所示 `calc_es`）。

#### 强制从指定轨道取电子（高级用法）

`dscf.dscf_calculation` 接受一个轨道列表。例如：

* `[[1.0, lumo, 1], [-1.0, pi, 0]]`
  表示：给 `lumo` 加 1 个电子，同时强制从 `pi`（另一 `AEOrbital` 实例）取走 1 个电子（`spin=0`）。

---

### 2.3 激发吸附体（Exciting an adsorbate）：强杂化轨道的典型场景

线性展开 ΔSCF 就是为**强杂化轨道**设计的：例如分子化学吸附（chemisorption）在过渡金属表面。此时传统 ΔSCF 可能失效，因为你想占据的“分子轨道”不再是单一 Kohn–Sham 态能描述的。

原文给出两个脚本（`homo.py` 和 `lumo.py`），用于 CO 吸附在 Pt(111) 顶位（on-top）的 HOMO / LUMO 相关计算。
脚本计算激发态总能量，激发能由基态与激发态能量之差得到。实际工作中通常从已优化结构（例如 `gs.gpw`）重启更常见。

#### `homo.py`（从 HOMO 取电子，构造“空穴”）

要点：

* 先对气相 CO 计算并保存 HOMO 的 pseudo-wavefunction 与 projector overlaps
* `Estart=-100.0, Eend=0.0` 表示只包含费米能以下的态（默认通常是费米能以上）

```python
from ase.visualize import view
from ase.build import fcc111, add_adsorbate
from gpaw import GPAW
from gpaw.mixer import MixerSum
import gpaw.dscf as dscf

filename = 'homo'

c_mol = GPAW(mode='fd',
             nbands=9,
             h=0.2,
             xc='RPBE',
             kpts=(8, 6, 1),
             spinpol=True,
             convergence={'energy': 100,
                          'density': 100,
                          'eigenstates': 1.0e-9,
                          'bands': 'occupied'},
             txt='CO_homo.txt')

calc = GPAW(mode='fd',
            nbands=80,
            h=0.2,
            xc='RPBE',
            kpts=(8, 6, 1),
            eigensolver='cg',
            spinpol=True,
            mixer=MixerSum(nmaxold=5, beta=0.1, weight=100),
            convergence={'energy': 100,
                         'density': 100,
                         'eigenstates': 1.0e-7,
                         'bands': -10},
            txt=filename + '.txt')

# Import Slab with relaxed CO
slab = fcc111('Pt', size=(1, 2, 3), orthogonal=True)
add_adsorbate(slab, 'C', 2.0, 'ontop')
add_adsorbate(slab, 'O', 3.15, 'ontop')
slab.center(axis=2, vacuum=4.0)

view(slab)

molecule = slab.copy()
del molecule[:-2]

# Molecule
molecule.calc = c_mol
molecule.get_potential_energy()

# Homo wavefunction
wf_u = [kpt.psit_nG[4] for kpt in c_mol.wfs.kpt_u]

# Homo projector overlaps
mol = range(len(slab))[-2:]
p_uai = [dict([(mol[a], P_ni[4]) for a, P_ni in kpt.P_ani.items()])
         for kpt in c_mol.wfs.kpt_u]

# Slab with adsorbed molecule
slab.calc = calc
orbital = dscf.AEOrbital(calc, wf_u, p_uai, Estart=-100.0, Eend=0.0)
dscf.dscf_calculation(calc, [[-1.0, orbital, 1]], slab)
slab.get_potential_energy()
```

#### `lumo.py`（占据 LUMO；处理 2π 简并导致的 band 对应问题）

要点：

* $2\pi$ 轨道在 CO 中通常简并（$2\pi_x/2\pi_y$）
* 你希望占据某个特定分量（如 $2\pi_y$），必须在每个 k-point 上判断它对应 band=5 还是 6
  原文做法：与参考 `lumo(k=0)` 做重叠比较，选择更大的那条 band

```python
from numpy import reshape, dot

from ase.visualize import view
from ase.build import fcc111, add_adsorbate
from gpaw import GPAW
from gpaw.mixer import MixerSum
import gpaw.dscf as dscf

filename = 'lumo'

c_mol = GPAW(mode='fd',
             nbands=9,
             h=0.2,
             xc='RPBE',
             kpts=(8, 6, 1),
             spinpol=True,
             convergence={'energy': 100,
                          'density': 100,
                          'eigenstates': 1.0e-9,
                          'bands': -2},
             txt='CO_lumo.txt')

calc = GPAW(mode='fd',
            nbands=80,
            h=0.2,
            xc='RPBE',
            kpts=(8, 6, 1),
            eigensolver='cg',
            spinpol=True,
            mixer=MixerSum(nmaxold=5, beta=0.1, weight=100),
            convergence={'energy': 100,
                         'density': 100,
                         'eigenstates': 1.0e-7,
                         'bands': -10},
            txt=filename + '.txt')

# Import Slab with relaxed CO
slab = fcc111('Pt', size=(1, 2, 3), orthogonal=True)
add_adsorbate(slab, 'C', 2.0, 'ontop')
add_adsorbate(slab, 'O', 3.15, 'ontop')
slab.center(axis=2, vacuum=4.0)

view(slab)

molecule = slab.copy()
del molecule[:-2]

# Molecule
molecule.calc = c_mol
molecule.get_potential_energy()

# Find band corresponding to lumo
lumo = c_mol.get_pseudo_wave_function(band=5, kpt=0, spin=1)
lumo = reshape(lumo, -1)

wf1_k = [c_mol.get_pseudo_wave_function(band=5, kpt=k, spin=1)
         for k in range(c_mol.wfs.kd.nibzkpts)]
wf2_k = [c_mol.get_pseudo_wave_function(band=6, kpt=k, spin=1)
         for k in range(c_mol.wfs.kd.nibzkpts)]

band_k = []
for k in range(c_mol.wfs.kd.nibzkpts):

    wf1 = reshape(wf1_k[k], -1)
    wf2 = reshape(wf2_k[k], -1)
    p1 = abs(dot(wf1, lumo))
    p2 = abs(dot(wf2, lumo))
    if p1 > p2:
        band_k.append(5)
    else:
        band_k.append(6)

# Lumo wavefunction
wf_u = [kpt.psit_nG[band_k[kpt.k]] for kpt in c_mol.wfs.kpt_u]

# Lumo projector overlaps
mol = range(len(slab))[-2:]
p_uai = [dict([(mol[a], P_ni[band_k[kpt.k]]) for a, P_ni in kpt.P_ani.items()])
         for kpt in c_mol.wfs.kpt_u]

#   Slab with adsorbed molecule
slab.calc = calc
orbital = dscf.AEOrbital(calc, wf_u, p_uai)
dscf.dscf_calculation(calc, [[1.0, orbital, 1]], slab)
slab.get_potential_energy()
```

> **实用说明（建议你在中文文档里强调）**
> dSCF/ΔSCF 属于“约束激发态”的自洽计算，收敛往往比基态更难：
>
> * 常需要更稳的 mixer（如 `MixerSum`）
> * 可能需要 `eigensolver='cg'` 或更严格的本征态收敛
> * 需要仔细检查 `nbands`、k 点、以及是否出现占据漂移
>   这些都应当做独立收敛测试，而不是只照搬示例参数。

---

## 3) External potential（外加势 / 外场）

外加势（external potentials）会作用在**所有带电粒子**上，即电子与原子核（electrons and nuclei）。

在脚本中通常是先构造一个外场对象，然后用 `external=` 传给 `GPAW(...)`。

---

### 3.1 恒定电场：`gpaw.external.ConstantElectricField`

示例：沿 z 方向施加 2.5 V/Å 的电场：

```python
from gpaw.external import ConstantElectricField
calc = GPAW(external=ConstantElectricField(2.5, [0, 0, 1]), ...)
```

类签名（原文）：

**class `gpaw.external.ConstantElectricField(strength, direction=[0, 0, 1], tolerance=1e-07)`**
恒定外电场（constant electric field）

* **strength: float**

  * 场强，单位 **V/Å**
* **direction: vector**

  * 极化方向（方向向量）
* **tolerance**

  * 数值容差

> **实用提醒**
> 外电场与周期边界条件的兼容性、以及是否需要配合偶极修正（dipole correction）与真空收敛，依赖体系与求解器设置。
> 遇到 warning 时务必以 `txt` 输出为准，并做 vacuum/kpts 收敛测试。

---

### 3.2 点电荷势：`gpaw.external.PointChargePotential`

示例：两个点电荷（单位为 |e|）：

```python
from gpaw.external import PointChargePotential
pc = PointChargePotential([-1, 1], [[4.0, 4.0, 0.0], [4.0, 4.0, 10.0]])
calc = GPAW(external=pc, ...)
```

类签名（原文）：

**class `gpaw.external.PointChargePotential(charges, positions=None, rc=0.2, rc2=inf, width=1.0)`**

* **charges: list of float**

  * 电荷量，单位 **|e|**
* **positions: (N, 3) array-like**

  * 点电荷位置，单位 **Å**
* **rc: float**

  * 库仑势的内截断（inner cutoff），单位 Å
* **rc2: float**

  * 外截断（outer cutoff），单位 Å
* **width: float**

  * 外截断平滑函数的宽度

数值处理（原文保留）：

* 当 `r < rc` 时，用关于 $r^2$ 的三阶多项式替代 $1/r$，并匹配值、1/2 阶导数与积分。
* 当 `rc2 - width < r < rc2` 时，对 $1/r$ 乘上平滑截断函数（三阶多项式）。

你也可以给 `rc` 负值，此时对所有 r 使用：

```text
(r^4 - rc^4) / (r^5 - |rc|^5)
```

并且 **不在 rc2 处截断**。

---

### 3.3 z 方向台阶势：`gpaw.external.StepPotentialz`

示例：z 方向 step potential：

```python
from gpaw.external import StepPotentialz
step = StepPotentialz(zstep=10, value_right=-7)
calc = GPAW(external=step, ...)
```

类签名（原文）：

**class `gpaw.external.StepPotentialz(zstep, value_left=0, value_right=0)`**

* **zstep: float**

  * 将空间分成左右两侧的 z 值（Å）
* **value_left: float**

  * 左侧（z < zstep）势能值（eV），默认 0
* **value_right: float**

  * 右侧（z >= zstep）势能值（eV），默认 0

---

### 3.4 恒定磁场：`gpaw.bfield.BField`

类签名（原文）：

**class `gpaw.bfield.BField(field)`**
恒定磁场（constant magnetic field）

* **field**

  * 磁场向量，单位为 **eV/bohr-magneton（玻尔磁子）**

---

### 3.5 叠加多个外场：`gpaw.external.PotentialCollection`

可以用 `PotentialCollection` 同时施加多个外加势：

```python
from gpaw.external import ConstantElectric, PointChargePotential
from gpaw.external import PotentialCollection
ext1 = ConstantElectricField(1)
ext2 = PointChargePotential([1, -5], positions=((0, 0, -10), (0, 0, 10)))
collection = PotentialCollection([ext1, ext2])
calc = GPAW(external=collection, ...)
```

类签名（原文）：

**class `gpaw.external.PotentialCollection(potentials)`**

* **potentials: list**

  * 外加势对象列表

> **小提示（兼容性）**
> 上面示例里导入了 `ConstantElectric`，但前文使用的是 `ConstantElectricField`。不同版本可能存在别名或文档笔误。
> 实际使用时请以你安装的 GPAW 版本中可成功 `import` 的类名为准（通常 `ConstantElectricField` 更常见）。

---

## 4) Grids（网格）

本节解释 GPAW 在不同边界条件下，波函数（wave functions）在网格（grid）上的表示方式，以及为什么你设置了 `gpts=(8,8,8)` 但取出的波函数数组可能是 `(7,7,7)`。

---

### 4.1 一个简单设定

假设有一个 `Atoms` 对象，位于边长为 `L` 的立方晶胞：

```python
L = 2.0
atoms = Atoms(cell=(L, L, L), pbc=True)
```

并使用网格间距 `h=0.25 Å`，这等价于 `gpts=(8, 8, 8)`（因为 $L/h = 2.0/0.25 = 8$）。

---

### 4.2 周期边界（`pbc=True`）下的网格索引

当使用周期边界条件时，x 轴（y、z 类似）可以理解为：

```text
 0 1 2 3 4 5 6 7 0 1 2 3 4 5 6 7 0 1 2 3 4 5 6 7 0
-+---------------+---------------+---------------+-
-L               0               L              2*L
```

要点：

* 波函数在 `8×8×8` 网格上表示
* 网格点编号为 `0..7`
* 由于周期性，坐标 `x=L` 与 `x=0` 等价（所以图上会出现重复的 0..7）

---

### 4.3 零边界（`pbc=False`）下的网格索引与数组形状

如果使用零边界条件（Dirichlet，波函数在边界为 0），则 x 轴类似：

```text
  0 1 2 3 4 5 6
+---------------+
0               L
```

要点（原文含义）：

* 在 `x=0 Å` 与 `x=L Å` 处，波函数 **严格为 0**
* 因此只需要存储内部非零点：形成 `7×7×7` 的网格（编号 0..6）

这解释了为什么：

* 你在输入中指定 `gpts=(8,8,8)`
* 但 `get_pseudo_wave_function()` 取出的数组形状可能是 `(7,7,7)`

> **关于文档中的 “padding” 提示（原文：Update this XXX how about padding?）**
> GPAW 在内部数值实现中可能会为了有限差分 stencil、MPI 域分解等分配额外的“padding/ghost”区域；但用户通过接口取出的波函数/密度通常是物理区域（不含内部实现细节）。
> 如果你在做底层开发或需要严格核对网格布局，应查看 `txt` 输出以及 GPAW 的 grid descriptor（如 `calc.wfs.gd` 等内部对象）。

---

### 4.4 示例：pbc=False 时波函数 shape 的变化

原文示例（保留原意；注意部分方法名可能是旧接口风格）：

```python
L = 2.0
atoms = Atoms(...,
              cell=(L, L, L),
              pbc=False)
calc = GPAW(..., gpts=(8, 8, 8))
atoms.SetCalculator(calc)
e = atoms.get_potential_energy()
wf = calc.get_pseudo_wave_function(band=0)
wf.shape
(7, 7, 7)
calc.GetGridSpacings()
array([ 0.25,  0.25])
```

> **实用说明（版本兼容）**
>
> * 新版 ASE 通常推荐：`atoms.calc = calc`（替代 `atoms.SetCalculator(calc)`）
> * 某些 GPAW 方法名在不同版本可能是 `GetGridSpacings()` / `get_grid_spacings()` 这类大小写差异。若报错，请根据 `dir(calc)` 或当前版本 API 调整。

---

## LCAO 模式（LCAO Mode）

GPAW 支持一种替代性的计算模式：**LCAO mode（Linear Combination of Atomic Orbitals，原子轨道线性组合）** [LCAO-article]。在 LCAO 模式中，波函数不再主要用实空间网格（grid-based wave functions）表示，而是用一组“类似原子轨道（atomic orbital-like）”的基函数（basis functions）展开。

* **优点（why use LCAO）**：计算通常显著更便宜（更少自由度、内存占用更低、对大体系更友好）。
* **局限（main limitation）**：精度受限于所选 **basis set（基组）** 的质量与完备性（basis incompleteness）。因此必须做**基组收敛测试（basis convergence）**，尤其是能量差（吸附能、反应能垒、相对稳定性）这类敏感量。

下文将说明：LCAO 的基本原理、如何生成/管理基组、如何运行与优化，以及与材料计算相关的性能与误差（如 BSSE）的处理要点。

> **提示**
> LCAO 模式也可用于 TD-DFT（Time-dependent DFT），通过 `LCAOTDDFT` 模块。

---

### 1. 原理简介（Introduction）

在 LCAO 模式中，Kohn–Sham 赝波函数（Kohn-Sham pseudo wave functions）$\tilde{\psi}*n(\mathbf r)$ 展开在一组“原子型”基函数 $\Phi*\mu(\mathbf r)$ 上，形式类似 SIESTA 方法 [Siesta]：

$$\tilde{\psi}*n(\mathbf r) = \sum*\mu c_{\mu n} \Phi_{\mu}(\mathbf r)$$

这些基函数一般构造成“数值径向函数（numerical radial functions）× 球谐函数（spherical harmonics）”：

$$\Phi_{nlm}(\mathbf{r})
= \Phi_{nlm}(\mathbf r^a + \mathbf R^a)
= \varphi_{nl}(r^a) Y_{lm}(\hat{\mathbf{r}}^a)$$

其中：

* $\mathbf R^a$：原子核 $a$ 的位置
* $\mathbf r^a = \mathbf r - \mathbf R^a$：相对该原子核的坐标

在该近似下，变分自由度不再是实空间波函数本身，而是展开系数 $c_{\mu n}$。因此本征值问题变为广义本征方程：

$$\sum_\nu H_{\mu\nu} c_{\nu n} = \sum_{\nu} S_{\mu\nu} c_{\nu n} \epsilon_n$$

* $H_{\mu\nu}$：基函数表象下的 Hamiltonian 矩阵
* $S_{\mu\nu}$：overlap（重叠）矩阵
* 这通常可通过对（广义）哈密顿矩阵做**直接对角化（direct diagonalization）**来求解（计算量随体系规模通常呈**立方标度（cubic scaling）**）。

> **材料计算中的直观理解**
>
> * FD/PW 模式：自由度 ~ 网格点/平面波数（通常很大，但可系统收敛）。
> * LCAO 模式：自由度 ~ 每个原子使用的基函数数目（通常少很多，因此更快），但需要验证基组是否足够好。

---

### 2. 基组生成（Basis-set generation）

要进行 LCAO 计算，你需要为体系中**每一种元素**准备对应的基组文件（basis-set file）。GPAW 提供 `gpaw basis` 子命令来生成基组。

例如：

```bash
$ gpaw basis --search H.PBE Cl.PBE
```

这会基于已安装的 PAW 数据集（setup）文件 `H.PBE.gz` 与 `Cl.PBE.gz`，生成：

* `H.dzp.basis`
* `Cl.dzp.basis`

其中 `dzp` 表示 **double zeta polarized（双ζ + 极化）**，也是默认基组类型。

**放置位置：**基组文件通常应放在与 GPAW setups 相同的目录中（参见 *Installation of PAW datasets* / `GPAW_SETUP_PATH` 的说明），这样 GPAW 才能自动找到它们。

> **说明与建议**
>
> * 很多较新的 GPAW 安装已经自带常用元素的标准 basis 文件；但如果你想进一步自定义（例如改变 polarization function），或缺失某些元素的 basis 文件，就可以用 `gpaw basis` 重新生成。
> * 基组生成与 PAW setup 的选择密切相关（特别是 `xc`、valence 选取等）。必要时也可参考 `gpaw dataset` 子命令（用于生成 basis 所依赖的 setup 文件）。

查看完整参数列表：

```bash
$ gpaw basis --help
```

---

### 3. 运行 LCAO 计算（Running a calculation）

运行 LCAO 计算时，需要在计算器中同时指定：

* `mode='lcao'`
* `basis=...`（基组名称或映射）

示例：

```python
>>> calc = GPAW(mode='lcao',
>>>             basis='dzp',
>>>             ...)
```

`basis` 的写法与 `setups` 非常类似：可以是字符串，也可以是字典，为不同元素/不同原子指定不同基组。例如：

* `basis={'H': 'dzp', 'O': 'mine', 'C': 'szp'}`

> **大体系性能强烈建议（非常重要）**
> 对较大的 LCAO 体系，**务必启用 ScaLAPACK**，将“立方标度的对角化步骤”并行化，同时分布多矩阵。
> 若条件允许，建议安装并启用 **ELPA** [Elpa]，进一步加速本征求解。相关设置见 `parallel` keyword（并行参数）。

---

### 4. 示例：用 LCAO 优化水分子，并切换到 FD 继续优化（Example）

下面示例用 LCAO 计算器优化（relax）水分子结构。`QuasiNewton` 优化器会使用 LCAO 局域基组计算得到的 forces。

```python
from ase import Atoms
from ase.optimize import QuasiNewton
from gpaw import GPAW

a = 6
b = a / 2

mol = Atoms('H2O',
            [(b, 0.7633 + b, -0.4876 + b),
             (b, -0.7633 + b, -0.4876 + b),
             (b, b, 0.1219 + b)],
            cell=[a, a, a])

calc = GPAW(nbands=4,
            h=0.2,
            mode='lcao',
            basis='dzp')

mol.calc = calc
dyn = QuasiNewton(mol, trajectory='lcao_h2o.traj')
dyn.run(fmax=0.05)
```

优化完成后，你可以**无缝切换到 FD（finite-difference）模式**继续优化，以提升最终精度：

```python
>>> calc.set(mode='fd')
>>> dyn.run(fmax=0.05)
```

> **推荐工作流（材料计算中很常用）**
>
> 1. 用 LCAO（dzp/szp 等）快速预优化/粗收敛（省时间）
> 2. 切换到 FD 或 PW 做精修（forces/energy 更可控地收敛）
> 3. 对关键性质（吸附能、能垒）用更严格设置做收敛测试

---

### 5. 基组的更多说明（More on basis sets）

#### 5.1 zeta 与极化（polarization）

* **最小基组（minimal basis）**：对原子每个价态（valence state）提供一个“原子轨道型函数”。
* 为了提高基组张成空间（span），会给每个价态增加额外径向函数（radial functions）：

  * single-zeta（sz）
  * double-zeta（dz）
  * triple-zeta（tz）…
    zeta 数表示每个价态使用多少个径向函数。

通常还希望加入一个对应于“最低未占据角动量量子数”的基函数，称为 **polarization function（极化函数）**。

一般情况下：

* **dzp（double-zeta polarized）** 往往是获得“合理精度”的最低推荐配置
* GPAW 默认生成的也是 dzp

#### 5.2 以 N 的 dzp 为例：函数数量如何数？

对氮（N）而言，价态为 2s 与 2p：

* 2s：2 个径向函数
* 2p：2 个径向函数
* 极化函数：d 型（一个径向函数）

因此共有 5 个“不同的径向函数”。每个角动量 $l$ 的简并度是 $2l+1$：

* s（l=0）：1 个
* p（l=1）：3 个
* d（l=2）：5 个

所以 GPAW 对一个 N 原子最终会使用总计 **13 个**基函数。

对过渡金属，价态通常包含 s 与 d，则极化函数会是 p 型，因此总基函数数目常见为 **15**（原文给出的典型情况）。

#### 5.3 可视化/检查基组：`gpaw-analyse-basis`

要绘制已生成的基组函数：

```bash
$ gpaw-analyse-basis -f H.dzp.basis O.dzp.basis
```

* 加 `-f`：绘制你指定的文件
* 不加 `-f`：会在 `GPAW_SETUP_PATH` 中寻找第一个匹配文件，而不一定是你指定的那个文件路径

更多选项：

```bash
$ gpaw-analyse-basis --help
```

> **实用建议**
>
> * 做严肃材料计算前，建议至少对关键元素（如吸附原子/过渡金属）检查基函数是否过于“硬/软”，并对 dzp/tzp、是否加更扩展（diffuse）函数做对比收敛测试。
> * LCAO 的误差往往体现在相对能量（能量差）上更敏感，因此“基组收敛”比单点能量更重要。

---

### 6. Ghost atoms 与 BSSE（basis set superposition error）

在表面吸附体系中，表面拥有大量基函数。吸附体（adsorbate）在靠近表面时，可能“借用”表面的基函数自由度，使得总能量虚假降低，从而得到**过强的结合能**。这就是 **BSSE（basis set superposition error，基组叠加误差）**。

一种常用消除方式是在“孤立吸附体”的计算中加入 **ghost atoms**：

* ghost atom **具有基函数**（提供相同的变分自由度）
* 但除此之外不参与计算：**没有 projectors、compensation charges 等**
* 从而保证不同计算中波函数可用的自由度一致（类似 counterpoise 思路）

#### 6.1 如何设置 ghost atoms？

计算方式与普通计算相同：ASE `Atoms` 对象仍包含所有原子（包括 ghost），唯一差别是 ghost 原子的 setup type 设置为 `'ghost'`。

强调：普通原子与对应 ghost 原子的唯一区别是 setup type。

例：Cu 为 ghost，O/H 为普通原子：

```python
>>> GPAW(setups={'Cu' : 'ghost', 'O' : 'paw', 'H' : 'paw'},
         basis='dzp',
         mode='lcao',
         ...)
```

例：除了 Atoms 对象中的第 17 与第 42 个原子用普通 setups，其它全是 ghost：

```python
>>> GPAW(setups={'default': 'ghost', 17: 'paw', 42:'paw'},
         basis='dzp',
         mode='lcao',
         ...)
```

> **实用提醒**
>
> * BSSE 对吸附能、弱相互作用（vdW）、以及分子-分子结合能影响可能非常显著。
> * 如果你用 LCAO 做吸附能，建议把“BSSE 校正（ghost / counterpoise）”作为标准流程之一，至少对代表性构型做对比。

---

### 7. 性能注意事项（Notes on performance）

对于较大的 LCAO 计算：

* **关键：ScaLAPACK**
* **推荐：ELPA**（更快的 eigensolver）
  详见 ScaLAPACK 专章。

#### 7.1 为什么 LCAO 的“其它部分”更重要？

FD（grid-based finite-difference）与 LCAO 的主要区别是赝波函数表示方式；但**密度与势（density/potential）仍然用实空间网格方法处理**。因此在 LCAO 模式下：

* 波函数相关计算更轻，反而密度/势相关计算会占用更高的 CPU 比例
* 所以要注意这些操作的并行与性能（例如 `augment_grids`）

#### 7.2 一个性能配置示例（大体系建议模板）

下面示例展示 LCAO 获得良好性能的一些关键参数。示例体系本身太小，不足以充分展示并行收益，但参数对大体系很有价值：

```python
from ase.build import molecule
from ase.optimize import QuasiNewton
from gpaw import GPAW

atoms = molecule('CH3CH2OH', vacuum=4.0)
atoms.rattle(stdev=0.1)  # displace positions randomly a bit

calc = GPAW(mode='lcao',
            basis='dzp',
            nbands='110%',
            parallel=dict(band=2,  # band parallelization
                          augment_grids=True,  # use all cores for XC/Poisson
                          sl_auto=True,  # enable ScaLAPACK parallelization
                          use_elpa=True))  # enable Elpa eigensolver
atoms.calc = calc

opt = QuasiNewton(atoms, trajectory='opt.traj')
opt.run(fmax=0.05)
```

参数含义（实用解释）：

* `nbands='110%'`：band 数为占据带数的 110%，给出一定空带用于收敛与占据展宽
* `parallel['band']`：band 并行
* `augment_grids=True`：让所有核用于 XC/Poisson 等网格相关工作（对 LCAO 很重要）
* `sl_auto=True`：自动启用 ScaLAPACK 并行
* `use_elpa=True`：启用 ELPA 求解器（若可用）

> **Note（原文提示，了解即可）**
> 原文下一段讨论旧的 `FDPoissonSolver`。该求解器已被 `FastPoissonSolver` 替代，后者一般表现良好，因此那段经验对默认求解器不再适用。
> 但若你出于特殊原因使用旧求解器/自定义 Poisson solver，仍可参考：
>
> * 尽量让每个轴向网格点数可被 8 整除（如 `gpts=(96,96,96)`）
> * Gauss-Seidel 通常比 Jacobi 更快（Poisson 迭代步数可减少约 40%）

---

### 8. 高级基组生成（Advanced basis generation）

`gpaw basis` 命令的后端是类 `gpaw.atom.basis.BasisMaker`。当你需要命令行无法设置的特殊参数（例如修改 valence 选取、增强球半径、特殊的极化策略等），可直接用它生成基组。

**class `gpaw.atom.basis.BasisMaker(valence_data, name=None, run=True, gtxt='-', non_relativistic_guess=False, xc='PBE', save_setup=False)`**
用于创建原子基函数的类。

下面示例：为 Au 生成一个 RPBE 的 double-zeta 基组，并将本来空的 p 态视作 valence，同时使用非标准增强球半径（augmentation sphere size）。

```python
from gpaw.atom.generator import Generator
from gpaw.atom.configurations import parameters
from gpaw.atom.basis import BasisMaker

symbol = 'Au'
args = parameters[symbol]  # Dictionary of default setup parameters
args['rcut'] = 2.6  # Set cutoff of augmentation sphere

generator = Generator(symbol, 'RPBE', gpernode=2 * Generator.default_gpernode)
setup = generator.run(write_xml=False, **args)

bm = BasisMaker.from_setup_and_generator(
    setup, generator, name='special', run=False)

# Create double-zeta RPBE basis where p orbital is considered a valence state
# (ordinary dzp basis would use a smaller p-type Gaussian for polarization)
# The list jvalues indicates which states should be included, and the
# ordering corresponds to the valence states in the setup.
basis = bm.generate(zetacount=2, polarizationcount=0,
                    energysplit=0.1, jvalues=[0, 1, 2],
                    rcutmax=12.0)

basis.write_xml()  # Dump to file 'Au.special.dz.basis'
```

> **实用提醒**
> 基组生成依赖 setup 生成器（setup generator），因此任何适用于 setup 的参数往往也会影响 basis。做这种“特殊基组”时，务必记录参数并做验证（例如与标准 dzp、或与 FD/PW 的对照）。

---

### 9. 其它说明（Miscellaneous remarks）

#### 9.1 FD/PW 模式也会用一次 LCAO 初始化

在 FD 或 PW 模式下，GPAW 会先做一次 LCAO 迭代来初始化波函数与密度。你可以在 FD/PW 模式下指定 `basis` 来提高初猜质量，但**不会影响后续迭代的表示方式**：

```python
>>> calc = GPAW(mode='fd', basis='dzp', ...)
```

#### 9.2 若不指定 basis，默认用什么？

如果在任意模式下你不指定 `basis`，GPAW 会使用 setup 中的 pseudo partial waves $\tilde \phi_i^a(\mathbf r)$，并在半径 **8 Bohr** 处平滑截断，作为基函数。

* 这在多数情况下大致相当于 single-zeta（sz）
* 取决于 setup 中定义的未占据态，有时对某些元素可能接近 single-zeta polarized（szp）

---

### 10. Local Orbitals（局域轨道，LO）

在 LCAO 模式中，可以构造一组更“局域化”的轨道（**Local Orbitals, LOs**），用于定义有效 tight-binding Hamiltonian。

与 Wannier functions（WFs）不同：

* LO 不基于 Kohn–Sham 态投影（projection）
* 不需要为 WFs 提供物理输入（如初猜、投影窗口等）
* LOs 通过对 LCAO Hamiltonian 做**子块对角化（sub-diagonalization）**直接得到

构造方式（要点）：

* 对系统中任意原子，取其原子轨道（AOs）对应的 Hamiltonian 子块进行子对角化
* 得到一组“原子型”的 LOs：天然以原子为中心
* 同一原子上的 LOs 相互正交（orthogonal），但不同原子之间不保证正交
* 还可以只对子系统的一部分原子做 LO（例如只关心量子结中的分子桥、或基底上的吸附体），其余仍保持原 AO 表象

更多细节见 `Local Orbitals` tutorial。

---

### 11. Augmented Basis Set Generation（增强基组：NGTO）

默认 LCAO 基组可以用 **numerical Gaussian-type orbitals（NGTOs）**进行增强。

> Note
> GPAW 的 LCAO 基组是 **PAW basis sets**；而像 Basis Set Exchange 这样的库主要提供全电子（all-electron）基组参数。
> 因此通常只能“借用”其中的 **diffuse orbitals（弥散函数）**参数，因为它们是更平滑的 Gaussians，适合与 PAW/NAO 形式结合。

示例代码（原文保留）：

```python
from gpaw.lcao.generate_ngto_augmented import \
    create_GTO_dictionary as GTO, generate_nao_ngto_basis

gtos_atom = {'H': [GTO('s', 0.02974),
                   GTO('p', 0.14100)],
             'C': [GTO('s', 0.04690),
                   GTO('p', 0.04041),
                   GTO('d', 0.15100)],
             'O': [GTO('s', 0.07896),
                   GTO('p', 0.06856),
                   GTO('d', 0.33200)],
             }

for atom in ['H', 'C', 'O']:
    generate_nao_ngto_basis(atom, xc='PBE', name='aug',
                            nao='dzp', gtos=gtos_atom[atom])
```

这里的 Gaussian 参数对应于 aug-cc-pvdz 基组中的 diffuse functions（示例说明来自 `Basis Set Exchange`）。

> **实用建议**
>
> * 增强基组（尤其加 diffuse）可能改善弱相互作用、负离子态、Rydberg-like 状态等描述，但也可能带来数值线性相关（linear dependence）风险。
> * 使用 NGTO 增强后，建议监控：SCF 收敛、条件数/本征值异常、以及能量差随基组增强的收敛性。

---


## XC Functionals（交换-相关泛函）

本节介绍 GPAW 中交换-相关（exchange-correlation, **XC**）泛函的使用与高级特性，包括：

* **libxc** 泛函调用与已知问题（特别是 **MGGA**）
* **Exact exchange（EXX）/Hybrid functionals（杂化泛函）**：PW 与 FD 两套实现、以及非自洽（post-processing）计算
* **Range separated functionals（RSF，范围分离杂化泛函）**：原理、实现限制与实用技巧（含 γ/omega 调参）
* **vdW-DF / BEEF-vdW**：非局域范德华泛函、自洽/非自洽、以及 `libvdwxc` 加速后端
* **van der Waals correction**：TS09 等色散修正（dispersion correction）
* **QNA**（Quasi-non-local approximation）功能示例

> **实用总览（建议先看这一段）**
>
> * **常规材料计算**：优先用 semilocal（LDA/GGA）如 `PBE`，必要时用 vdW（如 `BEEF-vdW`、`vdW-DF`）或 D3/D4/TS09 等修正方案（视体系而定）。
> * **Hybrid/EXX**：周期体系通常首选 **PW 后端**（支持 k 点、对称性、力），分子/大真空盒子且只需 Γ 点时可考虑 FD 后端（但功能与速度限制明显）。
> * **RSF（LCY-/CAMY-）**：在 GPAW 的实现继承了 FD-EXX 的限制（Γ 点、非周期、RMMDIIS 等），对真空要求很高，适合“分子/团簇/大盒子”特定任务。
> * **MGGA**：务必注意 libxc 的 **FHC（Fermi hole curvature）**设置与某些 MGGA 需要 Laplacian 的问题，否则可能报错或严重不收敛。

---

### Exchange and correlation functionals（交换与相关泛函）

#### Libxc

GPAW 提供对 `libxc` 中各类 XC 泛函的访问接口。通常你可以直接用 `xc='...'` 指定 libxc 关键字，或者用 `xc={'name': ..., ...}` 传入带参数的泛函设置（例如带 `stencil` 等）。

---

#### Known Problems（已知问题与排查）

##### 1) MGGAs：某些 MGGA 需要 Laplacian（目前不提供）

部分 MGGA（例如使用 Becke-Roussel 89 交换的：`MGGA_X_BR89+MGGA_C_TPSS`）需要密度的 **laplacian（拉普拉斯项）**。在本文撰写时 GPAW 尚未提供该量，因此使用这类泛函会直接抛出异常（exception）。

> **实用建议**：遇到 MGGA 报错时，优先检查 libxc 泛函是否需要 Laplacian；必要时换用不需要 Laplacian 的 MGGA 或退回 GGA。

##### 2) MGGAs：libxc 默认启用 FHC 导致赝势代码错误与收敛问题

libxc 默认强制 **Fermi hole curvature（FHC）**，这会在使用赝势/PAW 的代码中引入错误结果与收敛问题。

* 在 **libxc > 7.0**：该行为可以（并且将会）在运行时关闭
* 在 **libxc < 7.0**：必须在编译时关闭（安装 libxc 时加 `--disable-fhc`）

你可以用下面的脚本检查当前 libxc 是否以 `--disable-fhc` 编译（该脚本的目的是：尝试跑一个 TPSS MGGA，并检查 SCF 是否能收敛）：

```python
"""check is libxc is compiled with --disable-fhc (needed for mggas)"""
from ase.build import molecule
from gpaw import GPAW, KohnShamConvergenceError
from gpaw.utilities.adjust_cell import adjust_cell

vacuum = 4.0
h = 0.3


def test_mgga_lxc_fhc():
    cluster = molecule('CO')
    adjust_cell(cluster, border=vacuum, h=h)
    calc = GPAW(xc='MGGA_X_TPSS+MGGA_C_TPSS',
                mode='fd',
                h=h,
                maxiter=14,
                convergence={
                    'energy': 0.5,
                    'density': 1.0e-1,
                    'eigenstates': 4.0e-1})
    cluster.calc = calc
    try:
        cluster.get_potential_energy()
    except KohnShamConvergenceError:
        pass
    assert calc.scf.converged


if __name__ == '__main__':
    test_mgga_lxc_fhc()
```

> **补充说明（实用角度）**
>
> * 这类检查脚本更像“环境一致性测试”，适合在集群/新环境中部署 GPAW 时先跑一遍。
> * MGGA 本身往往更敏感：即便环境正确，也建议提高计算稳定性（更稳的 mixer、合理展宽、更严格/更合适的 SCF 参数）。

---

### Exact exchange（精确交换）与 Hybrid functionals（杂化泛函）

目前 GPAW 有两套 exact exchange（EXX）实现：

1. **Plane-wave（PW）模式实现**：支持 k 点、利用对称性、可计算 forces（力）。
2. **Finite-difference（FD）模式实现**：见 `hybrid.py`。可自洽处理 **仅 Γ 点**体系（分子与大盒子）；**不提供 forces**（并且功能/速度限制较多）。

> **实用建议**
>
> * 周期材料（需要 k 点、应力、力优化等）优先选 **PW 后端**。
> * 分子/团簇、只算 Γ 点、主要关心能量与轨道性质时才考虑 FD-EXX（并接受其限制）。

---

#### Self-consistent plane-wave implementation（PW 自洽实现）

PW 自洽实现可用的 `xc` 名称包括：`EXX`, `PBE0`, `HSE03`, `HSE06`, `B3LYP`。

一个 PW 模式的 `HSE06` 示例：

```python
atoms.calc = GPAW(
    mode={'name': 'pw', 'ecut': ...},
    xc='HSE06',
    ...)
```

或者用字典获得更细控制（full control）：

```python
atoms.calc = GPAW(
    mode={'name': 'pw', 'ecut': ...},
    xc={'name': 'HYB_GGA_XC_HSE06',
        'fraction': 0.25,
        'omega': 0.11,
        'backend': 'pw'},
    ...)
```

> **Note（重要）**
> 目前 `EXX`, `PBE0`, `B3LYP` 需要写成 `xc={'name': ..., 'backend': 'pw'}` 的形式！

---

#### Non self-consistent plane-wave implementation（PW 非自洽实现 / 事后修正）

这类方法基于一个已经完成的自洽 DFT 计算（通常是 semilocal，如 PBE/LDA），再**非自洽**地评估杂化泛函能量或本征值修正。常用于：

* 快速得到 hybrid 能量（避免完整自洽 hybrid 的高成本）
* 计算 hybrid 对能带/能级的修正（例如在 PBE 的基础上做 PBE0 修正）

##### 1) 非自洽 hybrid 能量：`non_self_consistent_energy`

**`gpaw.hybrids.energy.non_self_consistent_energy(calc, xcname, ftol=1e-09)`**
用于计算杂化泛函的非自洽能量贡献。对占据数小于 `ftol` 的态，EXX 积分会跳过（节省计算）。

示例：

```python
energies = non_self_consistent_energy('<gpw-file>',
                                      xcname='HSE06')
e_hyb = energies.sum()
```

对自洽能量的修正量是：`energies[1:].sum()`。

返回的能量分解（单位 eV）依次为：

1. DFT 总自由能（free energy，**未外推到零展宽**）
2. −（DFT XC 能量）
3. Hybrid 的半局域 XC 能量（semi-local part）
4. EXX core-core 能量
5. EXX core-valence 能量
6. EXX valence-valence 能量

> **实用提醒**
>
> * 如果你比较能量差（如吸附能），要确保各体系使用相同的 `ftol`、展宽与相同的基态参考（同样的 `gpw` 生成策略），否则误差可能不相消。

##### 2) 非自洽 hybrid 本征值：`non_self_consistent_eigenvalues`

**`gpaw.hybrids.eigenvalues.non_self_consistent_eigenvalues(calc, xcname, n1=0, n2=0, kpt_indices=None, snapshot=None, ftol=1e-09)`**
计算杂化泛函的非自洽本征值贡献。

* 基于一个自洽 DFT 计算（calc）
* 只计算 `n1` 到 `n2-1` 的本征值（默认 all bands）
* `kpt_indices` 默认 all k-points（IBZ）
* 对占据数 < `ftol` 的态跳过 EXX 积分
* `snapshot='name.json'`：可输出每个 k 点完成时的快照，便于长计算监控/恢复

返回三个数组（形状 `(nspins, nkpts, n2-n1)`，单位 eV）：

```python
nsceigs = non_self_consistent_eigenvalues
eig_dft, vxc_dft, vxc_hyb = nsceigs('<gpw-file>', xcname='PBE0')
eig_hyb = eig_dft - vxc_dft + vxc_hyb
```

可参考教程：`PBE0 calculations for bulk silicon`。

##### 3) 示例：氢原子 EXX（自洽 vs 在 LDA 上的非自洽）

```python
"""EXX hydrogen atom.

Compare self-consistent EXX calculation with non self-consistent
EXX calculation on top of LDA.
"""
from ase import Atoms
from ase.units import Ry
from gpaw import GPAW, PW
from gpaw.hybrids.eigenvalues import non_self_consistent_eigenvalues
from gpaw.hybrids.energy import non_self_consistent_energy

atoms = Atoms('H', magmoms=[1.0])
atoms.center(vacuum=5.0)

# Self-consistent calculation:
atoms.calc = GPAW(mode=PW(600),
                  xc={'name': 'EXX', 'backend': 'pw'})
eexx = atoms.get_potential_energy() + atoms.calc.get_reference_energy()

# Check energy
eexxref = -1.0 * Ry
assert abs(eexx - eexxref) < 0.001

# ... and eigenvalues
eig1, eig2 = (atoms.calc.get_eigenvalues(spin=spin)[0] for spin in [0, 1])
eigref1 = -1.0 * Ry
eigref2 = ...  # ?
assert abs(eig1 - eigref1) < 0.03
# assert abs(eig2 - eigref2) < 0.03

# LDA:
atoms.calc = GPAW(mode=PW(600),
                  xc='LDA')
atoms.get_potential_energy()

# Check non self-consistent eigenvalues
result = non_self_consistent_eigenvalues(atoms.calc,
                                         'EXX',
                                         snapshot='h-hse-snapshot.json')
eiglda, vlda, vexx = result
eig1b, eig2b = (eiglda - vlda + vexx)[:, 0, 0]
assert abs(eig1b - eig1) < 0.04
assert abs(eig2b - eig2) < 1.1

# ... and energy
energies = non_self_consistent_energy(atoms.calc, 'EXX')
eexxb = energies.sum() + atoms.calc.get_reference_energy()
assert abs(eexxb - eexx) < 0.03
```

---

#### Self-consistent finite-difference implementation（FD 自洽实现 / `hybrid.py`）

当前 FD-EXX（hybrid.py）实现存在以下功能缺失或限制：

* **周期体系支持不足**：代码可能不报错，但结果未充分测试
* **不支持 k-point sampling**：不考虑多 k 点，甚至不支持复波函数，因此肯定不可用
* **Forces**：含 Fock operator 的力实现“有但未测试”
* **速度极慢**：瓶颈在 Fock operator 的 Poisson 积分求解；目前每个 SCF 周期都用零初猜迭代求解，仍需进一步优化

加速思路之一：用**粗网格（coarse grid，用于 wave functions）**替代**细网格（fine grid，用于 densities）**计算 Fock potentials，理论上可加速约 8 倍。这可以通过在 `xc` 关键字中设置（示例 `test_coarse.py`，原文提及）。

并行方面：基于 domain decomposition 的并行是完全支持的。

FD-EXX 可用于：

* `PBE0`（杂化）
* `EXX`（Hartree-Fock 型精确交换）

支持的泛函（原文表）：

| Functional | Type   | Reference |
| ---------- | ------ | --------- |
| EXX        | Global | [AB98]    |
| PBE0       | Global | [Ba94]    |
| B3LYP      | Global | [SZ12]    |
| CAMY-B3LYP | RSF-LR | [AT08]    |
| CAMY-BLYP  | RSF-LR | [SZ12]    |
| CAMY-B3LYP | RSF-LR | [SZ12]    |
| LCY-BLYP   | RSF-LR | [SZ12]    |
| LCY-PBE    | RSF-LR | [SZ12]    |

这里：

* **Global**：全局杂化（global hybrid），在空间中各处都用固定比例的 Hartree-Fock exchange（EXX）
* **RSF-SR / RSF-LR**：范围分离（range-separated）

  * SR（short-range）：EXX 比例随距离增大而减小
  * LR（long-range）：EXX 比例随距离增大而增大
    更多信息见 `Range separated functionals (RSF)`。

##### 关于未占据态与光学激发（非常关键的“坑”）

该实现**缺少 optimized effective potential（OEP）**，因此使用 EXX 的未占据态往往更像（激发）电子亲和能（electron affinity）的近似。结果是：

* 用含 Hartree-Fock exchange 的轨道直接做 **lrTDDFT** 光学激发通常不是好选择

作为补救，GPAW 实现了 **Improved Virtual Orbitals（IVOs）** [HA71]。可通过关键字选择激发基（excitation basis）：

* `excitation=...`
* `excited=...`：从 HOMO 向下计数的态索引

示例：用 RSF + IVO 计算 NaCl 的激发能，并进一步做 lrTDDFT：

```python
"""Calculate the excitation energy of NaCl by an RSF using IVOs."""
from ase.build import molecule
from ase.units import Hartree
from gpaw import GPAW, setup_paths
from gpaw.mpi import world
from gpaw.occupations import FermiDirac
from gpaw.test import gen
from gpaw.eigensolvers import RMMDIIS
from gpaw.utilities.adjust_cell import adjust_cell
from gpaw.lrtddft import LrTDDFT

h = 0.3  # Gridspacing
e_singlet = 4.3
e_singlet_lr = 4.3

if setup_paths[0] != '.':
    setup_paths.insert(0, '.')

gen('Na', xcname='PBE', scalarrel=True, exx=True, yukawa_gamma=0.40)
gen('Cl', xcname='PBE', scalarrel=True, exx=True, yukawa_gamma=0.40)

c = {'energy': 0.005, 'eigenstates': 1e-2, 'density': 1e-2}
mol = molecule('NaCl')
adjust_cell(mol, 5.0, h=h)
calc = GPAW(mode='fd', txt='NaCl.txt',
            xc='LCY-PBE:omega=0.40:excitation=singlet',
            eigensolver=RMMDIIS(), h=h, occupations=FermiDirac(width=0.0),
            spinpol=False, convergence=c)
mol.calc = calc
mol.get_potential_energy()
(eps_homo, eps_lumo) = calc.get_homo_lumo()
e_ex = eps_lumo - eps_homo
assert abs(e_singlet - e_ex) < 0.15
calc.write('NaCl.gpw')

lr = LrTDDFT(calc, txt='LCY_TDDFT_NaCl.log',
             restrict={'istart': 6, 'jend': 7}, nspins=2)
lr.write('LCY_TDDFT_NaCl.ex.gz')
if world.rank == 0:
    lr2 = LrTDDFT.read('LCY_TDDFT_NaCl.ex.gz')
    lr2.diagonalize()
    ex_lr = lr2[1].get_energy() * Hartree
    assert abs(e_singlet_lr - e_singlet) < 0.05
```

IVO 在 lrTDDFT 中的支持参考 Berman 和 Kaldor [BK79] 的工作。

##### RMM-DIIS 限制与“先 PBE 再 EXX”的工作流

如果计算所需 `nbands` 超过数据集（dataset）提供的 band 数，GPAW 会随机初始化缺失的 bands。并且在该 Hartree-Fock exchange 路线下只能使用 **RMM-DIIS** eigensolver，因此：

* 态可能不会收敛到能量最低的态（会“卡在”不理想的子空间）

一种常见规避方法：先用 semilocal（如 PBE）收敛得到波函数，然后在其基础上切换到含 EXX 的设置（例如 PBE0 + IVO）。示例：

```python
"""Test calculation for unoccupied states using IVOs."""

from ase.build import molecule
from gpaw.utilities.adjust_cell import adjust_cell
from gpaw import GPAW, KohnShamConvergenceError, FermiDirac
from gpaw.eigensolvers import CG, RMMDIIS

calc_parms = [
    {'xc': 'PBE0:unocc=True',
     'eigensolver': RMMDIIS(niter=5),
     'convergence': {
         'energy': 0.005,
         'bands': -2,
         'eigenstates': 1e-4,
         'density': 1e-3}},
    {'xc': 'PBE0:excitation=singlet',
     'convergence': {
         'energy': 0.005,
         'bands': 'occupied',
         'eigenstates': 1e-4,
         'density': 1e-3}}]


def calc_me(atoms, nbands):
    """Do the calculation."""
    molecule_name = atoms.get_chemical_formula()
    atoms.set_initial_magnetic_moments([-1.0, 1.0])
    fname = '.'.join([molecule_name, 'PBE-SIN'])
    calc = GPAW(mode='fd',
                h=0.25,
                xc='PBE',
                eigensolver=CG(niter=5),
                nbands=nbands,
                txt=fname + '.log',
                occupations=FermiDirac(0.0, fixmagmom=True),
                convergence={
                    'energy': 0.005,
                    'bands': nbands,
                    'eigenstates': 1e-4,
                    'density': 1e-3})
    atoms.calc = calc
    try:
        atoms.get_potential_energy()
    except KohnShamConvergenceError:
        pass
    if calc.scf.converged:
        for calcp in calc_parms:
            calc.set(**calcp)
            try:
                calc.calculate(system_changes=[])
            except KohnShamConvergenceError:
                break
        if calc.scf.converged:
            calc.write(fname + '.gpw', mode='all')


loa = molecule('NaCl')
adjust_cell(loa, border=6.0, h=0.25, idiv=16)
loa.center()
loa.translate([0.001, 0.002, 0.003])
nbands = 25
calc_me(loa, nbands)
```

---

### Range separated functionals (RSF)（范围分离杂化泛函）

#### Introduction（原理简介）

范围分离泛函（RSF）是杂化泛函的一类。与 PBE0/B3LYP 这类**固定比例**（global）混合 EXX 与 DFT exchange 不同，RSF 按两点间距离 $r_{12}$ 混合 EXX 与 DFT exchange：

* 例如 PBE0：1/4 EXX + 3/4 DFT exchange（空间中处处相同）
* RSF：混合比例随 $r_{12}$ 变化，通过软函数 $\omega_\mathrm{RSF}(\gamma, r_{12})$

为实现该点，Hartree-Fock 交换积分中的库仑核：

$$\frac{1}{r_{12}} = \frac{1}{|r_1 - r_2|}$$

被拆成短程（SR）与长程（LR）两部分：

$$\frac{1}{r_{12}} = \underbrace{\frac{1 - [\alpha + \beta ( 1 - \omega_\mathrm{RSF} (\gamma, r_{12}))]}{r_{12}}}*{\text{SR, DFT}} + \underbrace{\frac{\alpha + \beta ( 1 - \omega*\mathrm{RSF} (\gamma, r_{12}))}{r_{12}}}_{\text{LR, HFT}}$$

* SR 部分由 LDA/GGA（如 PBE）的 exchange 处理
* LR 部分由 Hartree-Fock exchange（HFT/EXX）处理
* $\alpha,\beta$ 是混合参数：

  * $\alpha\neq 0,\beta=0$ 类似传统 global hybrid
  * $\alpha=0,\beta\neq 0$ 常标记为 `LC` + semilocal 名称（如 LC-PBE）
  * $\alpha\neq 0,\beta\neq 0$ 常标记为 `CAM` + 名称（如 CAM-BLYP）

分离函数 $\omega_\mathrm{RSF}$ 常见两类：

* 互补误差函数：$\omega_\mathrm{RSF}=\mathrm{erfc}(\gamma r_{12})$
* Slater 函数：$\omega_\mathrm{RSF}=e^{(-\gamma r_{12})}$

在 Gaussian 基组代码中 erfc 计算上更方便；而 Slater 函数在 Rydberg 态、charge transfer excitations 上往往更好。命名上：

* 用 Slater 函数的 RSF 会在标记中加 “Y”：如 **LCY-PBE**, **CAMY-B3LYP**
* 用 erfc 的保持 **LC-PBE**, **CAM-B3LYP**

另外，分离函数还依赖 screening factor（屏蔽因子）$\gamma$。$\gamma$ 的最佳取值仍在讨论中，部分观点认为它应与密度有关。多数 RSF 给了默认 $\gamma$，但也允许针对电离能、CT 激发、双自由基阳离子键长曲线等任务**调参（tune γ）**。

---

#### Implementation（在 GPAW 中的实现）

GPAW 的 RSF 实现由两部分组成：

* semilocal 部分：通过 libxc 实现
* Hartree-Fock exchange：通过 `hybrid.py` 实现

GPAW 使用 Slater 分离函数：$\omega_\mathrm{RSF}=e^{(-\gamma r_{12})}$。除可自定义 $\gamma$ 外，实现了以下泛函：

| Functional | α    | β    | $a_0^{-1}$ | Reference |
| ---------- | ---- | ---- | ---------- | --------- |
| CAMY-BLYP  | 0.2  | 0.8  | 0.44       | [AT08]    |
| CAMY-B3LYP | 0.19 | 0.46 | 0.34       | [SZ12]    |
| LCY-BLYP   | 0.0  | 1.0  | 0.75       | [SZ12]    |
| LCY-PBE    | 0.0  | 1.0  | 0.75       | [SZ12]    |

由于 RSF 基于 FD-EXX（hybrid.py）实现，继承其优缺点，简述如下：

* 支持自洽 RSF 计算
* **只能 Γ 点**
* **只能非周期边界条件**
* **eigensolver 只能用 RMMDIIS**

> **Important（非常重要：真空盒子大小）**
> RSF 的一个主要优势是保留 exchange potential 的 $1/r$ 渐近行为，因此对中性或带负电体系必须用**很大的盒子**。
> 原文建议：每个原子周围至少 **6 Å** 真空；对阴离子（anions）需要更大。

更多实现与 RSF 细节参考 [WW18]，以及 [Wu16] 的详细讨论（原文引用）。

---

#### Simple usage（最简用法）

选择 RSF 泛函通常很简单，例如：

```python
"""First example for using RSF."""
from ase import Atoms
from gpaw import GPAW
from gpaw.eigensolvers import RMMDIIS
from gpaw.occupations import FermiDirac

h = 0.30
co = Atoms('CO', positions=[(0, 0, 0), (0, 0, 1.15)])
co.center(5)
# c = {'energy': 0.005, 'eigenstates': 1e-4}  # Usable values
c = {'energy': 0.1, 'eigenstates': 3, 'density': 3}  # Values for test
calc = GPAW(mode='fd', txt='CO.txt', xc='LCY-PBE', convergence=c,
            eigensolver=RMMDIIS(), h=h,
            occupations=FermiDirac(width=0.0), spinpol=False)
co.calc = calc
co.get_potential_energy()
```

从这个例子可以直接看到三个关键点：

1. `xc='LCY-PBE'` 即可选择 RSF
2. 必须用 `eigensolver=RMMDIIS()`
3. 收敛标准通常要“放松一点”（降低收敛要求），否则很容易卡 SCF

---

#### Improving results（提高结果可靠性与数值稳定性）

RSF 计算中有两个常见“额外复杂性”：

1. 自洽时需要 core 电子贡献，这必须在生成 PAW datasets 时就计算并写入（即 datasets 必须带 EXX/RSF 所需信息）

2. 在笛卡尔网格上算 exchange 需要数值求解（screened）Poisson 方程。对带电体系（例如态与自身交换）需要用一个 Gaussian 表示“过量电荷”并先中和，再求解 Poisson，再把 Gaussian 的解加回。若要移除的电荷是 off-center，则中和电荷中心应匹配过量电荷中心以避免人为偶极。
   GPAW 提供了使用电荷中心的 Poisson solver 选项：`PoissonSolver(use_charge_center=True)`。

示例（带专用 datasets 与 Poisson solver）：

```python
"""Calculations using RSF with dedicated datasets and Poisson-solver."""
from ase import Atoms
from gpaw import GPAW, setup_paths
from gpaw.poisson import PoissonSolver
from gpaw.eigensolvers import RMMDIIS
from gpaw.occupations import FermiDirac
from gpaw.test import gen

if setup_paths[0] != '.':
    setup_paths.insert(0, '.')

for atom in ['C', 'O']:
    gen(atom, xcname='PBE', scalarrel=True, exx=True,
        yukawa_gamma=0.75)

h = 0.30
co = Atoms('CO', positions=[(0, 0, 0), (0, 0, 1.15)])
co.center(5)

# c = {'energy': 0.005, 'eigenstates': 1e-4}  # Usable values
c = {'energy': 0.1, 'eigenstates': 3, 'density': 3}  # Values for test

calc = GPAW(mode='fd', txt='CO.txt', xc='LCY-PBE', convergence=c,
            eigensolver=RMMDIIS(), h=h,
            poissonsolver=PoissonSolver(use_charge_center=True),
            occupations=FermiDirac(width=0.0), spinpol=False)
co.calc = calc
co.get_potential_energy()
```

PAW datasets 也可用命令生成（原文）：

* `gpaw-setup -f PBE -x --gamma=0.75 C O`

---

#### Tuning γ（调参：优化 omega/gamma）

原文指出，γ 的最佳值仍在讨论中。一个常见调参策略是：通过调 γ，使得 **HOMO 的负本征值**匹配计算得到的电离能（IP）。在 GPAW 中可以通过 `omega` 传入 γ：

```python
"""Calculation utilizing RSF with optimized gamma."""
from ase import Atoms
from gpaw import GPAW, setup_paths
from gpaw.poisson import PoissonSolver
from gpaw.eigensolvers import RMMDIIS
from gpaw.occupations import FermiDirac
from gpaw.test import gen

# IP for CO using LCY-PBE with gamma=0.81 after
# dx.doi.org/10.1021/acs.jctc.8b00238
IP = 14.31

if setup_paths[0] != '.':
    setup_paths.insert(0, '.')

for atom in ['C', 'O']:
    gen(atom, xcname='PBE', scalarrel=True, exx=True,
        yukawa_gamma=0.81)

h = 0.30
co = Atoms('CO', positions=[(0, 0, 0), (0, 0, 1.15)])
co.center(vacuum=5)

# c = {'energy': 0.005, 'eigenstates': 1e-4}  # Usable values
c = {'energy': 0.1, 'eigenstates': 3, 'density': 3}  # Values for test

calc = GPAW(mode='fd', txt='CO.txt', xc='LCY-PBE:omega=0.81', convergence=c,
            eigensolver=RMMDIIS(), h=h,
            poissonsolver=PoissonSolver(use_charge_center=True),
            occupations=FermiDirac(width=0.0), spinpol=False)
co.calc = calc
co.get_potential_energy()
(eps_homo, eps_lumo) = calc.get_homo_lumo()
assert abs(eps_homo - -IP) < 0.35
```

---

### vdW-DF and BEEF-vdW（非局域 vdW 泛函与 BEEF-vdW）

GPAW 支持 vdW-DF 家族泛函，方式包括：

* GPAW 内置接口（built-in interface）
* 外部 `libvdwxc` 库（更快更可扩展）

> **注意**：两者使用不同的 kernel 参数化，因此结果会有轻微差异。

GPAW 自洽实现了一系列 vdW-DF 类型泛函，以及 `BEEF-vdW`。常见 vdW-DF 变体包括：

* `vdW-DF`, `vdW-DF2`, `vdW-DF-cx`, `optPBE-vdW`, `optB88-vdW`, `C09-vdW`

其中：

* `vdW-DF-cx` 仅在启用 `libvdwxc` 时可用
* vdW-DF 泛函的自旋极化推广（spin-polarized generalization）也仅在 `libvdwxc` 下可用

实现细节（原文要点）：

* 自洽实现使用 Perez-Soler FFT 算法来计算 Rutgers-Chalmers 非局域相关（nonlocal correlation）的能量与势（该相关项在实空间原本是 6 维积分）
* 也提供一种非自洽方法：直接在实空间对积分求和（较慢）

---

#### Doing a vdW-DF calculation（如何进行 vdW-DF 计算）

* 强烈推荐用 **自洽 FFT 方法**，优于实空间方法
* 许多情况下 vdW-DF 的密度与普通 GGA（如 PBE）密度非常相似，因此常见策略是：
  **先做 GGA 自洽 → 用 FFT 方法非自洽评估 vdW-DF 总能差**
* 但 **vdW-DF forces** 需要自洽势，因此做结构优化时仍应自洽跑 vdW-DF

> **大体系建议**：vdW 泛函可能很贵，建议使用 `libvdwxc`，通常能把 vdW 评估效率提高一个数量级，并可扩展到任意大体系并行。

---

#### Selfconsistent vdW-DF calculations（自洽 vdW-DF）

```python
from ase import *
from gpaw import GPAW
vdw = 'vdW-DF'
atoms = ...
calc = GPAW(xc=vdw, ...)
atoms.calc = calc
e = atoms.get_potential_energy()
```

---

#### Perturbative vdW-DF calculations（非自洽/微扰式评估）

```python
from gpaw import GPAW
xc = 'vdW-DF'
calc = GPAW('input.gpw')
GGA_energy = calc.get_potential_energy()
vdWDF_diff = calc.get_xc_difference(xc)
vdWDF_energy = GGA_energy + vdWDF_diff
```

你可以把 `vdW-DF` 替换成 `vdW-DF2`, `vdW-DF-cx`（需 libvdwxc）, `optPBE-vdW`, `optB88-vdW`, `C09-vdW` 等。

要显式使用 `libvdwxc` 后端，例如：

* `xc={'name': 'vdW-DF', 'backend': 'libvdwxc'}`

> 注意：`libvdwxc` 使用不同 kernel 参数化，计算值会略有差异。

---

#### Non-default FFT parameters（非默认 FFT 插值参数）

vdW-DF 非局域核的 spline 插值由一组参数决定。若要改变这些参数，需要显式初始化 vdW-DF 基类。示例：在 `vdW-DF2` 中把 `Nalpha` 改为 24：

```python
from ase import *
from gpaw import GPAW
from gpaw.xc.vdw import VDWFunctional
vdw = VDWFunctional('vdW-DF2', Nalpha=24)
atoms = ...
calc = GPAW(xc=vdw, ...)
atoms.calc = calc
e = atoms.get_potential_energy()
```

---

#### Real-space method vdW-DF（实空间方法：很慢，仅建议很小体系）

实空间非自洽方法可能只对极小体系有意义。使用方式：

```python
from gpaw.xc.vdw import VDWFunctional
vdw = VDWFunctional('vdW-DF', fft=False, nspins=1, ncut=0.0005)
```

* `nspins=1`：非自旋极化（spin-paired）
* `nspins=2`：自旋极化
* `ncut`：密度 cutoff，小于该值的密度不纳入 6D 积分

---

#### BEEF-vdW functional（BEEF-vdW 与误差估计）

`BEEF-vdW` 使用 `vdW-DF2` 的非局域相关能与势，自洽实现于 GPAW。它的一个重要特点是：可对用 BEEF-vdW 自洽得到的“目标量”（通常是**相对能量**）给出一个**误差估计（ensemble error estimate）**：

* 在 BEEF-vdW 的电子密度上，非自洽地应用一组 XC 泛函“集合”（ensemble）
* 通过这些预测的方差来估计不确定性（注意：主要针对相对能量，而非绝对总能）

示例：计算 H₂ 的结合能与其 ensemble 误差估计：

```python
from ase import *
from gpaw import GPAW
from ase.dft.bee import BEEFEnsemble
xc = 'BEEF-vdW'
h2 = Atoms('H2',[[0.,0.,0.],[0.,0.,0.75]])
h2.center(vacuum=3)
cell = h2.get_cell()
calc = GPAW(mode='fd', xc=xc)
h2.calc = calc
e_h2 = h2.get_potential_energy()
ens = BEEFEnsemble(calc)
de_h2 = ens.get_ensemble_energies()
del h2, calc, ens
h = Atoms('H')
h.set_cell(cell)
h.center()
calc = GPAW(mode='fd', xc=xc)
h.calc = calc
e_h = h.get_potential_energy()
ens = BEEFEnsemble(calc)
de_h = ens.get_ensemble_energies()
E_bind = 2*e_h - e_h2
dE_bind = 2*de_h[:] - de_h2[:]
dE_bind = dE_bind.std()
```

补充要点（原文）：

* `BEEFEnsemble` 模块最近从 GPAW 移到了 ASE
* 默认 ensemble XC 泛函数目为 **2000**，以保证误差估计收敛

  * `ens.get_ensemble_energies(N)` 可改为 N
* `de_h2` 和 `de_h` 是 2000 个能量扰动值数组
* **非常重要**：不要对“总能扰动数组本身”分别取标准差；应该对“相对能扰动”（如结合能扰动）取标准差（示例中正确做法是先算 `dE_bind` 再 `.std()`）

---

### libvdwxc（更快的 vdW-DF 家族实现）

`libvdwxc` 提供 vdW-DF 家族非局域 vdW 泛函的快速、可扩展实现。要使用它：

1. 安装 `libvdwxc`
2. 用它重新编译 GPAW

依赖：

* FFTW3 与 FFTW3-MPI
* 对超大体系可安装 PFFT 以更好扩展性
  对“现实大小体系”，FFTW3-MPI 通常已经高效，甚至可能比 PFFT 更快。

运行时通过 backend 指定，例如：

```python
from ase.build import molecule
from gpaw import GPAW

atoms = molecule('H2O')
atoms.center(vacuum=3.0)
# There are these functionals: vdW-DF, vdW-DF2, vdW-DF-cx, optPBE-vdW,
# optB88-vdW, C09-vdW, BEEF-vdW, and mBEEF-vdW.
# There are three modes: serial, mpi, and pfft. Default is auto.
calc = GPAW(mode='fd',
            xc={'name': 'BEEF-vdW', 'backend': 'libvdwxc', 'mode': 'mpi'})
atoms.calc = calc
atoms.get_potential_energy()
```

并行提示（原文）：

* `libvdwxc` 会在 domain decomposition 上自动用尽可用核
* 若你还并行了 k 点或 bands，特别是 PW 模式，务必传 `parallel={'augment_grids': True}`（或等价写法）让所有核也参与 XC/Poisson/vdW 的网格工作（见 `Parallel runs`）

> **注意**：`libvdwxc 0.4` 没有 stress term（应力项）实现。

---

### van der Waals correction（色散修正：TS09 等）

Tkatchenko 与 Scheffler 提出了在 PBE 基础上的修正（通常称 TS09）。该方法几乎不增加额外计算成本，但表现很好（原文给出 S26 测试集对 CCSD 参考的误差比较，单位 meV）：

| ·                       | PBE | TPSS | vdW-DF | vdW-DF2 | TS09 | Grimme D4 |
| ----------------------- | --- | ---- | ------ | ------- | ---- | --------- |
| Mean absolute deviation | 115 | 154  | 76     | 48      | 16   | 14        |
| RMS deviation           | 108 | 128  | 60     | 42      | 21   | 14        |

说明（原文）：

* 误差来自 S26 测试集相对 CCSD 能量
* GPAW 计算使用 `h=0.18` 且至少 `4 Å` 真空
* TS09 与 FHI-aims 的结果一致性很好
* Grimme D4 在 `github` 可用（原文）

#### 计算 S26 测试集示例（原文代码）

```python
import sys

from ase import Atoms
from ase.build.connected import connected_indices
from ase.calculators.vdwcorrection import vdWTkatchenko09prl
from ase.data.s22 import data
from ase.parallel import paropen

from gpaw import GPAW, FermiDirac
from gpaw.analyse.hirshfeld import HirshfeldPartitioning
from gpaw.analyse.vdwradii import vdWradii
from gpaw.utilities.adjust_cell import adjust_cell

try:
    from dftd4 import D4_model
except ModuleNotFoundError:
    pass

h = 0.18
box = 4.

xc = 'TS09'
if len(sys.argv) > 1:
    xc = sys.argv[1]

f = paropen('energies_' + xc + '.dat', 'w')
print('# h=', h, file=f)
print('# box=', box, file=f)
print('# molecule E[1]  E[2]  E[1+2]  E[1]+E[2]-E[1+2]', file=f)
for molecule in data:
    print(molecule, end=' ', file=f)
    ss = Atoms(data[molecule]['symbols'],
               data[molecule]['positions'])
    # split the structures
    s1 = ss[connected_indices(ss, 0)]
    s2 = ss[connected_indices(ss, -1)]
    assert len(ss) == len(s1) + len(s2)
    if xc == 'TS09' or xc == 'TPSS' or xc == 'M06-L' or xc == 'dftd4':
        c = GPAW(mode='fd', xc='PBE', h=h, nbands=-6,
                 occupations=FermiDirac(width=0.1))
    else:
        c = GPAW(mode='fd', xc=xc, h=h, nbands=-6,
                 occupations=FermiDirac(width=0.1))
    E = []
    for s in [s1, s2, ss]:
        s.calc = c
        adjust_cell(s, box, h=h)
        if xc == 'TS09':
            s.get_potential_energy()
            cc = vdWTkatchenko09prl(HirshfeldPartitioning(c),
                                    vdWradii(s.get_chemical_symbols(), 'PBE'))
            s.calc = cc
        elif xc == 'dftd4':
            s.get_potential_energy()
            cc = D4_model(xc='PBE', calc=c)
            s.calc = cc
        if xc == 'TPSS' or xc == 'M06-L':
            ene = s.get_potential_energy()
            ene += c.get_xc_difference(xc)
            E.append(ene)
        else:
            E.append(s.get_potential_energy())
    print(E[0], E[1], E[2], end=' ', file=f)
    print(E[0] + E[1] - E[2], file=f)
    f.flush()
f.close()
```

> **实用提示**
>
> * 这里 TS09 是通过 ASE 的 `vdWcorrection` 叠加实现的，依赖 Hirshfeld 分区与 vdW 半径。
> * 对大型体系（尤其周期材料），要评估其适用性与收敛性（真空、网格、展宽等），并确保不同体系的设置一致以便误差相消。

---

### Quasi-non-local exchange correlation approximation（QNA）

原文在 “Rationale” 处未给出具体说明（仅保留标题）。这里给出使用方式与接口说明。

#### Using QNA（使用 QNA）

示例（原文）：

```python
from gpaw import GPAW, PW
from ase.lattice.compounds import L1_2

QNA = {'alpha': 2.0,
       'name': 'QNA',
       'orbital_dependent': False,
       'parameters': {'Au': (0.125, 0.1), 'Cu': (0.0795, 0.005)},
       'setup_name': 'PBE',
       'type': 'qna-gga'}

atoms = L1_2(['Au','Cu'],latticeconstant=3.74)
calc = GPAW(mode=PW(300),
            xc = QNA,
            kpts=kpts,
            txt='AuCu3_QNA.txt')
atoms.get_potential_energy()
```

* `xc` 通过字典形式传入
* `parameters` 对不同元素给出参数
* `setup_name` 指明使用哪个 setup 作为基础（这里是 `PBE`）
* `type='qna-gga'` 指明 QNA 的类型

类接口（原文）：

**class `gpaw.xc.qna.QNA(atoms, parameters, qna_setup_name='PBE', alpha=2.0, override_atoms=None, stencil=2)`**

* **`get_description()`**

  * 返回该泛函的长描述字符串，或 `None`。
* **`todict()`**

  * 返回 XC 泛函的字典表示。
  * 该表示适用于 libxc kernels；其他类可能应覆盖该函数且不应依赖这一默认实现。

> **实用建议（扩展说明）**
> QNA 属于较专门的泛函形式，通常需要对目标体系/元素参数有明确来源与验证。若用于材料性质预测（晶格常数、弹性等），建议：
>
> * 明确记录参数来源
> * 与常用 GGA（如 PBE）及可能的实验/高精度参考做对照
> * 做 k 点、cutoff、应力/几何等收敛测试，避免把数值误差误认为“泛函改进”

---

## Custom convergence criteria（自定义收敛判据）

GPAW 的自洽场循环（SCF, self-consistent field）收敛判据通过 `convergence` 字典控制。除了默认判据（energy/density/eigenstates 等）之外，你还可以：

* 添加更多判据（例如 forces、work function、minimum iterations）
* 改变某些判据的检查方式（例如检查最近 3 次还是 4 次迭代）
* 以类（Criterion）方式编写你自己的收敛判据

> **实用建议**
>
> * 绝大多数材料计算用默认收敛判据即可。
> * 只有在遇到“SCF 看似收敛但物理量仍漂移”（如功函数、力不稳定）或需要特定策略（如固定迭代次数、强制检查 forces）时，才建议自定义。

---

### 1) 额外的 convergence 关键字（Additional convergence keywords）

`convergence` 字典除了默认字典中的键外，还支持额外键，例如：

* `'forces'`：力的收敛
* `'work function'`：功函数（work function）收敛
* `'minimum iterations'`：最少迭代次数
* （完整列表与参数见 `Built-in criteria`）

例如：要求功函数在最近三次 SCF 迭代中的变化不超过 0.001 eV：

```python
from gpaw import GPAW

convergence = {'work function': 0.001}

calc = GPAW(...,
            convergence=convergence)
```

注意：上面只显式写了 `'work function'`，但 **默认判据（energy、eigenstates、density）仍然会以默认阈值启用**。

如果你想“等效关闭”某个默认判据，可以把它设为 `np.inf`（无穷大），例如：

```python
import numpy as np
convergence = {'energy': np.inf}  # 基本等价于不检查能量收敛
```

> **实用提醒**
> 把默认判据设为 `np.inf` 会让 SCF 退出条件更依赖你额外设置的判据（或 maxiter）。请确保你的判据足够合理，否则容易“假收敛”或“永不退出”。

---

### 2) 修改判据行为（Changing criteria behavior）

除了用浮点数作为快捷写法，还可以用判据类来精细控制其行为。

例如默认 `convergence={'energy': 0.0005}` 表示检查**最近 3 次**（n_old=3）的能量变化是否小于 0.0005（单位见 `Energy` 判据定义）。如果你希望检查最近 4 次迭代，则可以这样写：

```python
from gpaw.convergence_criteria import Energy

convergence = {'energy': Energy(tol=0.0005, n_old=4)}
```

并且：

* `convergence={'energy': 0.0005}` 实际是 `convergence={'energy': Energy(0.0005)}` 的快捷写法
  即字典值 `0.0005` 会变成 `Energy` 的第一个位置参数。

---

### 3) 收敛 forces（Converging forces）

你可以直接要求 forces 收敛，例如：

```python
convergence = {'forces': 0.01}
```

其含义是：与前一次 SCF 迭代相比，每个原子的力变化（以向量差的模长衡量）中，**最大的那个变化**必须小于 `0.01 eV/Å`。

由于计算 forces 会额外消耗时间和内存，forces 判据默认 `calc_last=True`（先满足其他判据再开始检查 forces），以提高效率。

如果你希望每一步 SCF 都检查 forces（例如你在调试 forces 的稳定性），可以这样做：

```python
from gpaw.convergence_criteria import Forces

convergence = {'forces': Forces(0.01, calc_last=False)}
```

#### 3.1 相对 forces 收敛（按系统最大力比例）

你还可以把 forces 的收敛阈值设成“相对于当前系统最大力”的比例（rtol）。这在几何优化早期特别有用：如果体系离极小点很远，力很大，此时没必要把 SCF/forces 收敛得非常严格。

例如：把 forces 收敛到“最大力的 10%”：

```python
import numpy as np
from gpaw.convergence_criteria import Forces

# Converge forces to 10% of the highest force.
convergence = {'forces': Forces(atol=np.inf, rtol=0.1)}
```

如果同时提供 `atol` 和 `rtol`，则每个 SCF 循环中会取“更严格”的那个标准：

```python
from gpaw.convergence_criteria import Forces

# 在几何优化过程中：当力较大时，用 rtol=0.1；
# 当力降低后，atol=0.01 eV/Å 会变成更严格的标准。
convergence = {'forces': Forces(atol=0.01, rtol=0.1)}
```

> **实用建议（几何优化常见搭配）**
>
> * 早期优化：forces 判据可适当放松（rtol 方式很实用）。
> * 接近收敛结构：再切换到更严格的 `atol`（例如 0.01 eV/Å 或更小），保证最终 forces 的一致性与可重复性。

---

### 4) 示例：固定迭代次数（Example: fixed iterations）

你可以让 SCF “固定运行 N 次迭代”，常用于调试、生成中间态、或某些非标准工作流。

做法是：

1. 把默认判据都设为 `np.inf`（相当于关闭）
2. 用 `'minimum iterations'` 指定最小迭代次数
3. 同时确保 `maxiter` 大于该值，否则 SCF 会先被 maxiter 截停

例如：运行**恰好 10 次**迭代（实践中通常还需要配合 `maxiter=10` 才能“恰好”退出）：

```python
import numpy as np

convergence = {'energy': np.inf,
               'eigenstates': np.inf,
               'density': np.inf,
               'minimum iterations': 10}
```

`MinIter` 也可以与其他判据联合使用：它只负责保证“至少迭代 N 次”，即使其他判据提前满足，也不会提前退出。

---

### 5) 编写你自己的判据（Writing your own criteria）

如果你希望根据某个自定义物理量/数值量来判断 SCF 收敛，可以继承 `Criterion` 并实现 `__call__` 等方法。结构示例：

```python
from gpaw.convergence_criteria import Criterion


class MyCriterion(Criterion):
    name = 'my criterion'  # 必须是唯一名字（unique name）
    tablename = 'mycri'    # <= 5 个字符，会作为 SCF 表头显示
    calc_last = False      # True: 等其他判据都满足后才检查（适合昂贵判据）

    def __init__(self, ...):
        ...
        # self.description 会打印到输出日志顶部
        self.description = 'My custom criterion with tolerance ...'

    def __call__(self, context):
        ...
        # context 包含当前计算状态的引用（如 Hamiltonian、波函数等）
        converged = ...  # True/False：是否满足
        entry = ...      # <= 5 字符：写入 SCF 表格的一列
        return converged, entry

    def reset(self):
        ...  # 每当 SCF 重启时清理内部状态


calc = GPAW(...,
            convergence={'custom': [MyCriterion(0.01, 4)]}
           )
```

要点：

* 所有用户自写判据必须通过 `convergence` 字典中的特殊键 `'custom'` 注入
* `'custom'` 对应一个列表，你可以放多个判据实例

> **Note（重要）**
>
> 1. 如果你把计算写成重启文件（`calc.write('out.gpw')`），GPAW 重新打开 `out.gpw` 时**不会自动恢复你的自定义判据类**。你需要在重启后手动重新添加。
> 2. 如果你同时运行多个 GPAW calculator 实例，请确保每个 calculator 都拿到**独立的自定义判据实例**（built-in 判据内部会复制，不需要你操心；自定义的需要）。

---

### 6) 内置判据列表（Built-in criteria）

下表列出 GPAW 内置的判据类、在 `convergence` 字典中使用的快捷名（name attribute）、是否默认启用、是否 `calc_last`、以及是否会覆盖其它判据。

| class        | name attribute（字典键）  | default? | calc_last? | override_others? |
| ------------ | -------------------- | -------- | ---------- | ---------------- |
| Energy       | `energy`             | Yes      | No         | No               |
| Density      | `density`            | Yes      | No         | No               |
| Eigenstates  | `eigenstates`        | Yes      | No         | No               |
| Forces       | `forces`             | No       | Yes        | No               |
| WorkFunction | `work function`      | No       | No         | No               |
| MinIter      | `minimum iterations` | No       | No         | No               |
| MaxIter      | `maximum iterations` | No       | No         | Yes              |

下面给出各内置判据的类说明（按原文翻译并补充必要解释）。

---

#### `Energy(tol, *, n_old=3, relative=True)`：总能量收敛

**class `gpaw.convergence_criteria.Energy(tol, *, n_old=3, relative=True)`**
总能量收敛判据。

参数：

* **tol: float**

  * 收敛阈值：最近 `n_old` 次（外推的）总能量的最大波动（peak-to-peak difference）不得超过 tol。
* **n_old: int**

  * 比较的能量值数量。例如 n_old=3 表示比较“当前 + 前两次”的能量波动。
* **relative: bool**

  * 若 True：使用“每个价电子归一化”的能量（eV/(valence electron)）；
  * 若 False：使用总能量（eV）。

> **实用说明**
>
> * `relative=True` 对不同规模体系更稳健（能量阈值随电子数缩放）。
> * 对能量差极敏感任务（能垒、相对稳定性），建议做收敛测试而非只调阈值。

---

#### `Density(tol)`：电子密度收敛

**class `gpaw.convergence_criteria.Density(tol)`**
电子密度收敛判据。

参数：

* **tol: float**

  * 收敛阈值：密度变化的积分绝对值（integrated absolute value），并按价电子数归一化。单位：
    `electrons/(valence electron)`。

---

#### `Eigenstates(tol)`：本征态残差收敛

**class `gpaw.convergence_criteria.Eigenstates(tol)`**
本征态收敛判据。

参数：

* **tol: float**

  * 收敛阈值：Kohn–Sham 方程残差平方的积分，并按价电子数归一化。单位：
    `eV^2/(valence electron)`。

> **实用说明**
> 本判据在不同模式（例如 LCAO）下可能不完全起相同作用；若你发现 eigenstates 判据对收敛表现不敏感，需要结合能量/密度/forces 做综合判断。

---

#### `Forces(atol, rtol=inf, calc_last=True)`：力收敛

**class `gpaw.convergence_criteria.Forces(atol, rtol=inf, calc_last=True)`**
forces 收敛判据。

参数：

* **atol: float**

  * 绝对收敛阈值。每个原子与前一次迭代相比的力变化，用 l2-norm（欧氏距离）衡量；所有原子中最大的变化必须小于 atol。单位：`eV/Å`。
* **rtol: float**

  * 相对阈值。把“最大力变化”与 `rtol * (系统当前最大力)` 比较，满足则收敛。
* **calc_last: bool**

  * True：等其他判据先收敛后再检查 forces（更省）。
  * False：每一步 SCF 都检查 forces（更严格但更耗资源）。

---

#### `WorkFunction(tol=0.005, n_old=3)`：功函数收敛

**class `gpaw.convergence_criteria.WorkFunction(tol=0.005, n_old=3)`**
功函数收敛判据。

参数：

* **tol: float**

  * 最近 `n_old` 次功函数的最大波动不得超过 tol。单位：eV。
* **n_old: int**

  * 比较的功函数值数量。

---

#### `MinIter(n)`：最少迭代次数

**class `gpaw.convergence_criteria.MinIter(n)`**
强制 SCF 至少运行 n 次迭代。

参数：

* **n: int**

  * SCF 在退出前必须完成的最少迭代次数。

---

#### `MaxIter(n)`：最多迭代次数（可覆盖其它判据）

**class `gpaw.convergence_criteria.MaxIter(n)`**
强制 SCF 最多运行 n 次迭代，达到则退出（并可能覆盖其它判据）。

参数：

* **n: int**

  * SCF 在退出前允许的最大迭代次数。

---

---

## Pipek-Mezey Wannier Functions（Pipek–Mezey Wannier 函数，PMWF）

### 1) Introduction（简介）

**Pipek-Mezey Wannier functions（PMWF）** 是一种与 **Maximally Localized Wannier Functions（MLWF，最大局域化 Wannier 函数；常见为 Foster–Boys 局域化）**不同的局域化轨道构造方法。

PMWF 的特点：

* 具有更强的“化学直觉”（chemical intuition）
* 能较好地保持 **σ（sigma）** 与 **π（pi）** 类型轨道的区分
* 以 spread function（展宽函数）衡量时，PMWF 通常可达到与 MLWF **相当的局域性**
* 但 MLWF 在实践中经常会把化学上不同的轨道混合（mix chemically distinct orbitals），而 PMWF 对此更“克制”

> **实用理解**
> 如果你的目标是做化学键分析、局域轨道分析、或者希望得到“σ/π 更清晰”的局域轨道表示，PMWF 往往比 MLWF 更符合直觉。

---

### 2) Localization（适用范围与约束）

PMWF 适用于：

* 三种计算模式：**LCAO、PW、FD**
* 两类边界条件：**open boundary（非周期/开边界）**与 **periodic boundary（周期边界）**

但对周期体系有一个明确要求：

* 必须使用 **uniform Monkhorst-Pack grid（均匀 Monkhorst–Pack k 点网格 saying）**

> **实用建议**
>
> * 周期体系中，优先使用 `kpts=(n1, n2, n3)` 这类均匀 MP 网格设置。
> * 如果你之前使用的是 `kpts={'density': ...}` 等方式，请确认它生成的是均匀 MP 网格，并且与 PMWF 的要求兼容（必要时显式写 size）。

---

---

## The Perdew-Zunger Self-Interaction Correction（PZ-SIC，自相互作用校正）

Perdew–Zunger Self-Interaction Correction（**PZ-SIC**）是一种用于修正 DFT 中 **self-interaction error（自相互作用误差）**的方案。其实现属于高级功能，通常计算更昂贵，也更敏感。

原文要点：

* 该实现支持三种模式：**PW、FD、LCAO**
* 使用“direct minimization（直接极小化）”方法（见原文所指 `here`）
* 由于该泛函 **不是 unitary invariant functional（不具幺正不变性）**：
  为了找到最低能量态，需要使用 **complex orbitals（复数轨道）**

> **为什么需要复数轨道？（直观解释）**
> 对非幺正不变的能量泛函，轨道的“旋转/混合”会改变能量；允许复数轨道可以扩大可变分的轨道空间，有助于避免陷入较高能的局部极小点，从而更接近真实最低能态。

---

### 示例 1：FD 模式下的 PZ-SIC（使用 PM_PZ 局域化）

```python
import numpy as np
from ase import Atoms
from gpaw import FD, GPAW
from gpaw.directmin.etdm_fdpw import FDPWETDM

# Water molecule:
d = 0.9575
t = np.pi / 180 * 104.51
H2O = Atoms('OH2',
            positions=[(0, 0, 0),
                       (d, 0, 0),
                       (d * np.cos(t), d * np.sin(t), 0)])
H2O.center(vacuum=5.0)

calc = GPAW(mode=FD(force_complex_dtype=True),
            xc='PBE',
            occupations={'name': 'fixed-uniform'},
            eigensolver=FDPWETDM(localizationtype='PM_PZ',
                                 functional={'name': 'PZ-SIC',
                                             'scaling_factor':
                                                 (0.5, 0.5)},
                                 grad_tol_pz_localization=1.0e-4),
            mixer={'backend': 'no-mixing'},
            symmetry='off'
            )

H2O.set_calculator(calc)
H2O.get_potential_energy()
H2O.get_forces()
```

**参数/关键词的实用解释（建议写进中文参考文档）**：

* `mode=FD(force_complex_dtype=True)`：强制使用复数数据类型（complex dtype），以允许 complex orbitals。
* `occupations={'name': 'fixed-uniform'}`：固定占据（fixed occupations）。这常用于直接极小化方法，避免占据数在过程中漂移。
* `eigensolver=FDPWETDM(...)`：使用 ETDM（Exponential Transformation Direct Minimization）类直接极小化求解器（FD/PW 版本）。

  * `localizationtype='PM_PZ'`：局域化类型（与 PM/PZ 相关）。
  * `functional={'name': 'PZ-SIC', 'scaling_factor': (0.5, 0.5)}`：选择 PZ-SIC 泛函，并设置缩放因子（两元组；具体含义以实现为准，常用于缩放校正强度）。
  * `grad_tol_pz_localization=1.0e-4`：局域化梯度阈值（收敛控制）。
* `mixer={'backend': 'no-mixing'}`：关闭密度混合（no-mixing）。直接极小化路线通常不需要传统 SCF mixer。
* `symmetry='off'`：关闭对称性。对复数轨道与非幺正不变泛函，禁用对称性约束通常更稳妥。

> **ASE 接口小提示**
> 示例使用 `H2O.set_calculator(calc)`。在新 ASE 风格中也可写 `H2O.calc = calc`。两者效果相同，按你的 ASE 版本习惯即可。

---

### 示例 2：PW 模式下的 PZ-SIC（原文说明）

要使用 PW 模式：

* 导入 PW 模式（如 `from gpaw import PW`）
* 将示例中的 `FD` 替换成 `PW`

原文说明：**“To use PW mode, just import PW mode and replace FD with PW.”**

> **实用提醒**
> PW 下如果涉及 k 点/周期体系，建议先确认 PZ-SIC 的算法路径与设置是否适合你的体系规模与并行条件（该类高级泛函往往比常规 GGA 昂贵得多）。

---

### 示例 3：LCAO 模式下的 PZ-SIC（direct minimization）

```python
import numpy as np
from ase import Atoms
from gpaw import GPAW, LCAO
from gpaw.directmin.etdm_lcao import LCAOETDM

# Water molecule:
d = 0.9575
t = np.pi / 180 * 104.51
H2O = Atoms('OH2',
            positions=[(0, 0, 0),
                       (d, 0, 0),
                       (d * np.cos(t), d * np.sin(t), 0)])
H2O.center(vacuum=5.0)

calc = GPAW(mode=LCAO(force_complex_dtype=True),
            xc='PBE',
            occupations={'name': 'fixed-uniform'},
            eigensolver=LCAOETDM(localizationtype='PM_PZ',
                                 functional={'name': 'PZ-SIC',
                                             'scaling_factor':
                                                 (0.5, 0.5)}),
            mixer={'backend': 'no-mixing'},
            nbands='nao',
            symmetry='off'
            )

H2O.calc = calc
H2O.get_potential_energy()
H2O.get_forces()
```

这里相对 FD 示例的关键差异：

* `mode=LCAO(force_complex_dtype=True)`：LCAO 也启用复数轨道
* `eigensolver=LCAOETDM(...)`：LCAO 版本的 ETDM
* `nbands='nao'`：band 数等于原子轨道数（NAO 数），这是 LCAO 下常见设置，保证基组空间完整被纳入。

---

### PZ-SIC 使用时的实用建议（补充）

1. **先用普通泛函（如 PBE）得到合理结构/初猜**
   PZ-SIC 往往更难收敛、更敏感。建议先用 PBE（FD/PW/LCAO 均可）优化结构并保存 `.gpw`，再在其基础上做 PZ-SIC。

2. **对真空、网格/ecut、基组做更严格收敛测试**
   自相互作用校正改变能量与轨道的敏感性，数值误差更容易被放大。

3. **注意 symmetry 与 occupations**
   原文示例明确 `symmetry='off'` 与 `fixed-uniform`，通常是为了减少额外约束与不稳定因素。若你想开启对称性或使用其它占据方式，建议逐步验证。

---

# 理论（Theory）

本页翻译并整理 GPAW 文档中与 **SCF 收敛（density mixing）** 和 **恒电位电化学界面（Solvated Jellium Method, SJM）** 相关的理论内容。为便于实际使用，保留关键英文术语（如 *Pulay mixing*、*work function*、*jellium slab*、*implicit solvent* 等），并补充了一些常见使用建议与注意事项。

---

## 密度混合（Density Mixing）

密度混合（density mixing）是 DFT 自洽场（SCF, self-consistent field）迭代中最核心的数值加速/稳定手段之一。其目标是把“新得到的输出密度”与“历史输入密度”按一定规则组合，减少振荡、加速收敛。GPAW 中常用的是 **Pulay mixing**（历史密度的线性组合 + 近似最优更新）。

### 在 GPAW 中指定混合方案（Specifying a Mixing Scheme in GPAW）

在 GPAW 计算器中通过 `mixer` 关键字指定混合器（mixer）与其参数：

```python
from gpaw import GPAW, Mixer
calc = GPAW(..., mixer=Mixer(beta=0.05, nmaxold=5, weight=50.0), ...)
```

上面这组参数常被作为“当默认混合不收敛时”的推荐备选（更保守、更稳定）。

---

### Mixer 类型与物理含义（Mixer classes and what they mix）

GPAW 提供多种 mixer 类，主要差别在于“混合对象”是总密度、各自旋密度、还是密度/磁化密度的组合。Pulay mixing 可以基于：

1. **分别混合两个自旋通道的密度**：`Mixer`

   * **注意**：对于自旋极化（spin-polarized）体系，除非磁矩被固定（fixed magnetic moment），否则通常不适用/不稳定（原文提示：This will not work … unless the magnetic moment is fixed）。
2. **混合总电子密度（total density）**：`MixerSum2`
3. **密度矩阵按自旋分别混合，同时赝电子密度按自旋求和后混合**：`MixerSum`
4. **分别混合总密度与磁化密度（magnetization density）**：`MixerDif`

   * 这里磁化密度定义为两自旋密度之差：
     [
     m(\mathbf r)=n_\uparrow(\mathbf r)-n_\downarrow(\mathbf r)
     ]

> **经验建议（很实用）**
>
> * 自旋极化分子/团簇：`MixerDif` 往往是比较稳的默认选择。
> * 体相/周期体系：`MixerSum` 有时更好（尤其在某些 bulk 金属或磁性体系里）。
> * 若体系强振荡或“来回跳”：优先减小 `beta`（更保守），或减少历史步数 `nmaxold`（对某些原子/过渡金属体系可能有效）。

---

### Mixer 参数解释（beta / nmaxold / weight）

所有 mixer 类通常接受参数：

* `beta`：线性混合系数（linear mixing coefficient）

  * 越大收敛越快但越容易振荡；越小越稳定但可能更慢。
* `nmaxold`：参与 Pulay 历史混合的“旧密度”数量（history length）
* `weight`：度量（metric）相关的权重参数，用于控制不同波长密度误差在混合中的相对重要性（原文：long wavelength changes weighted higher）

`MixerDif` 额外还接受一组针对磁化密度的参数：

* `beta_m`、`nmaxold_m`、`weight_m`

> **快速调参套路（故障排查版）**
>
> 1. SCF 强烈振荡/发散：先把 `beta` 降到 0.02–0.05 区间；必要时 `nmaxold` 降到 1–3。
> 2. SCF 很慢但稳定：可以逐步把 `beta` 增大一点（例如从 0.02 → 0.05 → 0.1）。
> 3. 自旋体系不稳：优先尝试 `MixerDif`（并确认初始磁矩合理，必要时固定总磁矩）。

---

## Solvated Jellium（constant-potential electrochemistry，溶剂化果冻模型/恒电位电化学）

### 概览（Overview）

**Solvated Jellium Method（SJM）** 是一种在 DFT 中模拟电化学界面（electrochemical interfaces）的简化方法，完整理论描述见 [Kastlunger2018]。SJM 的核心能力是：通过调节体系电子数来控制模拟电极的电势（electrode potential），通常表现为界面上方的 **work function（功函数）** 或相关的电势参照量。

SJM 支持两种运行方式：

* **constant-charge（恒电荷）**：固定额外电子数（excess electrons）
* **constant-potential（恒电位）**：自动迭代调整电子数，使电势达到目标值

SJM 计算器（`SJM` calculator）使用方式与普通 `GPAW` 计算器类似：能返回能量与力（energy and forces），但可在“固定电势”条件下做到这一点。

> **关于能量的 Legendre transform（勒让德变换）**
> 恒电位条件下，“直接输出的能量”在热力学意义上通常对应于对电荷自由度做了适当变换后的势能形式（原文提示：Please see the note below on the Legendre-transform of the energy）。实际工作中，比较不同构型/反应路径时应确认你比较的是同一热力学势（比如固定电势下的势能差）。

**计算成本**：电势控制靠简单迭代技术实现。在做轨迹（trajectory）计算时（如结构弛豫或 NEB）：

* 第一张图像（first image）通常比常规 DFT 更慢（需要把电势“平衡”到目标值）
* 后续图像成本增加较小
* 文档给出的经验估计：整条轨迹的额外成本通常 **< 50%**

实践指南见：`Solvated Jellium Method (SJM)` tutorial。

---

### 理论背景（Theoretical background）

SJM 的设计哲学是：用尽量简单的模型捕捉关键物理，同时避免引入“伪效应（spurious effects）”。SJM 由两部分组成：

1. **jellium（果冻背景电荷）**
2. **implicit solvent（隐式溶剂）**

此外，在 SJM 中常常还会使用 **explicit solvent（显式溶剂，如水分子）**。原因是：隐式溶剂的主要目的通常不是精细描述每个物种的溶剂化结构，而是用于**屏蔽净电场（screen the net field）**；而表面反应物种的真实溶剂化更常由显式水分子承担。

下面分别说明。

---

### 果冻板（jellium slab）：充电与电荷补偿（The jellium slab: charging）

在周期体系（periodic system）中不能有净电荷（net charge）。因此当你向体系加入额外的分数电子（fractional electrons）时，必须有等量反号的补偿电荷（counter charge）。

在 GPAW 中可以用 `JelliumSlab` 方便地实现这一点：

* 在模拟胞（cell）的某个固定区域加入“涂抹开的背景电荷”（smeared-out background charge）
* 同时把体系电子数增加/减少，数量等于 slab 的总电荷

直观上，你希望看到：

* 额外电子主要**局域在金属表面一侧（electrode surface side）**
* 而不是跑到模拟体相的一侧（bulk side）

SJM 通过两点实现这种“单侧充电”：

1. **只把 jellium 区域放在模拟胞的一侧**
2. 使用 **dipole correction（偶极修正）**（SJM 默认包含），在静电上把模拟胞两侧“解耦”（electrostatically decouple）

> 原文中提到“4.4 V 与 4.3 V 两个模拟的差异”，这通常对应于图示：目标电势略变时，体系通过改变少量电子数使表面电荷分布/电势剖面发生细微变化。

#### jellium 的符号与区域形状

jellium 区域常被理解为：

* 一片“涂抹开的正背景电荷” + 相应数量的额外电子

但符号也可以反过来：

* 变成“涂抹开的负背景” + 同时减少总电子数

因此同一个工具既可向正方向也可向负方向扰动电子数，从而把电势调到目标值。

**重要特性**：jellium 区域不与原子重叠（does not overlap any atoms）。这与“在整个晶胞施加均匀背景电荷”的做法不同，后者可能与原子/分子产生不期望的伪相互作用（spurious interactions），从而扭曲电子结构。SJM 通过把补偿电荷限制在不含原子的区域来避免这一问题。

另外，jellium 区域不必一定是规则 slab 几何；它也可以跟随隐式溶剂的 cavity 形状（例如用 `jelliumregion` 关键字，详见 SJM 文档）。

---

### 溶剂化：屏蔽强电场（The solvation: screening）

如果只有“多余电子 + jellium 补偿电荷”，它们会在反应区形成一个人为偏大的电场（artificially high potential field）。为了屏蔽这一强电场，SJM 在显式溶剂层上方引入 **implicit solvent（隐式溶剂）**，并让其完全包围 jellium 补偿电荷所在区域。

SJM 使用 Held 与 Walter 的隐式溶剂模型 [Held2014]：本质上是改变真空区域的介电常数（dielectric constant），从而实现电场屏蔽。更多内容可参考：`Continuum Solvent Model (CSM)` tutorial。

#### 为什么既用显式溶剂又用隐式溶剂？

* 显式溶剂（如水分子）主要负责：真实溶剂化结构、氢键网络、局域溶剂化等
* 隐式溶剂主要负责：屏蔽净电场、提供远场介电响应

在典型 SJM 设置中：

* 隐式溶剂位于显式水层之上（可能也会对显式水分子提供一定稳定）
* 但隐式溶剂不应进入真正的反应区域，否则会出现“重复溶剂化（double-solvating）”

如果隐式溶剂不小心渗入了不希望的区域，可以通过加入 **ghost atoms** 来把溶剂排除（exclude）出特定空间区域。

---

### 恒内电势（CIP）DFT：Constant Inner Potential (CIP) DFT

除基于功函数（work function）控制电势外，SJM 也可以运行在 **constant inner potential（CIP）** 模式：控制电极内部的平均静电势（inner potential）。相关概念见文档中的 `The electrode potential` 部分。

CIP-DFT 的用法与标准 SJM 类似，但需要在 `sj` 字典中指定参考标度与 `cip` 子字典：

```python
sj_cip = {'target_potential': pot, # potential on the absolute potential scale
          'pot_ref': 'CIP',
          'cip': {'autoinner': {'nlayers': 4},
                  'mu_pzc': mu_pzc, # Fermi level at zero charge
                  'phi_pzc': phi_pzc # inner potential at zero charge,
                  }
          }
```

关键点解释：

* `pot_ref: 'CIP'`：把电势参考定义为 inner potential 相关标度
* `autoinner`：自动识别“电极内部区域（inner region）”，此处设置为 4 层金属 slab
* `mu_pzc` / `phi_pzc`：在 **PZC（potential of zero charge，零电荷电势）** 下的费米能级与内电势基准
* `target_potential`：在“绝对电势标度（absolute potential scale）”下的目标电势

原文给出一个与 **SHE（standard hydrogen electrode，标准氢电极）** 的换算例子：

* 若 PZC 为 **4.44 V vs SHE**，则 `mu_pzc = -4.44`
* 若目标电势为 **0.6 V vs SHE**，则 `target_potential = 3.84`

> **实用提醒**
> 这里的符号约定与“电子能量/电势参考”的定义相关，使用时建议保持同一套标度并在项目内统一记录（尤其在对接实验电势 vs SHE/RHE 时）。

---

#### 如何获得 `mu_pzc` 和 `phi_pzc`（推荐做法）

文档建议最简单的方法是先做一次“校准（calibration）”计算：

```python
sj_calib = {'excess_electrons':0,
    'pot_ref': 'CIP',
    'cip': {'autoinner': {'nlayers': 4,
            'threshold': 0.01}
            }}
calc = SJM(sj=sj_calib...)
atoms.calc = calc
atoms.get_potential_energy()
phi_pzc = calc.get_inner_potential(atoms)
mu_pzc = calc.get_fermi_level()
```

含义：

* `excess_electrons: 0`：零额外电子（对应 PZC 相关基准态）
* `threshold`：用于自动识别 inner 区域的判据阈值（细节依实现）
* 运行后用：

  * `calc.get_inner_potential(atoms)` 得到 `phi_pzc`
  * `calc.get_fermi_level()` 得到 `mu_pzc`

> **CIP-DFT 实用检查清单**
>
> * `autoinner` 是否正确识别了电极内部（比如层数 nlayers 是否合适）
> * 你的 slab 是否足够厚（否则“内区平均势”可能不稳定）
> * 真空与溶剂区域设置是否与目标电化学界面一致
> * 在做反应路径/NEB 时，建议先用第一张图像完成电势平衡，再启动整个轨迹计算（符合文档对计算成本的描述）

# 教程与练习（Tutorials and exercises）

本节翻译并整理 GPAW 文档中 **前两个主题**：

1. **Calculation of atomization energies**（原子化能/解离能的计算）
2. **Structure optimization**（结构优化：分子键长、晶格常数、应力张量与全弛豫示例）

为了便于直接用于材料计算中文参考文档，本文在翻译基础上做了结构化整理，并补充一些常见“坑点”与实用建议。代码块保持原样（必要处仅做解释，不随意改动 API 名称与参数）。

---

## 原子化能的计算（Calculation of atomization energies）

> **Warning（重要警告）**
> 主流 DFT（尤其是常见的 LDA/GGA）对“孤立原子（isolated atoms）”的电子态描述往往不够准确，**尤其是过渡金属原子（transition metals）**。实际表现通常是：SCF 很难收敛、或收敛到不物理的自旋态。
> 因此做原子化能/结合能时，要格外注意：**自旋态、初始磁矩、真空大小、收敛参数**。

### 目标：计算 H₂ 的原子化能

原子化能（atomization energy）常用定义为：

[
E_\text{atomization} = \sum_i E(\text{atom}_i) - E(\text{molecule})
]

对 H₂ 来说：

[
E_\text{atomization} = 2E(\text{H}) - E(\text{H}_2)
]

下面脚本将计算氢原子与氢分子的总能，并输出原子化能。

---

### 示例脚本：H₂ 原子化能（atomization.txt）

```python
# web-page: atomization.txt

from ase import Atoms
from ase.parallel import paropen as open
from gpaw import GPAW, PW

a = 10.  # Size of unit cell (Angstrom)
c = a / 2

# Hydrogen atom:
atom = Atoms('H',
             positions=[(c, c, c)],
             magmoms=[0],
             cell=(a, a + 0.0001, a + 0.0002))  # break cell symmetry

# gpaw calculator:
calc = GPAW(mode=PW(),
            xc='PBE',
            hund=True,
            txt='H.out')
atom.calc = calc

e1 = atom.get_potential_energy()
calc.write('H.gpw')

# Hydrogen molecule:
d = 0.74  # experimental bond length
molecule = Atoms('H2',
                 positions=([c - d / 2, c, c],
                            [c + d / 2, c, c]),
                 cell=(a, a, a))

calc = calc.new(hund=False,  # no hund rule for molecules
                txt='H2.out')

molecule.calc = calc
e2 = molecule.get_potential_energy()
calc.write('H2.gpw')

with open('atomization.txt', 'w') as fd:
    print(f'  hydrogen atom energy:     {e1:5.2f} eV', file=fd)
    print(f'  hydrogen molecule energy: {e2:5.2f} eV', file=fd)
    print(f'  atomization energy:       {2 * e1 - e2:5.2f} eV', file=fd)
```

---

### 脚本逐段解释（实用导读）

1. **构造氢原子 `Atoms` 对象**

* 把单个 H 放进一个很大的立方盒子（`a=10 Å`），用真空把原子“隔离开”
* `cell=(a, a+0.0001, a+0.0002)`：刻意轻微破坏晶胞对称性（break cell symmetry）
  这在某些情况下可以避免对称性导致的简并/数值问题（尤其对原子、开壳层体系更常见）。

> 代码里 `magmoms=[0]` 看似设置为非磁性，但这里同时启用了 `hund=True`。
> 在 GPAW 中 `hund=True` 会用 Hund’s rule 给出自旋极化的初始占据/磁矩，并且**会覆盖用户给的磁矩设置**（这是常见用法，尤其对原子更稳）。

2. **设置 GPAW 计算器**

* `mode=PW()`：平面波（plane-wave）模式
* `xc='PBE'`：PBE 交换相关泛函
* `hund=True`：对原子使用 Hund 规则给出合理的开壳层自旋初猜（对很多元素非常重要）
* `txt='H.out'`：输出写入文件

3. **计算原子能量并保存重启文件**

* `e1 = atom.get_potential_energy()` 得到氢原子能量
* `calc.write('H.gpw')` 保存波函数、密度等（以后可 restart）

4. **构造氢分子并复用计算器参数**

* `d=0.74 Å` 使用实验键长构型（只是起点）
* `calc = calc.new(hund=False, txt='H2.out')`：
  从原子计算器“复制并修改少量参数”得到新的计算器：对分子关闭 Hund 规则，输出写到 `H2.out`。

5. **输出原子化能**

* 写到 `atomization.txt`（用 `ase.parallel.paropen` 适配 MPI 并行安全写文件）

---

### 运行结果示例与解释

运行脚本后，输出可能类似：

```text
  hydrogen atom energy:     -1.08 eV
  hydrogen molecule energy: -6.65 eV
  atomization energy:        4.50 eV
```

解释：

* 原子化能 `4.50 eV` 表示 H₂ 相对于 2 个孤立 H 原子更稳定（能量降低 4.50 eV）。
* 文档中引用 Blaha 等人的全电子（all-electron）PBE 计算给出约 `4.54 eV`，与该 PAW 结果非常一致，说明该示例在数值上是合理的。

---

### 一个关键现象：为什么非自旋极化氢原子能量接近 0？

文档指出：自旋极化的氢原子能量约 `-1.09 eV`。如果对原子用 `magmom=0` 且 `hund=False`（即强制非自旋极化），会得到接近 `0 eV` 的结果。

原因（非常重要，实用性强）：

* GPAW 的总能往往相对于“参考原子（reference atom）”定义；参考原子通常是**中性、成对自旋（spin-paired）、球对称**的状态（也是生成 setup 的参考）。
* 对 H 这种简单体系，非自旋极化的 H 原子能量就是参考能，理论上在 **非常大盒子 + 足够高数值精度（更密网格/更高 ecut）** 下应精确趋近于 0。

> **实用提醒**
> 做原子化能时一定要保证“原子能量计算”的自旋态正确，否则会出现非常离谱的能量（特别是过渡金属，错误自旋态会导致 eV 级误差）。

---

### 原子能量计算的实用建议（强烈推荐）

* **大真空盒子**：至少保证原子与自身周期像之间相距足够远（常用 8–12 Å 甚至更大，视体系而定）。
* **自旋设置**：

  * 原子通常应 `spinpol=True`（或通过 `hund=True` 触发）
  * 过渡金属原子建议显式设置合适的初始磁矩（或 `hund=True`），否则 SCF 很容易不收敛。
* **收敛与数值参数**：

  * PW 模式需收敛 `ecut`；FD 模式需收敛 `h`；LCAO 需收敛 basis
* **对称性**：原子/分子偶尔会因对称性导致收敛怪异，破坏晶胞对称（如本例）或 `symmetry='off'` 有时更稳（但会更贵）。

---

## 结构优化（Structure optimization）

本主题包含：分子键长优化、晶格常数寻找、PW 模式下应力张量与晶胞优化，以及更通用的结构弛豫“配方”（recipes）。

---

### 分子结构优化示例：H₂ 键长优化

在上一节原子化能示例中，我们用实验键长 `0.74 Å` 计算 H₂ 的能量。现在我们让 **ASE optimizer（结构优化器）** 自动寻找能量最低的几何构型，使所有原子力（forces）都小于 `0.05 eV/Å`。

> **Note**
> 你必须先运行原子化能脚本，生成 `H2.gpw`，下面脚本才能 restart。

---

#### 示例脚本：用 QuasiNewton 优化 H₂（optimization.txt）

```python
# web-page: optimization.txt

from gpaw import restart
from ase.parallel import paropen as open
from ase.optimize import QuasiNewton


molecule, calc = restart('H2.gpw', txt='H2-relaxed.txt')

e2 = molecule.get_potential_energy()
d0 = molecule.get_distance(0, 1)

with open('optimization.txt', 'w') as fd:
    print('experimental bond length:', file=fd)
    print(f'hydrogen molecule energy: {e2:5.2f} eV', file=fd)
    print(f'bondlength              : {d0:5.2f} Ang', file=fd)

    # Find the theoretical bond length:
    relax = QuasiNewton(molecule, logfile='qn.log')
    relax.run(fmax=0.05)

    e2 = molecule.get_potential_energy()
    d0 = molecule.get_distance(0, 1)

    print(file=fd)
    print('PBE energy minimum:', file=fd)
    print(f'hydrogen molecule energy: {e2:5.2f} eV', file=fd)
    print(f'bondlength              : {d0:5.2f} Ang', file=fd)
```

---

#### 结果示例与解释

```text
experimental bond length:
hydrogen molecule energy: -6.65 eV
bondlength              :  0.74 Ang

PBE energy minimum:
hydrogen molecule energy: -6.65 eV
bondlength              :  0.75 Ang
```

解释：

* PBE 优化得到的键长约 `0.75 Å`，与实验 `0.74 Å` 很接近。
* 能量显示相同到两位小数是因为打印格式较粗（`{e2:5.2f}`）；如果用更多小数位会看到细微差别。

---

#### 约束（constraints）：只动一个原子以节省时间

如果你想节省优化时间，可以固定其中一个原子，只移动另一个原子。文档中给出约束示例：

```python
molecule.set_constraint(FixAtoms(mask=[0, 1]))
```

这里需要做一个**非常实用的澄清**（避免踩坑）：

* 在 ASE 中 `FixAtoms` 常用两种方式：

  1. 用 `indices=[...]` 指定要固定的原子序号
  2. 用 `mask=[True/False,...]`（布尔数组）指定每个原子是否固定

因此更常见、也更清晰的写法例如：

```python
from ase.constraints import FixAtoms

# 固定第 0 号原子
molecule.set_constraint(FixAtoms(indices=[0]))

# 或者用 mask（长度等于原子数）
molecule.set_constraint(FixAtoms(mask=[True, False]))
```

文档原句强调：`mask` 是“每个原子对应一个布尔值”，用于指示是否固定。更多约束示例见 `ase.constraints` 模块。

---

### 寻找晶格常数（Finding lattice constants）

晶格常数（lattice constant）的常见求法是：在一系列体积/晶格缩放下计算总能，然后拟合能量-体积曲线（例如用抛物线或 EOS 方程）。为了可靠，必须对数值参数（`ecut`、`kpts`、smearing 等）做收敛测试。

#### 示例：fcc 铝（Fcc Aluminium）

先对平面波截断能（plane-wave cutoff, `ecut`）做收敛测试：

```python
import numpy as np
from ase.build import bulk
from gpaw import GPAW, PW

a0 = 4.04
al = bulk('Al', 'fcc', a=a0)
cell0 = al.cell.copy()

for ecut in range(200, 501, 50):
    al.calc = GPAW(mode=PW(ecut),
                   xc='PBE',
                   kpts=(8, 8, 8),
                   basis='dzp',
                   txt=f'Al-{ecut}.txt')
    for eps in np.linspace(-0.02, 0.02, 5):
        al.cell = (1 + eps) * cell0
        al.get_potential_energy()
```

解释（实用导读）：

* 外层循环：`ecut` 从 200 到 500 eV，每 50 eV 一步
* 内层循环：把晶胞按比例 `(1+eps)` 缩放（±2%），计算每个体积下的总能
* 这样你可以从输出中判断：在你关心的体积范围内，总能随 `ecut` 是否已收敛到你需要的精度（比如 meV/atom）

然后固定一个足够大的 `ecut`（例如 400 eV）对 k 点做收敛：

```python
calc = al.calc.new(mode=PW(400))
for k in range(4, 17):
    al.calc = calc.new(kpts=(k, k, k),
                       txt=f'Al-{k:02}.txt')
    for eps in np.linspace(-0.02, 0.02, 5):
        al.cell = (1 + eps) * cell0
        al.get_potential_energy()
```

> **实用扩展**
> 上述脚本只是“计算数据”。实际确定晶格常数通常会把 E(V) 数据拟合：
>
> * 简单：用抛物线拟合极小点
> * 更标准：用 `ase.eos.EquationOfState` 拟合 Birch–Murnaghan 等 EOS
>   这一步非常推荐加入到你自己的工作流中（尤其是材料性质计算）。

---

### 平面波模式与应力张量（Plane wave mode and Stress tensor）

实空间网格（real-space grid, FD）模式的主要优势是：大体系并行效率高。但对小体系，平面波（plane-wave, PW）基组往往更快。

PW 模式下：

* 所有量用周期超胞（periodic supercell）上的 Fourier 表示
* 必须使用周期边界条件（periodic boundary conditions）
* 关键收敛参数：

  * FD：网格间距 `h`
  * PW：截断能 `E_cut`，对应最大倒格矢截断 `G_cut`
    [
    E_{cut} = \frac{G_{cut}^2}{2} \quad (\text{atomic units})
    ]

---

#### 平面波截断能收敛（Converging the plane wave cutoff）

示例：对体相 Si 测试总能随 `E_cut` 的收敛：

```python
from ase.build import bulk
from gpaw import GPAW, PW

a = 5.421
si = bulk('Si', 'fcc', a=a)
# or equivalently:
# b = a / 2
# from ase import Atoms
# si = Atoms('Si2', cell=[[0, b, b], [b, 0, b], [b, b, 0]], pbc=True,
#           scaled_positions=[[0, 0, 0], [0.25, 0.25, 0.25]])

for x in [100, 200, 300, 400, 500, 600, 700, 800]:
    # for x in [0.24, 0.22, 0.20, 0.18, 0.16, 0.14, 0.12, 0.1]:
    calc = GPAW(mode=PW(x),
                # h=x,
                xc='PBE',
                kpts=(4, 4, 4),
                txt=f'convergence_{x}.txt')

    si.calc = calc

    print(x, si.get_potential_energy())
```

文档说明：

* 这里用较粗的 k 点 `(4,4,4)` 和默认费米展宽 `0.1 eV`（Fermi smearing），只是为了加速测试。
* 这些参数在正式计算中也应分别收敛，但此处先固定不变。
* 在 `E_cut = 400 eV` 时，总能通常可收敛到 ~5 meV 以内（示例结论）。

**练习（原文问题翻译）**
请修改脚本，在网格（grid）模式下用不同 `h` 计算（提示：注释掉 `mode=PW(x)` 相关，改用 `h=x` 那套循环），回答：

1. 要把总能收敛到 5 meV 内，需要多小的 `h`？
2. 与 PW 模式 400 eV 相比，计算耗时如何？
3. 比较两种计算使用的网格点数与平面波数（分别查看 `convergence_400.txt` 与对应网格计算的输出文件）
4. `.txt` 文件里会记录时间信息，可用来对比耗时

> **实用提示**
>
> * 做“总能收敛”比较时，要确保各计算之间其它参数一致（k 点、展宽、收敛阈值、赝势/PAW setup 等）。
> * 对材料性质（晶格常数、弹性等）更推荐比较**能量差**或 EOS 拟合结果对参数的敏感性，而不是只盯单点总能。

---

#### 晶胞优化（Optimizing the unit cell）

> **Warning（重要警告）**
> 由于同时优化晶胞与原子位置存在数值与算法困难，`ase.filters.UnitCellFilter` 有时可能给出不正确结果。
> **务必**用以下方式验证：
>
> * 分开做晶胞优化（例如 `ase.filters.StrainFilter`）
> * 再做原子位置优化（`ase.optimize`）
>   并且建议使用比本教程更严格的 `fmax`（本教程用 0.05 eV/Å 偏松，仅用于演示）。

在前面的铝例子里，我们通过扫描不同晶格常数并拟合能量找到最优点。PW 模式的一个重要优势是可计算 **stress tensor（应力张量）**，因此可以直接用应力来优化周期体系的晶胞参数。

下面脚本对体相 Si 做晶胞优化：

```python
import numpy as np
from ase.build import bulk
from ase.optimize.bfgs import BFGS
from ase.filters import UnitCellFilter
from gpaw import GPAW
from gpaw import PW

si = bulk('Si', 'fcc', a=6.0)
# Experimental Lattice constant is a=5.421 A

si.calc = GPAW(xc='PBE',
               mode=PW(400, dedecut='estimate'),
               kpts=(4, 4, 4),
               # convergence={'eigenstates': 1.e-10},  # converge tightly!
               txt='stress.txt')

uf = UnitCellFilter(si)
relax = BFGS(uf)
relax.run(fmax=0.05)  # Consider much tighter fmax!

a = np.linalg.norm(si.cell[0]) * 2**0.5
print(f'Relaxed lattice parameter: a = {a} Ang')
```

补充说明（更实用）：

* `dedecut='estimate'` 通常用于与有限 `ecut` 相关的应力修正（常被称为 Pulay stress 相关处理），帮助提高 stress 的可靠性。
* 如果你想可视化优化过程，建议在优化器里加 `trajectory='stress.traj'`，然后用：

  * `ase gui stress.traj`
    来查看结构演化（比直接看 `stress.txt` 更稳妥）。

文档提到该计算大约用 12 次迭代找到最优晶格常数。
对于已知实验晶格常数的简单晶胞，用 5 个点拟合抛物线可能更快；但对多参数晶胞（例如三斜晶、多个晶格参数耦合），stress tensor 会非常有价值。

---

### 结构弛豫“配方”（Structure relaxation recipes）

前面我们介绍了结构弛豫的基础：力驱动的几何优化、体积扫描、stress 优化等。这里给出更可复用的脚本框架，展示如何在 GPAW 中组织结构优化流程，并逐步提高精度。

---

#### 固定晶胞的弛豫（Fixed-cell relaxation）

最简单情况：**只优化原子位置**，晶胞形状/体积固定。做法：

1. 用 GPAW 计算 forces
2. 交给 ASE 优化器（如 `BFGS`）迭代更新原子位置，直至 `|F| < fmax`

先导入模块：

```python
from ase.build import bulk
from ase.optimize import BFGS
from gpaw import GPAW
```

下面是一个通用 `relax()` 函数（原文代码），包含：

* 设置 DFT 计算器参数
* 设置初始磁矩（可选）
* 可选加入 **DFT-D3**（van der Waals）修正
* 可选择只优化位置（fixcell=True）或进行全弛豫（fixcell=False，使用 cell filter）
* 输出 logfile 与 trajectory

```python
def relax(atoms, calculator_params,
          fmax=0.01, d3=False, fixcell=True,
          logname='opt.log',
          trajname='opt.traj'):

    # set DFT calculator
    calc_dft = GPAW(**calculator_params)

    # magnetize atoms
    atoms.set_initial_magnetic_moments(len(atoms) * [1])
    # non-magnetic calculation:
    # atoms.set_initial_magnetic_moments(len(atoms) * [0])

    # optionally include van der Waals DFT-D3
    if d3:
        from ase.calculators.dftd3 import DFTD3
        calc = DFTD3(dft=calc_dft)
    else:
        calc = calc_dft

    # set calculator
    atoms.calc = calc

    # set configuration to be optimized
    if fixcell:
        # only optimize positions of the atoms
        opt_conf = atoms
    else:
        # setup full relaxation
        # set unit cell filter
        opt_conf = FrechetCellFilter(atoms)

    # setup optimizer
    # specify logfile and trajectory file names
    opt = BFGS(opt_conf, logfile=logname, trajectory=trajname)
    # run the optimization until forces are smaller than fmax
    opt.run(fmax=fmax)

    return atoms
```

输出：

* 历史写到 `opt.log`
* 轨迹写到 `opt.traj`
* 函数返回最终构型 `atoms`

---

#### “先快后准”的参数组织方式（fast → refined）

文档建议把计算参数写到统一格式的参数文件中（如 `params_forces.py`），便于从粗略快速优化逐步收敛到高精度。

下面是一个“快速 forces”参数集（原文）：

```python
# fast forces
calculator_params = {
    "xc": "PBE",
    "basis": "dzp",
    "mode": {"name": "lcao"},
    "kpts": {"size": [1, 1, 1],
             "gamma": True},
    "convergence": {"density": 1e-3,
                    "forces": 1e-2},
    "occupations": {"name": "fermi-dirac",
                    "width": 0.05},
    "mixer": {"method": "fullspin",
              "backend": "pulay"},
    "txt": "rlx.txt",
}
```

提升精度的建议（原文）：

* 更密的 k 点采样（每个方向 5 个或更多）
* 更严格的收敛：例如
  `"convergence": {"density": 1e-6, "forces": 1e-4}`

> **实用提醒**
>
> * LCAO 快速结构优化很常用，但最终能量/力的“绝对精度”需要 FD 或 PW 再验证（取决于你的研究目标）。
> * 对金属体系，展宽（smearing）与 k 点收敛对力的稳定性影响很大。

---

#### 全弛豫（Full relaxation：同时优化晶胞与原子位置）

若要同时优化晶胞形状/体积与原子位置，可以用 cell filters（例如 `FrechetCellFilter`）：

* 需要先导入：`from ase.filters import FrechetCellFilter`
* 把 `Atoms` 包装成 filter 后直接传给优化器

在上面的 `relax()` 中，把：

```python
opt_conf = atoms
```

替换为：

```python
opt_conf = FrechetCellFilter(atoms)
```

即可。

原文强调：**应力（stress）计算需要 GPAW 的 PW 模式**，并且应使用更高精度参数（例如 `params_stresses.py`）。示例“accurate stress”参数集：

```python
# accurate stress
calculator_params = {
    "xc": "PBE",
    "basis": "dzp",
    "mode": {"name": "pw",
             "ecut": 600},
    "kpts": {"size": [5, 5, 5],
             "gamma": True},
    "convergence": {"density": 1e-6,
                    "forces": 1e-4},
    "occupations": {"name": "fermi-dirac",
                    "width": 0.05},
    "mixer": {"method": "fullspin",
              "backend": "pulay"},
    "txt": "rlx.txt",
}
```

文档最后提到可下载完整脚本 `relax.py`（此处保留原意：你可以把以上片段整理成完整脚本用于项目模板）。

---

### 小结：本节最重要的实践要点

* 原子化能计算中，“孤立原子”是数值最难的部分：**自旋态与收敛是关键**。
* 分子结构优化通常用 ASE 优化器（QuasiNewton/BFGS 等），收敛标准 `fmax` 要根据任务调整。
* 晶格常数建议通过 **收敛测试 + EOS 拟合**获得；PW 应力张量可直接做晶胞优化，但要警惕 filter 同时优化带来的误差并做交叉验证。
* 结构优化推荐“先快后准”：先用较松的参数快速靠近极小点，再用更严格参数精修与验证。


## Energetics（能量学）

本节涵盖 GPAW 中与**能量差**相关的一些典型工作流：内聚能（cohesive energy）、DFT+U、带电缺陷形成能（含 FNV 校正）、RPA 相关能，以及基于 TDDFT 的相关能（rALDA / rAPBE 等）。内容偏“实战”，会穿插一些容易踩坑的细节说明。

---

### 1. 体相 FCC Pt 的内聚能（Cohesive energy of bulk FCC Pt）

计算内聚能（cohesive energy）通常需要两类能量：

* **孤立原子能量** (E_\text{atom})（需要足够大真空、合理的自旋态）
* **体相能量** (E_\text{bulk})（通常要换算到 *每原子*）

一般定义（按 eV/atom）：
[
E_\text{cohesive} = E_\text{atom} - \frac{E_\text{bulk}}{N_\text{atoms in cell}}
]

> ⚠️ 实用提醒
> ASE 的 `bulk()` 返回的晶胞里**含几个原子**取决于结构/参数（是否 primitive/cubic 等）。因此在做内聚能时，务必检查 `len(atoms)` 并按需要除以原子数。下面示例直接 `e_atom - e_bulk`，隐含前提是体相能量是按 1 原子计的（或返回的是 1 原子原胞）。

对于某些元素（尤其过渡金属原子），**孤立原子的正确磁态**可能难以用常规 SCF（如 Davidson）稳定收敛。此时常用策略是使用 **Direct Minimization Methods（直接最小化）**，并配合固定占据、关闭混合等，使体系更稳定。

下面是 Pt 的示例（使用平面波 `PW` 模式 + ETDM 最小化）：

```python
from ase import Atoms
from ase.build import bulk
from gpaw import GPAW, GPAW_NEW

atom = Atoms('Pt')
atom.center(vacuum=6.0)
atom.calc = GPAW(
    xc='PBE',
    mode={'name': 'pw', 'ecut': 600.0},
    nbands=-2,
    mixer={'backend': 'no-mixing'},
    occupations={'name': 'fixed-uniform'},
    hund=True,
    eigensolver={'name': 'etdm-fdpw', 'converge_unocc': not GPAW_NEW},
    symmetry='off',
    txt='pt-atom.txt')
e_atom = atom.get_potential_energy()
atom.calc.write('pt-atom.gpw')

bulk = bulk('Pt', 'fcc', a=3.92)
k = 8
bulk.calc = GPAW(
    xc='PBE',
    mode={'name': 'pw', 'ecut': 600.0},
    kpts=(k, k, k),
    txt='pt-bulk.txt')
e_bulk = bulk.get_potential_energy()

e_cohesive = e_atom - e_bulk
print(e_cohesive, 'eV')
```

#### 1.1 如何检查“原子态是否算对了”？

强烈建议仔细查看 `pt-atom.txt`，确认：

* 收敛（SCF/ETDM 收敛信息）
* 自旋态/占据是否符合预期
* 费米能、占据是否稳定

下面给出一个更深入的分析脚本：它读取 `.gpw` 并对每个自旋通道的本征态做一个“轨道成分（character）”的粗略标注，输出到 CSV，便于你用 Excel/Origin 做进一步检查。

```python
from gpaw import GPAW
from gpaw.sphere.spherical_harmonics import names
from ase.units import Ha

calc = GPAW('pt-atom.gpw')
setup = calc.setups[0]
labels = [f'{n}{"spd"[l]}({names[l][m]})'
          for n, l in zip(setup.n_j, setup.l_j)
          for m in range(2 * l + 1)]

lines = ['#,eig. [eV],occ,character,eig. [eV],occ,character']
for n in range(11):
    line = str(n)
    for spin in [0, 1]:
        kpt = calc.wfs.kpt_qs[0][spin]
        P_i = kpt.P_ani[0][n]
        i = abs(P_i).argmax()
        label = labels
        eig = kpt.eps_n[n] * Ha
        occ = kpt.f_n[n]
        line += f',{eig:.3f},{occ:.1f},{labels[i]}'
    lines.append(line)

# Write csv-file:
with open('pt-atom.csv', 'w') as fd:
    print('\n'.join(lines), file=fd)
```

> ✅ 实用建议
> 对“难算的孤立原子”，常见稳态组合包括：
>
> * `symmetry='off'`（避免对称性约束误导）
> * `mixer={'backend': 'no-mixing'}` + `occupations={'name': 'fixed-uniform'}`
> * `hund=True`（按 Hund’s rule 猜初始磁矩/占据）
> * 直接最小化（如 `etdm-fdpw`）
>   同时，真空要足够大（例如 6 Å 甚至更大），并避免晶胞完美对称导致的偶发问题（有时会轻微打破 cell 对称性）。

---

### 2. DFT+U 理论与 GPAW 实现（DFT+U theory）

#### 2.1 为什么需要 DFT+U？

主流 LDA/GGA 往往不能正确描述**局域电子（localized electrons）**的强 onsite 相互作用，典型是过渡金属的 d 电子、稀土的 f 电子（有时 p 轨道也会）。DFT+U 通过额外的 Hubbard-like 项补偿这部分缺失。

常见两大分支：

* Liechtenstein 等人的形式：U 与 J 分开进入（更一般）[Liechtenstein]
* Dudarev / Anisimov 的形式：只用一个有效参数
  [
  U_\text{eff} = U - J
  ]
  GPAW 实现的是这一支 [Dudarev]

在 GPAW 中，DFT+U 的总能量写为：
[
E_\text{DFT+U} = E_\text{DFT} + \sum_a \frac{U_\text{eff}}{2}
\mathrm{Tr}(\rho^a - \rho^a \rho^a),
]
其中 (\rho^a) 是原子轨道占据矩阵（atomic orbital occupation matrix）。这个项可理解为对 (\rho) 施加惩罚项，使其更趋近幂等（idempotent），即偏向“满占据/空占据”的极限。

#### 2.2 GPAW 中如何设置 U（最常见写法）

在 GPAW 里，通常通过 `setups` 关键字给某元素指定 U：

* 对 Mn 的 d 轨道加 (U_\text{eff}=6) eV：

  ```python
  setups={'Mn': ':d,6.0'}
  ```
* 对 N 的 p 轨道加 (U_\text{eff}=6) eV：

  ```python
  setups={'N': ':p,6.0'}
  ```

#### 2.3 示例：NiO 中对 Ni-d 加 U

```python
from ase import Atoms
from gpaw import GPAW, PW, FermiDirac

# Setup up bulk NiO in an antiferromagnetic configuration:
a = 4.19  # lattice constants
b = a / 2**0.5
m = 2.0
atoms = Atoms('Ni2O2',
              pbc=True,
              cell=(b, b, a),
              positions=[(0, 0, 0),
                         (b / 2, b / 2, a / 2),
                         (0, 0, a / 2),
                         (b / 2, b / 2, 0)],
              magmoms=[m, -m, 0, 0])

k = 2  # number of k-points
atoms.calc = GPAW(mode=PW(400),
                  occupations=FermiDirac(width=0.05),
                  setups={'Ni': ':d,6.0'},  # U=6 eV for Ni d orbitals
                  txt='nio.txt',
                  kpts=(k, k, k),
                  xc='PBE')
e = atoms.get_potential_energy()
```

> ✅ 实用建议
> DFT+U 对能带/磁矩/带隙很敏感。建议：
>
> * 明确磁有序（如反铁磁）与初始磁矩 `magmoms`
> * 金属体系或有部分占据时，smearing（`FermiDirac(width=...)`）要收敛测试
> * U 值建议有来源（cRPA、linear-response 或文献/经验），并做敏感性分析。

---

#### 2.4 Hubbard 投影归一化（Normalization）与 GPAW/VASP 差异

占据矩阵 (\rho^a) 的投影在 GPAW 中通常在 augmentation sphere 半径内截断（setup 的投影子定义导致），因此投影一般小于 1。

此时有两种选择：

* **归一化（normalized）**：把投影按 augmentation sphere 内的积分归一
* **不归一化（not normalized）**：直接用截断后的投影

GPAW 默认是 **归一化**。但一些其他 PAW 代码（例如 VASP）常见设置是不归一化。GPAW 允许关闭归一化：在 setup 字符串末尾加 `,0`。

例如 N(p) 上 U=6 eV，关闭归一化：

```python
setups={'N': ':p,6.0,0'}
```

对 d/f 轨道影响通常不大（它们绝大部分波函数在 augmentation sphere 内），但对 p 轨道可能差异显著。

下面给出一个极端示例：原子 N 的 p 轨道在自旋分裂上，对“归一化/不归一化”的差异非常明显。

```python
from ase import Atoms
from gpaw import GPAW

n = Atoms('N', magmoms=[3])
n.center(vacuum=3.5)

# Calculation with no +U correction:
n.calc = GPAW(mode='lcao',
              basis='dzp',
              txt='no_u.txt',
              xc='PBE')
e1 = n.get_potential_energy()
n.calc.write('no_u.gpw')

# Calculation with a correction U=6 eV normalized:
n.calc = n.calc.new(setups={'N': ':p,6.0'}, txt='normalized_u.txt')
e2 = n.get_potential_energy()
n.calc.write('normalized_u.gpw')

# Calculation with a correction U=6 eV not normalized:
n.calc = n.calc.new(setups={'N': ':p,6.0,0'}, txt='not_normalized_u.txt')
e3 = n.get_potential_energy()
n.calc.write('not_normalized_u.gpw')
```

给出的 2p splitting（单位 eV）：

| no U             | 4.232  |
| ---------------- | ------ |
| normalized U     | 11.010 |
| not normalized U | 4.772  |

---

#### 2.5 对同一元素设置多个 U（多轨道同时修正）

可以给同一元素的不同轨道同时加 U，例如：

```python
setups={'Ni':':d,4.0,0;p,2.0,0'}
```

> ⚠️ Warning
> 该功能“理论上可用”，但实际应用测试较少；多数代码不支持同元素多 U 修正，因此横向对比参考可能有限。

---

### 3. 带电缺陷形成能：FNV 校正示例（Calculating the formation energies of charged defects）

带电缺陷（charged defects）在有限超胞中会产生显著的**长程库仑有限尺寸效应**，导致能量随超胞尺寸收敛很慢。FNV（Freysoldt–Neugebauer–Van de Walle）类方案通过：

* **模型电荷（model charge）**与介电屏蔽
* **势能对齐（potential alignment）**

来修正这一误差。

#### 3.1 示例：GaAs 中的 Ga 空位（(V_\mathrm{Ga}^{-3})）

我们考虑 (q=-3) 的 Ga vacancy（Ga 空位），使用 NxNxN 的 GaAs 超胞（含 (8N^3) 个原子）。先计算 pristine 与 defect 的总能，并保存 `.gpw` 供后续电静力校正使用。

```python
import sys
from ase.build import bulk
from gpaw import GPAW
from pathlib import Path

# Script to get the total energies of a supercell
# of GaAs with and without a Ga vacancy

N = int(sys.argv[1])  # NxNxN supercell
label = f'GaAs_{N}x{N}x{N}'
prs_path = Path(f'{label}_prs.gpw')
def_path = Path(f'{label}_def.gpw')

a0 = 5.628      # lattice parameter
charge = -3     # defect charge

params = {'mode': {'name': 'pw', 'ecut': 400},
          'xc': 'LDA',
          'kpts': {'size': (2, 2, 2), 'gamma': False},
          'occupations': {'name': 'fermi-dirac', 'width': 0.01}}

calc_charged = GPAW(charge=charge, **params)
calc_neutral = GPAW(charge=0, **params)

prim = bulk('GaAs', crystalstructure='zincblende', a=a0, cubic=True)
pristine = prim * (N, N, N)
pristine.calc = calc_neutral
pristine.get_potential_energy()
pristine.calc.write(prs_path)

defect = pristine.copy()
defect.pop(0)  # make a Ga vacancy
defect.calc = calc_charged
defect.get_potential_energy()
defect.calc.write(def_path)
```

未校正的能量差：
[
(E[X^q] - E_0)_\mathrm{uncorrected}
]
原文给出示例值为 21.78 eV（以 N=2 为例）。

#### 3.2 计算 FNV 校正：`charged_defect_corrections()`

`charged_defect_corrections()` 提供了对 `ElectrostaticCorrections()` 的高层接口。你需要提供：

* `epsilon`：介电常数（这里取 12.7，为 clamped-ion static limit）
* 模型电荷：高斯（Gaussian）模型电荷的位置与宽度

  * 例：中心在 (0,0,0)，标准差 0.72 Bohr（FWHM=2 Bohr）

```python
from ase.io.jsonio import write_json
from gpaw import GPAW
from gpaw.defects import charged_defect_corrections
from pathlib import Path

charge = -3
epsilon = 12.7
def_idx = 0
corrected = []
uncorrected = []
repeats = [1, 2, 3, 4]
for N in repeats:
    label = f'GaAs_{N}x{N}x{N}'
    prs_path = Path(f'{label}_prs.gpw')
    def_path = Path(f'{label}_def.gpw')

    calc_prs = GPAW(prs_path)
    calc_def = GPAW(def_path)

    elc = charged_defect_corrections(calc_pristine=calc_prs,
                                     calc_defect=calc_def,
                                     defect_index=def_idx,
                                     charge=charge,
                                     epsilon=epsilon)
    E_fnv = elc.calculate_correction()

    E_0 = calc_prs.get_potential_energy()
    E_X = calc_def.get_potential_energy()
    E_uncorr = E_X - E_0
    E_corr = E_uncorr + E_fnv

    if N == 2:
        profile = elc.calculate_potential_profile()

    corrected.append(E_corr)
    uncorrected.append(E_uncorr)

res = {'repeats': repeats, 'corrected': corrected,
       'uncorrected': uncorrected, 'profile': profile}

write_json('electrostatics.json', res)
```

该脚本会生成 `electrostatics.json`，其中包含：

* (\Delta V(z))（势能对齐函数）
* 模型势、以及 defect 与 pristine 的 planar average 势差

可用下面脚本绘图（调用 `plot_potentials()`）：

```python
# web-page: planar_averages.png
from ase.io.jsonio import read_json
from gpaw.defects.electrostatics import plot_potentials

data = read_json('electrostatics.json')

# obtain potential profile (only generated for N=2)
profile = data['profile']
plot_potentials(profile, png='planar_averages.png')
```

> ✅ 一致性检查（非常重要）
> 在缺陷远离区域（通常是晶胞中间）(\Delta V(z)) 应该趋于**平坦**。若不平坦，可能意味着：
>
> * 超胞太小/真空不足/介电模型不适用
> * 模型电荷中心或宽度不合适
> * 缺陷结构未充分弛豫、或存在长程极化未捕捉

原文给出示例：提取 (\Delta V=-0.14) eV，且 (E_\mathrm{lat}=-1.28) eV，于是：
[
E[V_\mathrm{Ga}^{-3}] - E_0
= \left[21.78  - (-1.28 ) + (-3)\times(-0.14)\right]\mathrm{eV}
= 23.5\ \mathrm{eV}
]

#### 3.3 形成能还缺什么？（Zhang–Northrup 公式的其他项）

要得到“缺陷形成能（formation energy）”，还需：

1. **价带顶位置** (\epsilon_v)（pristine 计算中得到）

   * GPAW 默认把平均静电势设为 0，因此与超胞 pristine 的参考通常已对齐，不必再额外对齐（但仍需你确认 band edge 是否被 k 点采到，比如 GaAs 需要包含 Γ 点）。
2. **化学势** (\mu_i) 的取值范围

   * 例如 Ga vacancy 需要 (\mu_\mathrm{Ga})，受限于 Ga-rich / As-rich 条件：
     [
     \mu_\mathrm{Ga}[\mathrm{bulk\ Ga}] > \mu_\mathrm{Ga}[\mathrm{GaAs}]>
     \mu_\mathrm{Ga}[\mathrm{bulk\ Ga}] + \Delta H_\mathrm{f}[\mathrm{GaAs}]
     ]
   * (\mu_\mathrm{Ga}[\mathrm{bulk\ Ga}]) 与 (\Delta H_\mathrm{f}[\mathrm{GaAs}]) 可由 Ga、As、GaAs 的总能计算获得。

---

### 4. RPA 相关能（Calculating RPA correlation energies）

RPA（Random Phase Approximation）给出一种非局域的基态相关能表达。它通常与 **exact exchange（EXX）**结合用于：

* van der Waals 相互作用（色散力）
* 吸附能、层间结合等

但代价非常高：需要大量未占据态（unoccupied bands）与较高的平面波截断（plane-wave cutoff）。

下面给出几个典型示例。

---

#### 4.1 示例 1：(N_2) 的 RPA 解离能（Atomization energy）

##### (1) 基态计算（PBE）并准备大量空带

关键点：

* 需要存波函数：`calc.write(..., mode='all')`
* 需要对完整哈密顿量对角化得到更多空带：`diagonalize_full_hamiltonian()`
* 为了后续 non-self-consistent HF/EXX（EXX@PBE）与 RPA，常用 `symmetry='off'`

```python
from ase.optimize import BFGS
from ase.build import molecule
from ase.parallel import paropen
from gpaw import GPAW, PW
from gpaw.hybrids.energy import non_self_consistent_energy as nsc_energy

# N
N = molecule('N')
N.cell = (6, 6, 7)
N.center()
calc = GPAW(mode=PW(600, force_complex_dtype=True),
            symmetry='off',
            nbands=16,
            maxiter=300,
            xc='PBE',
            hund=True,
            txt='N_pbe.txt',
            parallel={'domain': 1},
            convergence={'density': 1.e-6})

N.calc = calc
E1_pbe = N.get_potential_energy()

calc.write('N.gpw', mode='all')

E1_hf = nsc_energy('N.gpw', 'EXX').sum()

# calc.diagonalize_full_hamiltonian(nbands=4800)
calc.diagonalize_full_hamiltonian(nbands=None)
calc.write('N.gpw', mode='all')

# N2
N2 = molecule('N2')
N2.cell = (6, 6, 7)
N2.center()
calc = GPAW(mode=PW(600, force_complex_dtype=True),
            symmetry='off',
            nbands=16,
            maxiter=300,
            xc='PBE',
            txt='N2_pbe.txt',
            parallel={'domain': 1},
            convergence={'density': 1.e-6})

N2.calc = calc
dyn = BFGS(N2)
dyn.run(fmax=0.05)
E2_pbe = N2.get_potential_energy()

calc.write('N2.gpw', mode='all')

E2_hf = nsc_energy('N2.gpw', 'EXX').sum()

with paropen('PBE_HF.dat', 'w') as fd:
    print('PBE: ', E2_pbe - 2 * E1_pbe, file=fd)
    print('HF: ', E2_hf - 2 * E1_hf, file=fd)

# calc.diagonalize_full_hamiltonian(nbands=4800)
calc.diagonalize_full_hamiltonian(nbands=None)
calc.write('N2.gpw', mode='all')
```

> ✅ 提醒
> 对角化全哈密顿量是重计算；你需要评估内存/并行策略。脚本里 `nbands=None` 表示让代码自行决定。

##### (2) 频率积分收敛（frequency integration）

原文描述为“2000 个频点到 1000 eV”，但下面示例代码实际是**200 个频点，步长 0.5 eV，到约 100 eV**。实际使用时你可以据需要改大频率范围/点数。

```python
import numpy as np
from ase.parallel import paropen

from gpaw.xc.rpa import RPACorrelation

dw = 0.5
frequencies = np.array([dw * i for i in range(200)])
weights = dw * np.ones(len(frequencies))
weights[0] /= 2
weights[-1] /= 2

rpa = RPACorrelation('N2.gpw',
                     txt='frequency_equidistant.txt',
                     frequencies=frequencies,
                     weights=weights,
                     ecut=[50])
data = rpa.calculate_all_contributions()
Es_w = data.energy_wi[:, 0]

with paropen('frequency_equidistant.dat', 'w') as fd:
    for w, E in zip(frequencies, Es_w):
        print(w, E, file=fd)
```

Gauss–Legendre 参数方式（常用默认）示意：

```python
rpa = RPACorrelation(calc,
                     nfrequencies=16,
                     frequency_max=800.0,
                     frequency_scale=2.0)
```

* `nfrequencies`：点数
* `frequency_max`：最高频率（积分上限近似，理论应到无穷）
* `frequency_scale`：控制 (\omega\to 0) 附近的采样密度（若费米能附近 DOS 非零，(\omega\approx 0) 常更“敏感”）

##### (3) 对截断（ecut）外推到无穷（无限空带/无限平面波极限）

按 `RPA correlation energy` 的常用做法，对不同 ecut 计算相关能并外推：

```python
from ase.parallel import paropen
from ase.units import Hartree
from gpaw.xc.rpa import RPACorrelation

rpa = RPACorrelation('N.gpw', nblocks=8,
                     txt='rpa_N.txt', ecut=400)
E1_i = rpa.calculate()

rpa = RPACorrelation('N2.gpw', nblocks=8,
                     txt='rpa_N2.txt', ecut=400)
E2_i = rpa.calculate()
ecut_i = rpa.integral.ecut_i

f = paropen('rpa_N2.dat', 'w')
for ecut, E1, E2 in zip(ecut_i, E1_i, E2_i):
    print(ecut * Hartree, E2 - 2 * E1, file=f)
f.close()
```

外推与作图示例：

```python
# web-page: extrapolate.png, N2-data.csv
import numpy as np
import matplotlib.pyplot as plt
from gpaw.utilities.extrapolate import extrapolate
from pathlib import Path

a = np.loadtxt('rpa_N2.dat')
ext, A, B, sigma = extrapolate(a[:, 0], a[:, 1], reg=3, plot=False)
plt.plot(a[:, 0]**(-1.5), a[:, 1], 'o', label='Calculated points')
es = np.array([e for e in a[:, 0]] + [10000])
plt.plot(es**(-1.5), A + B * es**(-1.5), '--', label='Linear regression')

t = [int(a[i, 0]) for i in range(len(a))]
plt.xticks(a[:, 0]**(-1.5), t, fontsize=12)
plt.axis([0., 150**(-1.5), None, -4.])
plt.xlabel('Cutoff energy [eV]', fontsize=18)
plt.ylabel('RPA correlation energy [eV]', fontsize=18)
plt.legend(loc='lower right')
plt.savefig('extrapolate.png')

pbe, hf = (-float(line.split()[1])
           for line in Path('PBE_HF.dat').read_text().splitlines())
rpa = -A
Path('N2-data.csv').write_text(
    'PBE, HF, RPA, HF+RPA, Experimental\n'
    f'{pbe:.2f}, {hf:.2f}, {rpa:.2f}, {hf + rpa:.2f}, 9.89\n')
```

原文总结表（单位 eV）：

| PBE   | HF   | RPA  | HF+RPA | Experimental |
| ----- | ---- | ---- | ------ | ------------ |
| 10.60 | 4.84 | 4.91 | 9.75   | 9.89         |

> ✅ 经验提醒
> RPA 的优势主要来自“非局域相关”（vdW 等）。对小分子解离能，RPA 不一定总是优于 PBE。

---

#### 4.2 示例 2：石墨烯在 Cu(111) 表面的吸附（graphene adsorption）

该例演示在不同距离 d 下计算总能曲线，并在 `xc='RPA'` 时通过：

* PBE 自洽得到轨道
* non-self-consistent EXX（`nsc_energy(...,'EXX')`）
* `RPACorrelation` 得到相关能
  从而获得 HF+RPA 总能。

```python
from pathlib import Path

from ase import Atoms
from ase.build import fcc111
from gpaw import GPAW, PW, FermiDirac, MixerSum, Davidson
from gpaw.hybrids.energy import non_self_consistent_energy as nsc_energy
from gpaw.mpi import world
from gpaw.xc.rpa import RPACorrelation

# Lattice parametes of Cu:
d = 2.56
a = 2**0.5 * d
slab = fcc111('Cu', a=a, size=(1, 1, 4), vacuum=10.0)
slab.pbc = True

# Add graphite (we adjust height later):
slab += Atoms('C2',
              scaled_positions=[[0, 0, 0],
                                [1 / 3, 1 / 3, 0]],
              cell=slab.cell)


def calculate(xc: str, d: float) -> float:
    slab.positions[4:6, 2] = slab.positions[3, 2] + d
    tag = f'{d:.3f}'
    if xc == 'RPA':
        xc0 = 'PBE'
    else:
        xc0 = xc
    slab.calc = GPAW(xc=xc0,
                     mode=PW(800),
                     basis='dzp',
                     eigensolver=Davidson(niter=4),
                     nbands='200%',
                     kpts={'size': (12, 12, 1), 'gamma': True},
                     occupations=FermiDirac(width=0.05),
                     convergence={'density': 1e-5},
                     parallel={'domain': 1},
                     mixer=MixerSum(0.05, 5, 50),
                     txt=f'{xc0}-{tag}.txt')
    e = slab.get_potential_energy()

    if xc == 'RPA':
        e_hf = nsc_energy(slab.calc, 'EXX').sum()

        slab.calc.diagonalize_full_hamiltonian()
        slab.calc.write(f'{xc0}-{tag}.gpw', mode='all')

        rpa = RPACorrelation(f'{xc0}-{tag}.gpw',
                             ecut=[200],
                             txt=f'RPAc-{tag}.txt',
                             skip_gamma=True,
                             frequency_scale=2.5)
        e_rpac = rpa.calculate()[0]
        e = e_hf + e_rpac

        if world.rank == 0:
            Path(f'{xc0}-{tag}.gpw').unlink()
            with open(f'RPA-{tag}.result', 'w') as fd:
                print(d, e, e_hf, e_rpac, file=fd)

    return e


if __name__ == '__main__':
    import sys
    xc = sys.argv[1]
    for arg in sys.argv[2:]:
        d = float(arg)
        calculate(xc, d)
```

> ⚠️ 计算代价提醒
> 文中提到每个距离点可能要 ~50 CPU hours（取决于体系/参数/并行），并且需要良好的并行策略（k 点、bands、frequencies）。
> 另外，原文提到 smearing=0.01 eV，但示例代码用的是 `width=0.05`；如果你做精确吸附能曲线，建议对 smearing 也做收敛测试。

---

### 5. Si 的 RPA 内聚能练习（RPA calculation of the cohesive energy of Si）

此部分是一个“从 PBE → EXX@PBE → (RPA+EXX)@PBE”的循序练习。

内聚能写作：
[
E_\text{coh} = E_\text{at} - 0.5 E_\text{bulk}
]

> 提示：因 diamond Si 的 primitive cell 含 2 个 Si 原子，所以 `E_bulk`（每原胞）乘 0.5 得到每原子能量。

文中依次讨论：

* PBE 总能（便宜）
* EXX@PBE（非自洽 exact exchange，贵）
* RPA@PBE 相关能（更贵，需大量空带）

并给出了孤立 Si 原子在不同盒子 L 下的参考能量表：

| L(Å) | (E_\text{PBE}) (eV) |
| ---- | ------------------- |
| 6.0  | -0.664402266578     |
| 7.0  | -0.778484948334     |
| 8.0  | -0.82500272946      |
| 9.0  | -0.841856681349     |
| 10.0 | -0.848092042293     |
| 11.0 | -0.850367362642     |
| 12.0 | -0.85109735188      |

> ✅ 关键提醒
> GPAW 输出的 total energy 不是绝对能，而是相对 reference atoms 的能量；但做**内聚能**时 reference 项通常会相消（文中也提示可在输出里找 `reference = ...` 验证）。

文中还给了一个“孤立原子 PBE+EXX”脚本（很重，不建议作为练习跑），并提供结果表。你可据此计算 EXX@PBE 的内聚能，并与文献值对比。

---

### 6. 基于 TDDFT 的相关能（Correlation energies from TDDFT）

本节强调：RPA 是 TDDFT 相关能表达式在忽略 (f_{xc})（交换-相关核）时的近似。形式上：

[
E_c = -\int_0^{\infty}\frac{d\omega}{2\pi}\int_0^1d\lambda,
\mathrm{Tr}\left[v\chi^\lambda(\omega)-v\chi^{KS}(\omega)\right]
]

响应函数满足：
[
\chi^\lambda(\omega)=\chi^{KS}(\omega)+\chi^{KS}(\omega)\left[\lambda v+f_{xc}^\lambda(\omega)\right]\chi^\lambda(\omega)
]

* **RPA**：忽略 (f_{xc})
* 朴素的 adiabatic kernel 会引发发散/收敛问题
* 通过“renormalization”可得到 **rALDA**、**rAPBE** 等核（引入密度依赖的非局域性）

下面给出若干示例（都很耗资源，适合 HPC）。

---

#### 6.1 示例 1：H 原子的相关能（Correlation energy should be 0）

对单电子体系，精确相关能应为 0。示例流程：

1. 常规 DFT-LDA
2. RPA@LDA
3. rALDA

DFT-LDA 脚本：

```python
from ase import Atoms
from ase.parallel import paropen
from gpaw import GPAW
from gpaw import PW

resultfile = paropen('H.ralda.DFT_corr_energies.txt', 'w')
resultfile.write('DFT Correlation energies for H atom\n')

H = Atoms('H', [(0, 0, 0)])
H.center(vacuum=2.0)
calc = GPAW(mode=PW(400, force_complex_dtype=True),
            parallel={'domain': 1},
            hund=True,
            txt='H.ralda_01_lda.output.txt',
            xc='LDA')

H.calc = calc
E_lda = H.get_potential_energy()
E_c_lda = -calc.get_xc_difference('LDA_X')

resultfile.write(f'LDA correlation: {E_c_lda} eV')
resultfile.write('\n')

calc.diagonalize_full_hamiltonian()
calc.write('H.ralda.lda_wfcs.gpw', mode='all')
```

RPA：

```python
from gpaw.xc.rpa import RPACorrelation

rpa = RPACorrelation('H.ralda.lda_wfcs.gpw',
                     ecut=300,
                     txt='H.ralda_02_rpa_at_lda.output.txt')
rpa.calculate()
```

rALDA：

```python
from gpaw.xc.fxc import FXCCorrelation

fxc = FXCCorrelation('H.ralda.lda_wfcs.gpw',
                     xc='rALDA', txt='H.ralda_03_ralda.output.txt',
                     ecut=300)
fxc.calculate()
```

原文给出的（示例）结果（eV）：

| LDA   | RPA   | rALDA  |
| ----- | ----- | ------ |
| -0.56 | -0.56 | -0.032 |

以及 PBE/rAPBE：

| PBE   | RPA   | rAPBE  |
| ----- | ----- | ------ |
| -0.15 | -0.56 | -0.009 |

---

#### 6.2 示例 2：CO 的解离能（rAPBE 相对 RPA 的改进）

思路与 N2 类似：先做 PBE 得到波函数并存盘，再做 RPA/rAPBE 并外推。这里不再逐句重复，仅强调几个要点：

* `symmetry='off'` 可避免某些非自洽 HF/RPA 计算在对称性约束下出问题
* 需要 `diagonalize_full_hamiltonian()` + `mode='all'`
* 外推通常使用 (E_{cut}^{-1.5}) 标度
* 原文脚本中有 `if 0: # RAPBE not currently working!` 的注释，暗示某些版本/测试环境下 rAPBE 可能不稳定；实际使用请以你所用 GPAW/依赖库版本为准，并做小体系验证

原文总结表（eV）：

| PBE   | HF   | RPA   | rAPBE | Experimental |
| ----- | ---- | ----- | ----- | ------------ |
| 11.74 | 7.37 | 10.64 | 10.88 | 11.23        |

---

#### 6.3 示例 3：金刚石（diamond）的内聚能（RPA vs rAPBE）

要点：RPA 在“体相 vs 孤立原子”的对比里误差抵消较差，容易显著低估内聚能；rAPBE 通过改进核可缓解。

原文总结表（eV）：

| PBE  | HF   | RPA  | rAPBE | Experimental |
| ---- | ---- | ---- | ----- | ------------ |
| 7.75 | 5.17 | 7.06 | 6.65  | 7.55         |

---

#### 6.4 示例 4：不同 kernel 构造方式对 diamond 相关能的影响

kernel 的密度依赖处理存在自由度：

* 默认：two-point **density average**
* 可选：reciprocal-space **wavevector average**（有时更方便、更快）

示例脚本给出 rALDA（两种平均方式）与 RPA 的对比，并在某个未充分收敛的 ecut（131 eV）下给出相关能/电子（eV）：

| rALDA (dens. av.) | rALDA (wave. av) | RPA   |
| ----------------- | ---------------- | ----- |
| -1.16             | -1.13            | -1.40 |

并提到 fully converged RPA 约为 -1.781 eV/electron。

---

## 小结：做能量学计算的“通用检查清单”

为了让本节更可落地，给一个简短 checklist（尤其适合新手避免“算出来但不可信”）：

1. **能量差要统一基准**

   * cohesive / atomization / formation energies 都是差值；确保每个量的 cell、kpts、ecut/h、smearing、XC、setup 一致或可比。
2. **按“每原子/每原胞”明确归一**

   * cohesive energy 这类最容易因晶胞原子数搞错一个系数。
3. **孤立原子是高风险点**

   * 大真空、合理自旋/ Hund、必要时 direct minimization、必要时关对称性/混合。
4. **带电缺陷必须做电静力校正与一致性检查**

   * (\Delta V(z)) 是否平坦、epsilon 是否合理、模型电荷中心是否对齐。
5. **RPA/FXC 类计算要重视收敛与外推**

   * ecut、bands、kpts、frequency sampling 都要有收敛策略；能量绝对值更难收敛，能量差往往更快。

如果你希望我把这份 **Energetics** 中文版进一步“工程化”（例如：把每个示例补上命令行运行方式、输出文件关键字段解读、常见报错/内存估算建议、以及推荐的收敛顺序），我也可以在不改变原始代码块的前提下，追加一个“实战指南”附录。

## 电子结构（Electronic structure）

本节介绍在 GPAW/ASE 工作流中常用的电子结构分析：**态密度（DOS/PDOS/LDOS）**、**能带结构（band structure）**、**自旋轨道耦合（SOC）**、**非自洽 HSE06 能带**、**超胞能带展开（band unfolding）**、**GLLB‑sc 带隙**、**G0W0 带隙**以及**PBE0（非自洽）**相关示例。
我会在保留关键英文术语/接口名的同时，补充一些实用注意事项（尤其是能量参考、展宽、k 点与投影方法的适用范围）。

---

### 1. 态密度（Density of states, DOS）

GPAW 的 calculator 支持多种**投影态密度（projected DOS, PDOS）/局域态密度（local DOS, LDOS）**的获取方式，常见有四类：

* **Total DOS**：总态密度（全体系）。
* **Molecular orbital PDOS**：投影到分子气相轨道（molecular orbitals）的 PDOS（常用于吸附体系：看“分子轨道在表面上的谱分布”）。
* **Atomic orbital PDOS**：投影到原子轨道/投影子（projector）的 PDOS（例如 d-band center 分析）。
* **Wigner–Seitz LDOS**：基于 Wigner–Seitz 几何分区的局域 DOS（目前有一些限制）。

下面分别说明。

---

#### 1.1 总 DOS（Total DOS）

总态密度可以用 GPAW 的方法：

* `calc.get_dos(spin=0, npts=201, width=None)`

其中：

* `spin`：自旋通道（非自旋极化体系通常用 `spin=0`）。
* `npts`：能量网格点数（越大曲线越平滑，但文件/计算也更大）。
* `width`：Gaussian 展宽（单位 eV，常用 0.05–0.2 eV；金属体系通常需要合适的展宽来平滑尖峰，但也要避免“抹掉”精细结构）。

> 实用建议
>
> * **DOS 的可靠性高度依赖 k 点采样**：体相金属需要较密 k 网格，否则 DOS 会“锯齿/假峰”。
> * `e - e_f`（相对费米能级）是最常用的横轴表示方式。
> * 如果你看到 DOS 在费米能附近不稳定、随参数变化很大，优先检查：`kpts`、smearing（如 `FermiDirac(width=...)`）、以及是否真正收敛。

原文建议你查看 `dos.py` 脚本，并用它为 `bulk fcc aluminium` 作图（横轴为相对费米能级的能量）。

---

#### 1.2 分子轨道 PDOS（Molecular Orbital PDOS）

下面示例计算 CO 吸附在 Pt(111) 表面上的总 DOS，并将 DOS 投影到 **CO 气相轨道（gas phase orbitals）**上。
（`top.py` 用于生成所需的 `.gpw` 文件：如 `top.gpw`、`CO.gpw`。）

**PDOS 脚本（`pdos.py`）**：

```python
# web-page: pdos.png
from gpaw import GPAW, restart
import matplotlib.pyplot as plt

# Density of States
plt.subplot(211)
slab, calc = restart('top.gpw')
e, dos = calc.get_dos(spin=0, npts=2001, width=0.2)
e_f = calc.get_fermi_level()
plt.plot(e - e_f, dos)
plt.axis([-15, 10, None, 4])
plt.ylabel('DOS')

molecule = range(len(slab))[-2:]

plt.subplot(212)
c_mol = GPAW('CO.gpw')
for n in range(2, 7):
    print('Band', n)
    # PDOS on the band n
    wf_k = [kpt.psit_nG[n] for kpt in c_mol.wfs.kpt_u]
    P_aui = [[kpt.P_ani[a][n] for kpt in c_mol.wfs.kpt_u]
             for a in range(len(molecule))]
    e, dos = calc.get_all_electron_ldos(mol=molecule, spin=0, npts=2001,
                                        width=0.2, wf_k=wf_k, P_aui=P_aui)
    plt.plot(e - e_f, dos, label='Band: ' + str(n))
plt.legend()
plt.axis([-15, 10, None, None])
plt.xlabel('Energy [eV]')
plt.ylabel('All-Electron PDOS')
plt.savefig('pdos.png')
plt.show()
```

运行该脚本时，会对每个自旋与 k 点打印 (\int d\varepsilon \rho_i(\varepsilon))。
如果分子轨道 (\psi_i(\mathbf r)) 能够被 slab 的 Kohn–Sham 轨道较好展开，该积分应接近 1；因此它也可看作一种**“KS 空间完备性/可展开性”**的度量。

* 原文指出：band 7 和 8 较为离域（delocalized），无法被 slab 本征态良好展开。你可把 `range(2,7)` 改为 `range(2,9)`，会发现积分显著小于 1。

`calc.get_all_electron_ldos()` 的核心思想是：计算重叠（overlap）的模平方，并用指定宽度的归一化 Gaussian 做能量展宽。
注意：

* 能量单位是 eV，且相对于**平均势（average potential）**，因此示例中通常再减去 `e_f` 来以费米能为零点。
* 如果设置 `raw=True`，函数将返回**重叠和能量（Hartree）**，而不进行展宽与单位转换。

> 实用建议
> `.gpw`（带波函数）可能很大；如果你只想反复画图/换展宽参数，建议把原始 overlap 数据保存起来（例如 pickle），避免重复加载/处理庞大的波函数文件。

下面是“将 overlaps pickle 保存”的示例。

**Pickle 脚本（`p1.py`）**：

```python
from gpaw import GPAW, restart
import pickle

slab, calc = restart('top.gpw')
c_mol = GPAW('CO.gpw')
molecule = range(len(slab))[-2:]
e_n = []
P_n = []
for n in range(c_mol.get_number_of_bands()):
    print('Band: ', n)
    wf_k = [kpt.psit_nG[n] for kpt in c_mol.wfs.kpt_u]
    P_aui = [[kpt.P_ani[a][n] for kpt in c_mol.wfs.kpt_u]
             for a in range(len(molecule))]
    e, P = calc.get_all_electron_ldos(mol=molecule, wf_k=wf_k, spin=0,
                                      P_aui=P_aui, raw=True)
    e_n.append(e)
    P_n.append(P)
pickle.dump((e_n, P_n), open('top.pickle', 'wb'))
```

**从 pickle 画 PDOS（`p2.py`）**：

```python
from ase.units import Hartree
from gpaw import GPAW
from gpaw.utilities.dos import fold
import pickle
import matplotlib.pyplot as plt

e_f = GPAW('top.gpw').get_fermi_level()

e_n, P_n = pickle.load(open('top.pickle', 'rb'))
for n in range(2, 7):
    e, ldos = fold(e_n[n] * Hartree, P_n[n], npts=2001, width=0.2)
    plt.plot(e - e_f, ldos, label='Band: ' + str(n))
plt.legend()
plt.axis([-15, 10, None, None])
plt.xlabel('Energy [eV]')
plt.ylabel('PDOS')
plt.show()
```

---

#### 1.3 原子轨道/投影子 PDOS（Atomic Orbital PDOS）

该类 PDOS 可由 calculator 方法获取：

* `calc.get_orbital_ldos(a, spin=0, angular='spdf', npts=201, width=None)`

参数说明：

* `a`：原子索引（Atoms 对象中的 atom index）。
* `angular`：

  * 若为整数：指定该原子某个特定 projector（投影子）函数编号；
  * 若为字符串：由 `'s' 'p' 'd' 'f'` 组成（可以多个，如 `'pd'`），表示对该原子所有对应角动量的束缚态 projectors 求和。

如何知道整数 `angular` 对应哪个 projector？可用：

```python
from gpaw.utilities.dos import print_projectors
print_projectors('Au')
```

> 重要说明（非常容易误解）
> 原子 partial waves（部分波）**不是正交归一基（non-orthonormal basis）**，因此这里的“PDOS 投影”不具备严格正交投影的数学性质。
> 但是它仍然非常有用，通常可作为 DOS 的局域/轨道成分的**定性**指标（例如判断某能区主要是 d 态还是 sp 态）。

下面给出一个常见用法：计算 Au 的 d 态 PDOS，并进一步用积分定义 **d-band center** 与 **d-band width**（表面催化中非常常见的描述符）。

**生成 `au.gpw`（`atomic_orbital_gs.py`）**：

```python
from ase.build import bulk
from gpaw import GPAW

atoms = bulk('Au')
k = 8
atoms.calc = GPAW(mode='pw',
                  kpts=(k, k, k))
atoms.get_potential_energy()
atoms.calc.write('au.gpw')
```

**计算并画 d-PDOS 与 d-band center（`atomic_orbital_pdos.py`）**：

```python
import numpy as np
from scipy.integrate import trapezoid
import matplotlib.pyplot as plt
from gpaw import GPAW

calc = GPAW('au.gpw')
energy, pdos = calc.get_orbital_ldos(a=0, angular='d')
energy -= calc.get_fermi_level()
I = trapezoid(pdos, energy)
center = trapezoid(pdos * energy, energy) / I
width = np.sqrt(trapezoid(pdos * (energy - center)**2, energy) / I)
plt.plot(energy, pdos)
plt.xlabel('Energy (eV)')
plt.ylabel('d-projected DOS on atom 0')
plt.title(f'd-band center = {center:.1f} eV, d-band width = {width:.1f} eV')
# plt.show()
plt.savefig('au-ddos.png')
```

> ⚠️ Warning（强烈建议遵守）
> 你应该**先看 PDOS 曲线**，确认它在高能区没有不合理的“上翘尾巴（tail）”，再相信计算出的 center/width。
> 由于投影子很局域，有时会导致高能区出现人工尾巴，从而显著拉偏积分结果。
> 若出现这种情况，可考虑把 DOS 投影到 **LCAO orbitals**（局域原子轨道基）上，因为它们通常更“宽”，尾巴问题可能缓解。
> 但注意：LCAO 投影通常需要额外计算（grid 模式默认不会预先给出 LCAO 投影），而 projector 投影是“现成的”，几乎不额外耗时。

---

#### 1.4 Wigner–Seitz 局域 DOS（Wigner–Seitz LDOS）

> 注：目前该功能仅对 **Gamma 点**计算实现（即不使用 k 点采样）。

可通过：

* `calc.get_wigner_seitz_ldos(a, spin=0, npts=201, width=None)`

它代表在原子 (a) 处的局域 DOS 探针（基于 Wigner–Seitz 几何分区）。一些性质：

* 对所有原子求和可恢复总 DOS（但 `calc.get_dos` 更高效）。
* 对能量积分给出该分区中的电子数（但 `calc.get_wigner_seitz_densities(spin)` 更高效）。
* 分区完全基于**几何准则**（与电荷密度拓扑无关）。

若你需要更“物理化”的电荷划分，可用 `Bader` 算法（通常更准确但更耗时）；Wigner–Seitz 更快。

---

#### 1.5 投影到 LCAO 轨道的 PDOS（PDOS on LCAO orbitals）

DOS 也可以投影到 LCAO 的基函数上。计算时需要给出一个“轨道子空间”（subspace），本质是你要投影的 AO（原子轨道）索引集合。

例如：某个原子的 p 轨道在 AO 列表中对应索引 41、42、43，则想要该 p 子空间 PDOS，就传入 `[41, 42, 43]`。

下面给出示例脚本（对应 `lcaodos_gs.py` 与 `lcaodos_plt.py`）。
（该示例使用 `RestartLCAODOS` 从 `.gpw` 重启文件中提取并计算 subspace PDOS。）

```python
# web-page: lcaodos.png
import matplotlib.pyplot as plt
import numpy as np
from ase.io import read
from ase.units import Hartree

from gpaw import GPAW
from gpaw.utilities.dos import RestartLCAODOS, fold


name = 'HfS2'
calc = GPAW(name + '.gpw', txt=None)
atoms = read(name + '.gpw')
ef = calc.get_fermi_level()

dos = RestartLCAODOS(calc)
energies, weights = dos.get_subspace_pdos(range(51))
e, w = fold(energies * Hartree, weights, 2000, 0.1)

e, m_s_pdos = dos.get_subspace_pdos([0, 1])
e, m_s_pdos = fold(e * Hartree, m_s_pdos, 2000, 0.1)
e, m_p_pdos = dos.get_subspace_pdos([2, 3, 4])
e, m_p_pdos = fold(e * Hartree, m_p_pdos, 2000, 0.1)
e, m_d_pdos = dos.get_subspace_pdos([5, 6, 7, 8, 9])
e, m_d_pdos = fold(e * Hartree, m_d_pdos, 2000, 0.1)

e, x_s_pdos = dos.get_subspace_pdos([25])
e, x_s_pdos = fold(e * Hartree, x_s_pdos, 2000, 0.1)
e, x_p_pdos = dos.get_subspace_pdos([26, 27, 28])
e, x_p_pdos = fold(e * Hartree, x_p_pdos, 2000, 0.1)

w_max = []
for i in range(len(e)):
    if (-4.5 <= e[i] - ef <= 4.5):
        w_max.append(w[i])

w_max = np.asarray(w_max)
```

关于该脚本，原文解释要点如下（这里略作扩展）：

* 该计算中共有 51 个 basis functions，把 `range(51)` 全部拿来投影即可得到“总 DOS”（在 LCAO 轨道空间中的总投影）。
* 随后分别计算金属原子的 s/p/d 子空间 PDOS，以及硫族原子（chalcogen）s/p 子空间 PDOS。
* 在选择子空间时，这里“localized part”没有纳入，只选取“confined orbital part（较大 rc）”那部分轨道（这属于基组构造细节；具体实现与基组文件相关）。

> 实用补充：如何“自动获取某原子某角动量轨道的索引”？
> 示例里索引是硬编码的（如 `[2,3,4]`）。更稳健的做法通常是：从 LCAO basis 描述中解析每个 AO 的原子归属与 (n,l,m) 标签，再自动筛选。原文提到“更聪明的自动化方式稍后会介绍”。

最后是绘图部分（`lcaodos_plt.py` 的尾段）：

```python
plt.plot(e - ef, w, label='Total', c='k', lw=2, alpha=0.7)
plt.plot(e - ef, x_s_pdos, label='X-s', c='g', lw=2, alpha=0.7)
plt.plot(e - ef, x_p_pdos, label='X-p', c='b', lw=2, alpha=0.7)
plt.plot(e - ef, m_s_pdos, label='M-s', c='y', lw=2, alpha=0.7)
plt.plot(e - ef, m_p_pdos, label='M-p', c='c', lw=2, alpha=0.7)
plt.plot(e - ef, m_d_pdos, label='M-d', c='r', lw=2, alpha=0.7)

plt.axis(ymin=0., ymax=np.max(w_max), xmin=-4.5, xmax=4.5, )
plt.xlabel(r'$\epsilon - \epsilon_F \ \rm{(eV)}$')
plt.ylabel('DOS')
plt.legend(loc=1)
plt.savefig('lcaodos.png')
plt.show()
```

---

### 2. 电子能带结构计算（Calculation of electronic band structures）

本教程以 Si 为例，计算沿布里渊区高对称路径的 band structure。

#### 2.1 第一步：基态计算并保存 `.gpw`

对小体相体系，使用 **plane-wave mode（PW）**通常更合适。

```python
from ase.build import bulk
from gpaw import GPAW, PW, FermiDirac

# Perform standard ground state calculation (with plane wave basis)
si = bulk('Si', 'diamond', 5.43)
calc = GPAW(mode=PW(200),
            xc='PBE',
            kpts=(8, 8, 8),
            random=True,  # random guess (needed if many empty bands required)
            occupations=FermiDirac(0.01),
            txt='Si_gs.txt')
si.calc = calc
si.get_potential_energy()
ef = calc.get_fermi_level()
calc.write('Si_gs.gpw')
```

> 实用建议
>
> * `random=True` 通常用于需要较多空带时的初始猜测（避免某些情况下迭代本征求解器卡住）。
> * `FermiDirac(0.01)` 是较小 smearing；对半导体/绝缘体一般没问题，但仍需收敛测试。

#### 2.2 第二步：沿高对称路径计算本征值（fixed density）

使用 `fixed_density`：把密度/势固定在基态结果上，只计算指定 k 路径上的本征值（相当于非自洽能带）。

示例路径：`kpts={'path': 'GXWKL', 'npoints': 60}`
（FCC 晶格的特殊点定义参见 `ase.dft.kpoints.special_points`。）

同时为了计算所有路径点，关闭对称性：`symmetry='off'`。

```python
# Restart from ground state and fix potential:
calc = GPAW('Si_gs.gpw').fixed_density(
    nbands=16,
    symmetry='off',
    kpts={'path': 'GXWKL', 'npoints': 60},
    convergence={'bands': 8})
```

> 参数提示
>
> * `nbands`：要输出多少条能带（包含空带）。
> * `convergence={'bands': 8}`：固定密度下的本征值迭代收敛设置（这里含义是至少保证一定数量 bands 收敛；具体含义随版本/实现可能略有差异，但实践中常用于 bandstructure 的稳定收敛）。

#### 2.3 第三步：绘制能带

使用 ASE 的 `BandStructure` 工具：

```python
bs = calc.band_structure()
bs.plot(filename='bandstructure.png', show=True, emax=10.0)
```

完整脚本对应 `bandstructure.py`。

---

#### 2.4 自旋轨道耦合（Spin-orbit coupling, SOC）的影响

下面示例放大 Si 的 VBM 附近，看 SOC 引起的细小分裂。关键点：

* 需要**非共线（non-collinear）**计算。
* 通过 `experimental={'soc': True, ...}` 开 SOC（这是 GPAW 的实验接口/参数集）。
* 关闭对称性：`symmetry='off'`。

```python
from ase.build import bulk
from gpaw import GPAW, PW, FermiDirac
import numpy as np

# Non-collinear ground state calculation:
si = bulk('Si', 'diamond', 5.43)
si.calc = GPAW(mode=PW(400),
               xc='LDA',
               experimental={'magmoms': np.zeros((2, 3)),
                             'soc': True},
               kpts=(8, 8, 8),
               symmetry='off',
               occupations=FermiDirac(0.01))
si.get_potential_energy()

bp = si.cell.bandpath('LGX', npoints=100)
bp.plot()

# Restart from ground state and fix density:
calc2 = si.calc.fixed_density(
    nbands=16,
    basis='dzp',
    symmetry='off',
    kpts=bp,
    convergence={'bands': 8})

bs = calc2.band_structure()
bs = bs.subtract_reference()

# Zoom in on VBM:
bs.plot(filename='si-soc-bs.png', show=True, emin=-1.0, emax=0.5)
```

> 实用提醒
>
> * SOC 对轻元素（如 Si）效应较小，但对重元素（Au、Bi、W 等）以及拓扑材料非常关键。
> * SOC+非共线计算对并行/收敛更敏感，建议先用较小体系和较紧的 SCF 选项做 sanity check。

---

#### 2.5 非自洽 HSE06（Non self-consistent HSE06）

这一部分展示如何在一个半局域 DFT（示例为 LDA）基态之上，**非自洽**地计算 HSE06 的本征值（类似 “HSE06@LDA” 的概念）。

> 术语说明
>
> * **Non self-consistent HSE06**：不迭代更新密度/势，只在已有的 DFT 轨道/密度上评估 HSE06 的本征值修正。
> * 这种方法成本比全自洽 hybrid 小，但结果仍需要验证（尤其是能带排序/费米能附近态）。

**计算脚本（生成 `bs.pckl`）**：

```python
import pickle
from pathlib import Path

from ase.build import mx2
from gpaw.mpi import world
from gpaw.new.ase_interface import GPAW
from gpaw.new.pw.nschse import NonSelfConsistentHSE06


def mos2():
    """Do LDA calculation for MoS2 layer."""
    atoms = mx2(formula='MoS2', kind='2H', a=3.184, thickness=3.13,
                size=(1, 1, 1))
    atoms.center(vacuum=3.5, axis=2)
    k = 6
    atoms.calc = GPAW(mode={'name': 'pw', 'ecut': 400},
                      kpts=(k, k, 1),
                      txt='lda.txt')
    atoms.get_potential_energy()
    return atoms


def bandstructure(gs_calc, bp):
    """Calculate HSE06 bandstructure on top of LDA."""
    fermi_level = gs_calc.get_fermi_level()
    vacuum_level = gs_calc.dft.vacuum_level()
    N = 13 + 4  # 13 occupied + 4 empty
    bs_calc = gs_calc.fixed_density(
        kpts=bp,
        convergence={'bands': N},
        symmetry='off',
        txt='gmkg.txt')
    lda_skn = bs_calc.eigenvalues()
    hse = NonSelfConsistentHSE06.from_dft_calculation(
        gs_calc.dft, 'hse06.txt')
    hse_skn = hse.calculate(bs_calc.dft.ibzwfs, na=0, nb=N)
    # Return energies relative to vacuum level:
    return (lda_skn[0, :, :N] - vacuum_level,
            hse_skn[0] - vacuum_level,
            fermi_level - vacuum_level)


def run():
    atoms = mos2()
    bp = atoms.cell.bandpath('GMKG', npoints=50)
    lda_kn, hse_kn, fermi_level = bandstructure(atoms.calc, bp)
    if world.rank == 0:
        Path('bs.pckl').write_bytes(
            pickle.dumps((bp, lda_kn, hse_kn, fermi_level)))


if __name__ == '__main__':
    run()
```

**绘图脚本（`hse06.png`）**：

```python
# wep-page: hse06.png
import pickle
from pathlib import Path
import matplotlib.pyplot as plt


def plot(bp, lda_kn, hse_kn, fermi_level):
    ax = plt.subplot()
    x, xlabels, labels = bp.get_linear_kpoint_axis()
    labels = [label.replace('G', r'$\Gamma$') for label in labels]
    label = 'LDA'
    for y in lda_kn.T:
        ax.plot(x, y, color='C0', label=label)
        label = None
    label = 'HSE06@LDA'
    for y in hse_kn.T:
        ax.plot(x, y, color='C1', label=label)
        label = None
    ax.hlines(fermi_level, 0.0, x[-1], colors='black', label='Fermi-level')
    ax.legend()
    ax.set_xlim(0.0, x[-1])
    ax.set_ylim(fermi_level - 3.0, fermi_level + 3.0)
    ax.set_xticks(xlabels)
    ax.set_xticklabels(labels)
    ax.set_ylabel('eigenvalues relative to vacuum-level [eV]')
    # plt.show()
    plt.savefig('hse06.png')


if __name__ == '__main__':
    path = Path('bs.pckl')
    plot(*pickle.loads(path.read_bytes()))
```

该模块核心类与函数接口（保持原文）：

**class `gpaw.new.pw.nschse.NonSelfConsistentHSE06(ibzwfs, density, pot_calc, setups, relpos_ac, log='-')`**

* **calculate(ibzwfs, na=0, nb=0)**

  * 在多个 k 点上计算本征值。
  * 返回 DFT 与 HSE06 的本征值（单位 eV）。

* **calculate_one_kpt(psit2_nG, P2_ani, spin)**

  * 在单个 k 点计算本征值（单位 eV）。

* **classmethodfrom_dft_calculation(dft, log='-')**

  * 从已完成的 DFT 计算创建 HSE06 本征值计算器。

以及：

**`gpaw.new.pw.hybrids.truncated_coulomb(pw, omega=0.11, yukawa=False)`**
截断库仑（truncated Coulomb）的傅里叶变换。

Real space（实空间）：

$$\frac{erfc(ωr)}{r} \cdot$$

Reciprocal space（倒空间）：

$$\frac{4π}{(\mathbf{G}+\mathbf{k})^{2}} (1 - exp(-(\mathbf{G}+\mathbf{k})^{2} /(4 ω^{2} )))$$

（当 (\mathbf{G}+\mathbf{k}=0) 时极限为 (π/ω^2)。）

> 补充说明
> `erfc` 为 complementary error function（互补误差函数）。`omega` 是 screening 参数（HSE 系列常用值约 0.11 bohr⁻¹ 对应 HSE06）。

---

### 3. 超胞能带展开（Electronic Band Structure Unfolding for Supercell Calculations）

#### 3.1 简要理论背景（Brief theory overview）

超胞（supercell, SC）能带由于晶胞变大、能带折叠（folding）会变得“很乱”。但我们做超胞往往是为了研究：

* 缺陷（defects）
* 畸变（distortions）
* 合金/掺杂/吸附等轻微破坏周期性的问题

此时常希望回答：**缺陷引入后，原本原胞（primitive cell, PC）的能带特征在多大程度上仍可辨认？**
带展开（unfolding）就是把 SC 的能带信息映射/投影回 PC 的布里渊区，从而更直观地识别：

* 仍保留的“本征能带”部分（高权重）
* 缺陷态/局域态（往往表现为带隙内的谱权重）

---

#### 3.2 示例：MoS₂ 3×3 超胞 + 单个 S 空位

##### (1) Ground state：3×3 超胞含空位的基态计算

超胞较大时，使用 `mode='lcao'` 与 `basis='dzp'` 往往更省。

```python
from ase.build import mx2
from gpaw import GPAW, FermiDirac

structure = mx2(formula='MoS2', kind='2H', a=3.184, thickness=3.127,
                size=(3, 3, 1), vacuum=7.5)
structure.pbc = (1, 1, 1)

# Create vacancy
del structure[2]

calc = GPAW(mode='lcao',
            basis='dzp',
            xc='LDA',
            kpts=(4, 4, 1),
            occupations=FermiDirac(0.01),
            txt='gs_3x3_defect.txt')

structure.calc = calc
structure.get_potential_energy()
calc.write('gs_3x3_defect.gpw', 'all')
```

最后一行写出包含波函数的 `.gpw`，用于后续展开分析。

---

##### (2) 定义 PC 布里渊区路径，并找到对应 SC 的 (\vec{K})

核心：给定超胞变换矩阵 (M)，把 PC 的 (\vec{k}) 映射到 SC 的 (\vec{K})。

```python
from ase.build import mx2

from gpaw import GPAW
from gpaw.unfold import Unfold, find_K_from_k

a = 3.184
PC = mx2(a=a).get_cell(complete=True)
bp = PC.get_bravais_lattice().bandpath('MKG', npoints=48)
x, X, _ = bp.get_linear_kpoint_axis()

M = [[3, 0, 0], [0, 3, 0], [0, 0, 1]]

Kpts = []
for k in bp.kpts:
    K = find_K_from_k(k, M)[0]
    Kpts.append(K)
```

---

##### (3) 在这些 (\vec{K}) 点做非自洽 bands 计算（fixed density）

```python
calc_bands = GPAW('gs_3x3_defect.gpw').fixed_density(
    kpts=Kpts,
    symmetry='off',
    nbands=220,
    convergence={'bands': 200})

calc_bands.write('bands_3x3_defect.gpw', 'all')
```

---

##### (4) 展开并计算谱函数（spectral function）

先创建 `Unfold` 对象：

```python
unfold = Unfold(name='3x3_defect',
                calc='bands_3x3_defect.gpw',
                M=M,
                spinorbit=False)
```

计算 spectral function：

```python
unfold.spectral_function(kpts=bp.kpts, x=x, X=X,
                         points_name=['M', 'K', 'G'])
```

> 并行建议
> 该步骤可对 k 点并行（supercell 较大时很有必要）。

注意：该函数会输出两个 pickle：

* `weights_3x3_defects.pckl`：包含本征值 (\epsilon_{\vec{K}m}) 与权重 (P_{\vec{K}m}(\vec{k}))
* `sf_3x3_defects.pckl`：包含谱函数与能量数组

---

##### (5) 绘制谱函数（Spectral Function）

```python
# web-page: sf_3x3_defect_spec.png
from gpaw import GPAW
from gpaw.unfold import plot_spectral_function

calc = GPAW('gs_3x3_defect.gpw', txt=None)
ef = calc.get_fermi_level()

plot_spectral_function(filename='sf_3x3_defect',
                       color='blue',
                       eref=ef,
                       emin=-3,
                       emax=3)
```

会生成谱函数图，你能清楚看到带隙内的缺陷态（defect states）。

接口说明（保留原文）：

`**gpaw.unfold.plot_spectral_function(filename, color='blue', eref=None, emin=None, emax=None, scale=1)**`
用于绘制沿 k 路径的谱函数。

**class `gpaw.unfold.Unfold(name=None, calc=None, M=None, spin=0, spinorbit=None, theta=90, scale=1.0, phi=90, world=None)`**
用于将 SC 的 bands 展开到 PC。约定：大写变量常对应 SC，小写对应 PC。

---

### 4. 使用 GLLB‑sc 计算带隙（Calculating band gap using the GLLB-sc functional）

GLLB‑sc（`GLLBSC`）常用于半导体带隙计算：通过响应势（response potential）得到**导数不连续（derivative discontinuity）**，把 Kohn–Sham 带隙修正为更接近本征带隙（fundamental gap）。

#### 4.1 基态计算（以 Si 为例）

```python
from ase.build import bulk
from gpaw import GPAW, PW, FermiDirac

# Ground state calculation
atoms = bulk('Si', 'diamond', 5.431)
calc = GPAW(mode=PW(200),
            xc='GLLBSC',
            kpts=(7, 7, 7),  # Choose and converge carefully!
            occupations=FermiDirac(0.01),
            txt='gs.out')
atoms.calc = calc
atoms.get_potential_energy()
```

#### 4.2 计算不连续势与带隙

```python
# Calculate the discontinuity potential and the discontinuity
homo, lumo = calc.get_homo_lumo()
response = calc.hamiltonian.xc.response
dxc_pot = response.calculate_discontinuity_potential(homo, lumo)
KS_gap, dxc = response.calculate_discontinuity(dxc_pot)

# Fundamental band gap = Kohn-Sham band gap + derivative discontinuity
QP_gap = KS_gap + dxc

print(f'Kohn-Sham band gap:         {KS_gap:.2f} eV')
print(f'Discontinuity from GLLB-sc: {dxc:.2f} eV')
print(f'Fundamental band gap:       {QP_gap:.2f} eV')
```

> ⚠️ Warning（关键！）
> 计算不连续势需要 HOMO/VBM 与 LUMO/CBM 作为输入。要得到可靠结果，k 点采样必须包含真实的 VBM 与 CBM（否则你会拿到“错的边缘态”，从而带隙修正也错）。

原文示例：Si 得到约 1.05 eV（实验约 1.17 eV）。

---

#### 4.3 用 band structure calculator 精确定位 VBM/CBM（推荐策略）

一种稳妥做法是：用一个能带路径计算器 `bs_calc` 确保路径上包含 VBM/CBM，再用它返回更准确的 band edges。

```python
from ase.build import bulk
from gpaw import GPAW, PW, FermiDirac

# Ground state calculation
atoms = bulk('Si', 'diamond', 5.431)
calc = GPAW(mode=PW(200),
            xc='GLLBSC',
            kpts={'size': (8, 8, 8), 'gamma': True},
            occupations=FermiDirac(0.01),
            txt='gs.out')
atoms.calc = calc
atoms.get_potential_energy()

# Band structure calculation with fixed density
bs_calc = calc.fixed_density(nbands=10,
                             kpts={'path': 'LGXWKL', 'npoints': 60},
                             symmetry='off',
                             convergence={'bands': 8},
                             txt='bs.out')

# Plot the band structure
bs = bs_calc.band_structure().subtract_reference()
bs.plot(filename='bs_si.png', emin=-6, emax=6)
```

然后：

* 用 `bs_calc` 获取更准确的 `homo, lumo`
* 用基态 `calc` 计算 discontinuity potential
* 用 `bs_calc` 计算 discontinuity（确保包含 VBM/CBM）

```python
# Get the accurate HOMO and LUMO from the band structure calculator
homo, lumo = bs_calc.get_homo_lumo()

# Calculate the discontinuity potential using the ground state calculator and
# the accurate HOMO and LUMO
response = calc.hamiltonian.xc.response
dxc_pot = response.calculate_discontinuity_potential(homo, lumo)

# Calculate the discontinuity using the band structure calculator
bs_response = bs_calc.hamiltonian.xc.response
KS_gap, dxc = bs_response.calculate_discontinuity(dxc_pot)

# Fundamental band gap = Kohn-Sham band gap + derivative discontinuity
QP_gap = KS_gap + dxc

print(f'Kohn-Sham band gap:         {KS_gap:.2f} eV')
print(f'Discontinuity from GLLB-sc: {dxc:.2f} eV')
print(f'Fundamental band gap:       {QP_gap:.2f} eV')
```

> ⚠️ Warning
> 即使使用 band structure calculator，也要确保 bandpath 上确实包含 VBM 和 CBM。

> Note
> 原文指出：简单方法中如果把 `kpts=(7,7,7)` 改为 `kpts={'size': (8,8,8), 'gamma': True}` 可能导致结果差异很大。使用“第二种策略”（band structure 来定位边缘）通常能减少这种 k 网格差异带来的问题。

完整脚本：`gllbsc_si_simple.py` 与 `gllbsc_si_band_edges.py`。

#### 4.4 自旋极化 GLLB‑SC（Spin-polarized GLLB-SC）

原文说明：自旋极化 GLLB‑SC 目前在 svn trunk 中实现，但存在与 fermi smearing 和最高轨道参考能相关的收敛问题，且部分未充分测试。如要使用，建议先联系 `mikael.kuisma@tut.fi`。

---

### 5. G0W0 计算 Si 带隙（G0W0 calculation of the band gap of silicon）

> 背景理解（非常重要）
> DFT 的 Kohn–Sham 能带并不严格等同于电子加/去能（quasiparticle energies）。
> 带隙往往被 DFT（尤其是 LDA/GGA）低估，因为它没有完整包含加/去电子时的屏蔽与多体效应。
> GW 近似通过自能（self-energy）修正得到更接近准粒子谱（quasiparticle spectrum）的能量。

---

#### 5.1 基态计算与空带准备（必须使用平面波）

GW 在 GPAW 的实现需要 plane waves（PW），并且需要大量空带。常用流程：

1. 先做常规基态（收敛密度）
2. 再对最终哈密顿量做**完全对角化**得到大量空带：`diagonalize_full_hamiltonian()`
3. 保存包含波函数的 `.gpw`：`write(..., 'all')`

示意代码：

```python
calc = GPAW(mode=PW(300), ...) # Note we MUST use plane waves

atoms.calc = calc
atoms.get_potential_energy()   # Calculate ground state

calc.diagonalize_full_hamiltonian() # Does what it says ;)
calc.write('my_file.gpw', 'all')    # Saves calculator AND wavefunctions
```

注意：

* 要使用 `diagonalize_full_hamiltonian()` 必须 PW 模式：`mode=PW()` 并 `from gpaw import PW`。
* 目前 GPAW GW 实现不支持磁性体系：确保 `spinpol=False`。

原文练习建议：用 `ecut=200 eV`、`kpts={'size': (4,4,4), 'gamma': True}` 做一套基态并保存。
为何要 gamma-centered？因为很多半导体（如 Si）VBM/CBM 位置与 Γ 点相关，且对称性/采样也更容易控制。

---

#### 5.2 运行 G0W0（示例脚本）

```python
import numpy as np
from ase.parallel import parprint
from gpaw.response.g0w0 import G0W0

# We start by setting up a G0W0 calculator object
gw = G0W0('Si_gs.gpw',             # Path to groundstate gpw file
          filename='Si_g0w0_ppa',  # filename base for output files
          kpts=None,               # List of quasiparticle k-point indices
                                   # or None = all k-points
          bands=(3, 5),            # Range of quasiparticle bands - last
                                   # index is NOT included
          ecut=100.,               # Plane wave basis cut-off energy
          ppa=True)                # Use Plasmon Pole Approximation

# Perform the GW calculation. The results, ie. quasiparticle energies, as
# well as original Kohn-Sham eigenvalues, occupation numbers, DFT XC and
# self-energy contributions and renormalization factors are returned as a
# python dictionary object
result = gw.calculate()

ks_skn = result['eps']               # Get Kohn-Sham eigenvalues
ks_cbmin = np.amin(ks_skn[0, :, 1])  # DFT conduction band minimum
ks_vbmax = np.amax(ks_skn[0, :, 0])  # DFT valence band maximum
ks_gap = ks_cbmin - ks_vbmax         # DFT band gap

qp_skn = result['qp']                # GW quasiparticle energies
qp_cbmin = np.amin(qp_skn[0, :, 1])  # GW conduction band minimum
qp_vbmax = np.amax(qp_skn[0, :, 0])  # GW valence band maximum
qp_gap = qp_cbmin - qp_vbmax         # GW band gap

parprint(f'Kohn-Sham gap = {ks_gap:.3f}')
parprint(f'G0W0 gap = {qp_gap:.3f}')
```

* `ppa=True`：Plasmon Pole Approximation（PPA），用一个简单模型拟合介电函数频率依赖，速度快。
* `ppa=False`：在频率网格上显式计算介电函数，通常更慢。

---

#### 5.3 收敛性（Convergence）要点

GW 中关键收敛参数通常包括：

* **k 点采样**
* **频率网格密度**（`domega0`，默认 0.025；减半会显著增加计算量）
* **用于响应函数/屏蔽势的平面波截断**（`ecut`）
* **空带数 nbands**（通常与平面波数目同量级）

原文提示：

* GPAW 默认选择的 bands 数接近 plane waves 数，通常不必更多。直观解释：响应函数构造是对一组基（如 PW）上的态求和，超过该基维度的态并不能被更精确表示，继续加 bands 收益很小但代价巨大。
* GW 计算在 plane waves 数上的代价近似 ~(N_\text{PW}^3)（因为还会带上 bands 的增长），因此 ecut 提升很昂贵。
* 先收敛 DFT gap（k 点），再收敛 GW gap（k 点、ecut、频率网格等）是比较稳妥的路线。

---

### 6. 体相 Si 的 PBE0（非自洽）计算（PBE0 calculations for bulk silicon）

本教程在自洽 PBE 基态之上做**非自洽 PBE0**（通常用于得到更合理的带隙/能带修正，成本低于全自洽 hybrid）。

#### 6.1 PBE 与 PBE0 的 Γ‑Γ / Γ‑X 带隙

示例脚本（输出 `si-gaps.csv`）：

```python
# web-page: si-gaps.csv
from ase.build import bulk
from ase.parallel import paropen
from gpaw.hybrids.eigenvalues import non_self_consistent_eigenvalues
from gpaw import GPAW, PW

a = 5.43
si = bulk('Si', 'diamond', a)

fd = paropen('si-gaps.csv', 'w')

for k in range(2, 9, 2):
    name = f'Si-{k}'
    si.calc = GPAW(kpts={'size': (k, k, k), 'gamma': True},
                   mode=PW(200),
                   xc='PBE',
                   convergence={'bands': 5},
                   txt=name + '.txt')
    si.get_potential_energy()
    si.calc.write(name + '.gpw', mode='all')

    # Range of eigenvalues:
    n1 = 3
    n2 = 5

    ibzkpts = si.calc.get_ibz_k_points()
    kpt_indices = []
    for kpt in [(0, 0, 0), (0.5, 0.5, 0)]:  # Gamma and X
        # Find k-point index:
        i = abs(ibzkpts - kpt).sum(1).argmin()
        kpt_indices.append(i)

    # Do PBE0 calculation:
    epbe, vpbe, vpbe0 = non_self_consistent_eigenvalues(
        name + '.gpw',
        'PBE0',
        n1, n2,
        kpt_indices,
        snapshot=name + '.json')

    epbe0 = epbe - vpbe + vpbe0

    gg = epbe[0, 0, 1] - epbe[0, 0, 0]
    gx = epbe[0, 1, 1] - epbe[0, 0, 0]
    gg0 = epbe0[0, 0, 1] - epbe0[0, 0, 0]
    gx0 = epbe0[0, 1, 1] - epbe0[0, 0, 0]

    print(f'{k}, {gg:.3f}, {gx:.3f}, {gg0:.3f}, {gx0:.3f}', file=fd)
    fd.flush()

assert abs(gg - 2.559) < 0.01
assert abs(gx - 0.707) < 0.01
assert abs(gg0 - 3.873) < 0.01
assert abs(gx0 - 1.828) < 0.01
```

结果表（eV）：

| k-points | PBE(Γ-Γ) | PBE(Γ-X) | PBE0(Γ-Γ) | PBE0(Γ-X) |
| -------- | -------- | -------- | --------- | --------- |
| 2        | 2.451    | 0.602    | 4.188     | 2.152     |
| 4        | 2.542    | 0.691    | 3.952     | 1.911     |
| 6        | 2.556    | 0.705    | 3.891     | 1.847     |
| 8        | 2.559    | 0.707    | 3.873     | 1.828     |

> 实用解释
> `non_self_consistent_eigenvalues()` 返回三部分贡献：`eig_dft, vxc_dft, vxc_hyb`（这里命名为 `epbe, vpbe, vpbe0`）。
> 非自洽 hybrid 的本征值通常按：
> [
> \varepsilon^\text{hyb} = \varepsilon^\text{DFT} - v_{xc}^\text{DFT} + v_{xc}^\text{hyb}
> ]
> 对应代码里的 `epbe0 = epbe - vpbe + vpbe0`。

---

#### 6.2 晶格常数与体模量（Lattice constant and bulk modulus）

原文展示了用 ASE 数据库（`ase.db`）管理参数扫描，并对能量‑体积数据做 EOS 拟合（`ase.eos.EquationOfState`），得到晶格常数（以及拟合会给出 bulk modulus B）。

**计算并写入 `si.db`：**

```python
import ase.db
from ase.build import bulk
import numpy as np
from gpaw.hybrids.energy import non_self_consistent_energy as nsc_energy
from gpaw import GPAW, PW

a0 = 5.43

con = ase.db.connect('si.db')

for k in range(2, 9):
    for a in np.linspace(a0 - 0.04, a0 + 0.04, 5):
        id = con.reserve(a=a, k=k)
        if id is None:
            continue
        si = bulk('Si', 'diamond', a)
        si.calc = GPAW(kpts=(k, k, k),
                       mode=PW(400),
                       xc='PBE',
                       eigensolver='rmm-diis',
                       txt=None)
        si.get_potential_energy()
        name = f'si-{a:.2f}-{k}'
        si.calc.write(name + '.gpw', mode='all')
        epbe0 = nsc_energy(name + '.gpw', 'PBE0').sum()

        con.write(si, a=a, k=k, epbe0=epbe0)
        del con[id]
```

> 实用补充
>
> * `con.reserve(...)` 的写法适合并行/多进程跑参数扫描，避免重复计算同一参数点。
> * 实际工程中通常还会记录更多信息（如 ecut、smearing、收敛状态）以便追溯。

**拟合并绘图：**

```python
# web-page: si-a.png
import matplotlib.pyplot as plt
import ase.db
from ase.eos import EquationOfState


def lattice_constant(volumes, energies):
    eos = EquationOfState(volumes, energies)
    v, e, B = eos.fit()
    a = (v * 4)**(1 / 3)
    return a


con = ase.db.connect('si.db')
results = []
K = list(range(2, 9))
A = []
A0 = []
for k in K:
    rows = list(con.select(k=k))
    V = [row.volume for row in rows]
    E = [row.energy for row in rows]
    E0 = [row.epbe0 for row in rows]
    A.append(lattice_constant(V, E))
    A0.append(lattice_constant(V, E0))

print(K, A, A0)

plt.plot(K, A, label='PBE')
plt.plot(K, A0, label='PBE0')
plt.xlabel('number of k-points')
plt.ylabel('lattice constant [Ang]')
plt.savefig('si-a.png')
```

> 说明
> 这里 `a = (v * 4)**(1/3)` 的 “4” 来自 diamond 结构的体积与晶格常数换算关系（与所用晶胞/原胞定义有关）。若你改了结构/晶胞定义，这里的换算也需要相应调整。

---

## 总结：电子结构分析的常用“稳健做法”

为了让这份教程更适合做中文参考文档，这里给一份简短但实用的检查表：

1. **能量零点一致**：DOS/能带通常画 `E - E_F`；若用真空能级（vacuum level）对齐，确保体系确实有 vacuum（2D/slab）。
2. **DOS 先收敛 k 点，再谈细节**：尤其金属体相。
3. **投影方法要“知其局限”**：projector‑PDOS 可能有高能尾巴；需要时换 LCAO‑PDOS 或限制积分窗口。
4. **band structure 推荐 fixed_density**：先收敛基态，再沿路径求本征值；避免在 bandpath 上做全自洽（昂贵且不必要）。
5. **GLLB‑sc / GW / hybrid 的 band edges**：务必确保 VBM/CBM 在你的 k 网格或 bandpath 上“真的出现”，否则带隙会被误判。
6. **大方法（GW/RPA/unfolding）先做小体系 sanity check**：先确认脚本流程、对齐方式、输出含义，再上大体系/HPC。

如果你希望我把这一节再进一步“中文化工程整理”（例如：给出每种 DOS/PDOS 最常见的参数推荐区间、如何自动找 LCAO 轨道索引、如何判断 bandpath 是否包含 VBM/CBM 的程序化方法），我也可以在不改动原始示例代码的前提下，追加一个“实战附录”。

## 静电学与（电）动力学（Electrostatics and -dynamics）

> 本章聚焦于：**电势（electrostatic potential）**、**功函数（work function）**、**偶极层修正（dipole-layer correction / dipole correction）**、**均匀电子气（jellium）**、**氢原子裸库仑势（bare Coulomb potential）**、**μ子停驻位点（muon site）**、以及 **连续介质溶剂模型（Continuum Solvent Model, CSM）** 和 **溶剂化果冻模型（Solvated Jellium Method, SJM）** 的常电势计算工作流。
>
> 说明：文中保留了 GPAW/ASE 的关键英文术语与参数名（如 `pbc`, `poissonsolver`, `dipolelayer`, `PW-mode`, `fd`, `Fermi level` 等），便于直接对照代码与官方文档/源码。

---

### 1. GPAW 中的偶极层修正（Dipole-layer corrections）

#### 1.1 示例体系：Al(100) 双层薄膜 + 单侧 Na 吸附

构造一个 fcc(100) Al 的 2 层、$2\times2$ 表面薄膜（slab），在表面一侧放置一个 Na 吸附原子。示例使用 **实空间有限差分（`mode='fd'`）**。

```python
from ase.build import fcc100, add_adsorbate
from gpaw import GPAW

slab = fcc100('Al', (2, 2, 2), a=4.05, vacuum=7.5)
add_adsorbate(slab, 'Na', 4.0)
slab.center(axis=2)

slab.calc = GPAW(mode='fd',
                 txt='zero.txt',
                 xc='PBE',
                 setups={'Na': '1'},
                 kpts=(4, 4, 1))
e1 = slab.get_potential_energy()
slab.calc.write('zero.gpw')
```

**关键点解释：**

* `ase.build.fcc100()` 默认会在 **xy 平面周期（periodic）**，而在 z 方向通常用真空层隔开表面。
* 当体系在某方向 **非周期** 时（例如 `pbc=(True, True, False)`），GPAW（特别是在 real-space/fd 模式下）会对波函数在边界使用 **零边界条件（zero boundary conditions）**（可理解为 Dirichlet 边界：边界处波函数为 0），这会影响电势/电场在真空区的行为。
* 电势可以用 `get_electrostatic_potential()` 提取（见后续脚本示例）。

---

#### 1.2 全三维周期（3D PBC）会导致电势“被迫周期化”

如果把边界条件设为三维周期：

```python
slab.pbc = True
slab.calc = slab.calc.new(txt='periodic.txt')
e2 = slab.get_potential_energy()
```

此时：

* **静电势（electrostatic potential）必须周期**，并且其胞平均值通常会被归一为 0（“average to zero” 的常见表现）。
* 对于**非对称 slab（单侧吸附/带偶极）**，3D PBC 会引入“虚假的电场/电势斜坡”，使真空区电势难以变平，从而干扰功函数等量。

---

#### 1.3 为什么需要偶极修正：让真空区电场为零（平坦电势）

我们关心 slab 两侧的**功函数（work function, $\phi$）**，通常用：

[
\phi = V_{\text{vac}} - E_F
]

其中 $V_{\text{vac}}$ 是远离表面的真空电势平台（vacuum level），$E_F$ 是费米能级（Fermi level）。

要可靠读出 $V_{\text{vac}}$，需要真空区电势**平坦**（等价于真空区电场近似为 0）。对具有净偶极矩的 slab，可通过**偶极层修正（dipole correction）**实现：

```python
slab.pbc = (True, True, False)
slab.calc = slab.calc.new(poissonsolver={'dipolelayer': 'xy'},
                          txt='corrected.txt')
e3 = slab.get_potential_energy()
```

**参数含义：**

* `poissonsolver={'dipolelayer': 'xy'}`：在 **xy 平面**引入一个“偶极层”（dipole layer）用于抵消 slab 在 z 方向的偶极引起的平均电场，从而在真空区恢复平坦电势。
* 这是表面体系（尤其是不对称 slab）计算功函数、表面电势分布、带电表面等的常用处理。

---

#### 1.4 PW-mode（平面波）下的注意事项

在 **PW-mode**（`mode=PW(...)`）中，电势必须周期，因此修正后的电势形状会与 real-space 非周期边界的情况不同（可参考 [Bengtsson] 的讨论）。核心结论是：**PW-mode 的电势表现“必须周期”这一约束更强**，因此要格外注意功函数提取方式与参考零点。

---

#### 1.5 电势提取与作图脚本（教程原脚本，建议保留并复用）

> 下方脚本将：
> 1）读取 `zero/periodic/corrected/pwcorrected.gpw`；
> 2）取 `get_electrostatic_potential()` 并在 xy 平均；
> 3）与 `Fermi level` 一起画出随 z 的电势曲线；
> 4）对 corrected 情况在边界附近取真空电势点估算功函数；
> 5）额外输出 slab 的 POV-Ray 渲染。

```python
# web-page: zero.png, periodic.png, corrected.png, pwcorrected.png, slab.png
import numpy as np
import matplotlib.pyplot as plt
from ase.io import write
from gpaw import GPAW

# this test requires OpenEXR-libs

for name in ['zero', 'periodic', 'corrected', 'pwcorrected']:
    calc = GPAW(name + '.gpw', txt=None)

    efermi = calc.get_fermi_level()

    # Average over y and x:
    v = calc.get_electrostatic_potential().mean(1).mean(0)
    z = np.linspace(0, calc.atoms.cell[2, 2], len(v), endpoint=False)

    plt.figure(figsize=(6.5, 4.5))
    plt.plot(z, v, label='xy-averaged potential')
    plt.plot([0, z[-1]], [efermi, efermi], label='Fermi level')

    if name.endswith('corrected'):
        n = 6  # get the vacuum level 6 grid-points from the boundary
        plt.plot([0.2, 0.2], [efermi, v[n]], 'r:')
        plt.text(0.23, (efermi + v[n]) / 2,
                 r'$\phi$ = %.2f eV' % (v[n] - efermi), va='center')
        plt.plot([z[-1] - 0.2, z[-1] - 0.2], [efermi, v[-n]], 'r:')
        plt.text(z[-1] - 0.23, (efermi + v[-n]) / 2,
                 r'$\phi$ = %.2f eV' % (v[-n] - efermi),
                 va='center', ha='right')

    plt.xlabel(r'$z$, r$\AA$')
    plt.ylabel('(Pseudo) electrostatic potential, V')
    plt.xlim([0., z[-1]])
    if name == 'pwcorrected':
        title = 'PW-mode corrected'
    else:
        title = name.title()
    plt.title(title + ' boundary conditions')
    plt.savefig(name + '.png')

write('slab.pov',
      calc.atoms,
      rotation='-90x',
      show_unit_cell=2,
      povray_settings=dict(
          transparent=False,
          display=False)).render()
```

**实用提示（建议写入中文文档作为经验条目）：**

* `get_electrostatic_potential()` 得到的是 GPAW 定义下的（pseudo）电势；不同 setup、参考零点可能导致整体平移，但**功函数差（$V_{\text{vac}}-E_F$）**通常是稳健的。
* 真空电势平台建议取**远离表面**且**平坦区域**的平均，而不是单点；教程里用 `n=6` 只是演示。
* 对强偶极体系（吸附、外加场、带电 slab），偶极修正几乎是“必选项”。

---

### 2. 均匀电子气（Jellium）

本节复现实验性/经典理论工作 Lang 与 Kohn 的一些 **jellium** 计算（见 [Lang70]）。

#### 2.1 体相（Bulk）计算

目标：$r_s = 5$ Bohr。使用立方晶胞，晶格常数 $a=1.6$ Å，网格 `8*8*8`（对应 `h=0.2 Å`），k 点 `12*12*12`。

```python
import numpy as np
from ase import Atoms
from ase.units import Bohr
from gpaw.jellium import Jellium
from gpaw import GPAW, PW

rs = 5.0 * Bohr  # Wigner-Seitz radius
h = 0.2          # grid-spacing
a = 8 * h        # lattice constant
k = 12           # number of k-points (k*k*k)

ne = a**3 / (4 * np.pi / 3 * rs**3)
jellium = Jellium(ne)

bulk = Atoms(pbc=True, cell=(a, a, a))
bulk.calc = GPAW(mode=PW(400.0),
                 background_charge=jellium,
                 xc='LDA_X+LDA_C_WIGNER',
                 nbands=5,
                 kpts=[k, k, k],
                 h=h,
                 txt='bulk.txt')
e0 = bulk.get_potential_energy()
```

**解释：**

* $r_s$（Wigner–Seitz radius）定义了电子密度：$n = 3/(4\pi r_s^3)$（原子单位下）。
* `ne` 是该晶胞内的电子数（也是所需正背景电荷数），`background_charge=jellium` 相当于加入均匀正电背景以保证电中性。
* 输出中可看到晶胞包含约 `0.0527907` 个电子（取决于上述参数）。

---

#### 2.2 表面（Surface）计算

构造一个厚度 16.0 Å 的 jellium slab 放入盒子中，仅在 x/y 周期。盒子大小 `1.6*1.6*25.6 Å`，k 点 `12*12*1`：

```python
import numpy as np
from ase import Atoms
from ase.units import Bohr
from gpaw.jellium import JelliumSlab
from gpaw import GPAW, PW

rs = 5.0 * Bohr  # Wigner-Seitz radius
h = 0.2          # grid-spacing
a = 8 * h        # lattice constant
v = 3 * a        # vacuum
L = 10 * a       # thickness
k = 12           # number of k-points (k*k*1)

ne = a**2 * L / (4 * np.pi / 3 * rs**3)

eps = 0.001  # displace surfaces away from grid-points
jellium = JelliumSlab(ne, z1=v - eps, z2=v + L + eps)

surf = Atoms(pbc=(True, True, False),
             cell=(a, a, v + L + v))
surf.calc = GPAW(mode=PW(400.0),
                 background_charge=jellium,
                 xc='LDA_X+LDA_C_WIGNER',
                 kpts=[k, k, 1],
                 h=h,
                 convergence={'density': 0.001},
                 nbands=int(ne / 2) + 15,
                 txt='surface.txt')
e = surf.get_potential_energy()
surf.calc.write('surface.gpw')
```

**实用解释：**

* `JelliumSlab(ne, z1, z2)` 定义正背景只存在于 `z1~z2` 区间，从而形成“金属”区域；两侧留真空形成表面。
* `eps=0.001` 是数值技巧：把背景边界从网格点略微移开，避免表面正好落在网格点上引发数值振荡/收敛困难。

---

#### 2.3 表面能（surface energy）计算

表面能 $\sigma$ 的常见定义（对称双表面）：

[
\sigma = \frac{E_{\text{slab}} - N , E_{\text{bulk}}}{2A}
]

教程用等价写法（按几何比例换算 $N$）：

```python
from ase.io import read
e0 = read('bulk.txt').get_potential_energy()
e = read('surface.txt').get_potential_energy()
a = 1.6
L = 10 * a
sigma = (e - L / a * e0) / 2 / a**2
print('%.2f mev/Ang^2' % (1000 * sigma))
5.39 mev/Ang^2
print('%.1f erg/cm^2' % (sigma / 6.24150974e-5))
86.4 erg/cm^2
```

结果与 Lang–Kohn 的 `100 erg/cm^2` 相当接近（数量级一致、误差合理）。

---

#### 2.4 电子密度剖面（electron density profile）

```python
# web-page: fig2.png
import numpy as np
import matplotlib.pyplot as plt
from ase.units import Bohr
from gpaw import GPAW
rs = 5.0 * Bohr
calc = GPAW('surface.gpw')
density = calc.get_pseudo_density()[0, 0]
h = 0.2
a = 8 * h
v = 3 * a
L = 10 * a
z = np.linspace(0, v + L + v, len(density), endpoint=False)
# Position of surface is between two grid points:
z0 = (v + L - h / 2)
n = 1 / (4 * np.pi / 3 * rs**3)  # electron density
kF = (3 * np.pi**2 * n)**(1.0 / 3)
lambdaF = 2 * np.pi / kF  # Fermi wavelength
plt.figure(figsize=(6, 6 / 2**0.5))
plt.plot([-L / lambdaF, -L / lambdaF, 0, 0], [0, 1, 1, 0], 'k')
plt.plot((z - z0) / lambdaF, density / n)
# plt.xlim(xmin=-1.2, xmax=1)
plt.ylim(ymin=0)
plt.title(rf'$r_s={rs / Bohr:.1f}\ \mathrm{{Bohr}},\ '
          rf'\lambda_F={lambdaF / Bohr:.1f}\ \mathrm{{Bohr}}$')
plt.xlabel('DISTANCE (FERMI WAVELENGTHS)')
plt.ylabel('ELECTRON DENSITY')
plt.savefig('fig2.png')
```

---

#### 2.5 其它 jellium 几何

如果要做不同几何（非 slab / 非均匀分布等），需要继承 `Jellium` 并实现 `get_mask()` 方法来自定义背景电荷分布。

---

### 3. 氢原子裸库仑势（Bare Coulomb potential for hydrogen）

#### 3.1 使用平面波（plane waves）

GPAW 为氢提供一个特殊 setup：`ae`（all-electron）。它几乎不使用 PAW 的“重构魔法”，本质上就是裸库仑势 $-1/r$。

由于氢核处势能发散、波函数存在 cusp（尖点），用平面波展开时能量对截断能 `ecut` 的收敛会非常慢：

```python
from ase import Atoms
from gpaw import GPAW, PW
h = Atoms('H', cell=(5, 5, 5))
h.center()
for ecut in range(200, 1001, 100):
    h.calc = GPAW(setups='ae',
                  mode=PW(ecut),
                  txt=f'H-{ecut}-ae.txt')
    e = h.get_potential_energy()
```

查看收敛可用：

```bash
$ ase gui H.ae.txt
```

如果改用常规 PAW setup（例如实空间 fd）会快很多，把 `h.calc =` 替换为：

```python
h.calc = GPAW(mode='fd', txt='H.paw.txt')
```

---

#### 3.2 使用一维径向网格（1-d radial grid）

氢原子球对称，可用 1D 径向方程求解，命令行：

```bash
$ gpaw atom H -p
```

---

### 4. μ子停驻位点（Muon Site）

金属中植入的正 μ 子（positive muon）常停在对应**电子库仑势能（对电子而言）最大**的间隙位置。实践中常用 GPAW 计算得到的 **Hartree pseudo-potential**（在实现上常通过 `get_electrostatic_potential()` 得到的“伪静电势”近似）来估算停驻点：**寻找电势的极大值**。

本教程以 MnSi 为例，对比 A. Amato 等 [Amato14]：其通过 DFT 与实验分析得到 μ 子位点约为分数坐标 `(0.532, 0.532, 0.532)`。

---

#### 4.1 MnSi 计算（ASE + GPAW）

```python
from gpaw import GPAW, PW, MethfesselPaxton
from ase.spacegroup import crystal
from ase.io import write

a = 4.55643
mnsi = crystal(['Mn', 'Si'],
               [(0.1380, 0.1380, 0.1380), (0.84620, 0.84620, 0.84620)],
               spacegroup=198,
               cellpar=[a, a, a, 90, 90, 90])


for atom in mnsi:
    if atom.symbol == 'Mn':
        atom.magmom = 0.5

mnsi.calc = GPAW(xc='PBE',
                 kpts=(2, 2, 2),
                 mode=PW(800),
                 occupations=MethfesselPaxton(width=0.005),
                 txt='mnsi.txt')

mnsi.get_potential_energy()
mnsi.calc.write('mnsi.gpw')
v = mnsi.calc.get_electrostatic_potential()
write('mnsi.cube', mnsi, data=v)
```

这里输出 `mnsi.cube`（Gaussian cube）包含电势体数据（单位 eV），可用于 3D 可视化。

---

#### 4.2 用等值面（isosurface）寻找极大值（推荐的直观方法）

用外部可视化工具（例如 **Mayavi**；原文“majavi”通常指 mayavi）画等值面，取略低于最大值的等值：

```bash
$ python3 -m ase.visualize.mlab -c 11.1,13.3 mnsi.cube
```

`-c` 后是两条等值面电势值（最大值约 13.4 eV），也可用来识别次级局部极值。

---

#### 4.3 简化流程：画 2D 等高线辅助定位全局最大值

```python
# web-page: pot_contour.png
from gpaw import restart
import matplotlib.pyplot as plt
import numpy as np

mnsi, calc = restart('mnsi.gpw', txt=None)
v = calc.get_electrostatic_potential()
a = mnsi.cell[0, 0]
n = v.shape[0]
x = y = np.linspace(0, a, n, endpoint=False)

f = plt.figure()
ax = f.add_subplot(111)
cax = ax.contour(x, y, v[:, :, n // 2], 100)
cbar = f.colorbar(cax)
ax.set_xlabel('x (Angstrom)')
ax.set_ylabel('y (Angstrom)')
ax.set_title('Pseudo-electrostatic Potential')
f.savefig('pot_contour.png')
```

> 对比 [Amato14] 时要注意：示例为了省 CPU，k 点较少、截断能较低，但足够展示正确的极值位置；正式研究应提高收敛标准。

---

### 5. 连续介质溶剂模型（Continuum Solvent Model, CSM）

#### 5.1 理论背景（面向使用者的“够用版”理解）

GPAW 实现了一个 **连续介质溶剂模型（CSM）**。核心思想：

* 用实空间网格上的光滑空腔函数（smooth cavity）$g(\mathbf{r})$ 描述溶质/溶剂分区（“哪里是溶质、哪里是溶剂”不是硬边界，而是平滑过渡）。
* 溶剂的静电屏蔽用空间变化的介电常数（relative permittivity）$\epsilon(\mathbf{r})$ 表示：溶质区 $\epsilon$ 接近 1，溶剂区 $\epsilon$ 接近溶剂的体相介电常数（如水约 78）。

溶剂化自由能（solvation Gibbs energy）不仅来自静电，还包括：

* 空腔形成（cavity formation）
* 溶质-溶剂短程相互作用（short-range interactions）

因此模型允许加入额外的非静电相互作用项。

模型总的参数通常可归结为三类（以某套参数化为例）：

1. $\epsilon_{\infty}$：溶剂静态介电常数（直接来自实验）
2. $u_0$：与空腔/有效势相关，拟合实验体积
3. $\gamma$：与溶剂化自由能拟合相关（表面张力项等）

水（室温）参数见 Reference 1（文中提到 J. Chem. Phys. 141, 174108 (2014)）。DMSO 也有初步参数。经验上：中性分子水溶剂化自由能误差约 5 kJ/mol，阳离子约 13 kJ/mol；对阴离子往往不推荐直接使用该参数化（短程相互作用表征不足），除非你进一步改进/再参数化。

**可扩展性：**
GPAW 的 CSM 框架将空腔、介电函数、Poisson solver、非静电项等都实现为 Python 类，便于替换/扩展。但源代码里一些替代模型未充分测试，生产计算不建议直接用（可以作为二次开发参考）。

---

#### 5.2 用例：水中乙醇（Ethanol in Water）的溶剂化自由能

此处用“固定几何”的方式估算：

[
\Delta G_{\text{sol}} \approx E_{\text{solv}} - E_{\text{gas}}
]

（严格来说还可包含几何弛豫、热力学校正等；本教程以演示为主。）

```python
from gpaw import GPAW
from ase.build import molecule
from ase.units import mol, kJ, kcal, Pascal, m
from ase.parallel import parprint
from gpaw.solvation import (
    SolvationGPAW,             # the solvation calculator
    EffectivePotentialCavity,  # cavity using an effective potential
    Power12Potential,          # a specific effective potential
    LinearDielectric,  # rule to construct permittivity func from the cavity
    GradientSurface,  # rule to calculate the surface area from the cavity
    SurfaceInteraction  # rule to calculate non-electrostatic interactions
)

# all parameters on the user side of the solvation API follow the ASE
# unit conventions (eV, Angstrom, ...)

# non-solvent related DFT parameters
h = 0.2
vac = 5.0

# solvent parameters for water from J. Chem. Phys. 141, 174108 (2014)
u0 = 0.180  # eV
epsinf = 78.36  # dimensionless
gamma = 18.4 * 1e-3 * Pascal * m  # convert from dyne / cm to eV / Angstrom**2
T = 298.15  # Kelvin
atomic_radii = {'H': 1.09}


# create Atoms object for ethanol and add vacuum
atoms = molecule('CH3CH2OH')
atoms.center(vacuum=vac)

# perform gas phase calculation
atoms.calc = GPAW(mode='fd', xc='PBE', h=h, txt='gasphase.txt')
Egasphase = atoms.get_potential_energy()

# perform calculation with continuum solvent model from
# J. Chem. Phys. 141, 174108 (2014)
atoms.calc = SolvationGPAW(
    mode='fd', xc='PBE', h=h, txt='water.txt',
    cavity=EffectivePotentialCavity(
        effective_potential=Power12Potential(atomic_radii, u0),
        temperature=T,
        surface_calculator=GradientSurface()),
    dielectric=LinearDielectric(epsinf=epsinf),
    interactions=[SurfaceInteraction(surface_tension=gamma)])
Ewater = atoms.get_potential_energy()

# calculate solvation Gibbs energy in various units
DGSol_eV = Ewater - Egasphase
DGSol_kJ_per_mol = DGSol_eV / (kJ / mol)
DGSol_kcal_per_mol = DGSol_eV / (kcal / mol)

parprint('calculated Delta Gsol = %.0f meV = %.1f kJ / mol = %.1f kcal / mol' %
         (DGSol_eV * 1000., DGSol_kJ_per_mol, DGSol_kcal_per_mol))
```

期望结果：约 `-4.5 kcal/mol`；实验约 `-5.0 kcal/mol`。建议进一步阅读 `solvation module source code`（尤其是 `SolvationGPAW` 的参数传递与各组件类的含义）。

---

#### 5.3 快捷方式：直接获取 HW14 水参数

```python
import gpaw.solvation as solv

# ...

# convenient way to use HW14 water parameters:
calc = solv.SolvationGPAW(
    mode='fd', xc='PBE', h=0.2,  # non-solvent DFT parameters
    **solv.get_HW14_water_kwargs()
)

# ...
```

---

### 6. 溶剂化果冻模型（Solvated Jellium Method, SJM）

SJM 的典型目标：在**恒定电势（constant potential / grand-canonical）**条件下做表面反应、吸附、势垒（barrier）等计算，从而避免常电子数（canonical）计算中功函数随反应坐标变化 1–2 V 的问题。

---

#### 6.1 示例：Au(111) + 单个水分子覆盖层，在相对 SHE 为 -0.2 V 的电势

为了加速演示，用最小 slab + 单个水分子。这里要点是**电势输入必须用绝对标尺（absolute potential）**：
SHE 的绝对电势约 `4.4 V`，所以 `-0.2 V vs SHE` 对应 `4.2 V absolute`：

[
U_{\text{abs}} \approx 4.4\ \text{V} + U_{\text{vs SHE}}
]

脚本中使用 `sj={'target_potential': 4.2}` 即此含义。

同时，SJM 需要隐式溶剂来屏蔽电场，使用 GPAW 的 **CSM**（SJM 的计算器是 `SolvationGPAW` 的子类，溶剂参数直接传给计算器即可，不放入 `sj` 字典）。

```python
from ase.build import fcc111, molecule
from ase.units import Pascal, m

from gpaw.solvation.sjm import SJM, SJMPower12Potential
from gpaw.solvation import (
    EffectivePotentialCavity,
    LinearDielectric,
    GradientSurface,
    SurfaceInteraction)

# Build a tiny gold slab with a single water molecule above.
atoms = fcc111('Au', size=(1, 1, 3))
atoms.center(axis=2, vacuum=12.)
atoms.translate([0., 0., -4.])
water = molecule('H2O')
water.rotate('y', 90.)
water.positions += atoms[2].position + (0., 0., 4.4) - water[0].position
atoms.extend(water)

# Solvated jellium parameters.
sj = {'target_potential': 4.2}  # Desired potential

# Implicit solvent parameters (to SolvationGPAW).
epsinf = 78.36  # dielectric constant of water at 298 K
gamma = 18.4 * 1e-3 * Pascal * m
cavity = EffectivePotentialCavity(
    effective_potential=SJMPower12Potential(
        H2O_layer={'style': 'ghost_atoms'}),
    temperature=298.15,  # K
    surface_calculator=GradientSurface())
dielectric = LinearDielectric(epsinf=epsinf)
interactions = [SurfaceInteraction(surface_tension=gamma)]

# The calculator
calc = SJM(
    # General GPAW parameters.
    mode='fd',
    txt='Au111.txt',
    gpts=(16, 16, 136),
    kpts=(9, 9, 1),
    xc='PBE',
    maxiter=1000,
    # Solvated jellium parameters.
    sj=sj,
    # Implicit solvent parameters.
    cavity=cavity,
    dielectric=dielectric,
    interactions=interactions)
atoms.calc = calc

# Run the calculation.
atoms.get_potential_energy()
atoms.write('Au111.traj')
calc.write_sjm_traces(path='sjm_traces.out')  # .out for .gitignore
```

**你应该在输出 `Au111.txt` 中看到：**程序会调节 `excess electrons`（过剩电子数）直到目标电势在容差范围内（默认 `tol=0.01 V`）。首次通常需要多步，后续轨迹会更快（会复用电势-电荷斜率 `slope` 的信息）。示例中达到目标功函数仅需增加约 `0.007` 个电子（相对电中性体系）。

---

#### 6.2 为什么能量会变成 $\Omega = E - N\mu$（grand potential）

输出里会有类似信息：

```text
Legendre-transformed energies (grand potential, Omega = E - N mu)
 N (excess electrons):    +0.006526
 mu (-workfunction, eV):   -4.208990
 (Grand) free energy:    -23.628448
 (Grand) extrapolated:   -23.605878
```

* 这里的 $\Omega = E - N\mu$ 是 **Legendre 变换后的 grand potential**。
* `calc.get_potential_energy()`（以及 `atoms.get_potential_energy()`）默认返回的是这种 grand-canonical 能量（与恒电势体系的力一致）。
* 这使得它可以与 NEB、结构优化等方法兼容（前提是用对优化器，见后文）。
* 若你希望 `get_potential_energy()` 返回常电子数意义下的能量，可使用 `sj['grand_output']` 控制（具体含义以 SJM 文档/源码为准）。

---

#### 6.3 “sjm traces”：强烈建议检查（电势/背景电荷/空腔）

`calc.write_sjm_traces(...)` 会输出沿 z 的 $xy$ 平均剖面，帮助你确认：

* 溶剂空腔是否正确地从金属表面外开始
* jellium 背景电荷区域是否放置合理
* 电势边界与目标电势是否一致

生成的典型文件：

**potential.txt**

* $xy$ 平均电势，并以体系 `Fermi level` 为参考。外侧平台对应功函数；顶部边界电势应接近输入设定值。

**cavity.txt**

* $xy$ 平均空腔形状；与 `epsinf` 结合可理解为介电分布的剖面。

**background_charge.txt**

* $xy$ 平均的 jellium 背景电荷形状（归一化）。

> 备注：也可用 `calc.write_sjm_traces(style='cube')` 输出 3D cube 文件（非 1D trace）。

如果溶剂出现在不希望的区域，可用 ghost atoms 或 boundary plane 阻挡（见 `SJMPower12Potential` 的 docstring）；若溶剂“侵入金属区”，可按 docstring 增大原子半径等参数。

---

#### 6.4 恒电势结构弛豫：顺序（sequential）优化

在恒电势下做结构优化：示例在 Au(111) 上加 H 吸附，并固定部分原子：

```python
from ase import Atom
import ase.io
from ase.units import Pascal, m
from ase.optimize import BFGS
from ase.constraints import FixAtoms

from gpaw.solvation.sjm import SJM, SJMPower12Potential
from gpaw.solvation import (
    EffectivePotentialCavity,
    LinearDielectric,
    GradientSurface,
    SurfaceInteraction)

# Add an H adsorbate.
atoms = ase.io.read('Au111.traj')
atoms.append(Atom('H', atoms[2].position + (0., 0., 1.5)))

# Fix some atoms.
atoms.set_constraint(FixAtoms(indices=[0, 1]))

# Solvated jellium parameters.
sj = {'target_potential': 4.2}

# Implicit solvent parameters (to SolvationGPAW).
epsinf = 78.36  # dielectric constant of water at 298 K
gamma = 18.4 * 1e-3 * Pascal * m
cavity = EffectivePotentialCavity(
    effective_potential=SJMPower12Potential(
        H2O_layer={'style': 'ghost_atoms'}),
    temperature=298.15,  # K
    surface_calculator=GradientSurface())
dielectric = LinearDielectric(epsinf=epsinf)
interactions = [SurfaceInteraction(surface_tension=gamma)]

# The calculator
calc = SJM(
    # General GPAW parameters.
    mode='fd',
    txt='Au111-H-seq.txt',
    gpts=(16, 16, 136),
    kpts=(9, 9, 1),
    xc='PBE',
    maxiter=1000,
    # Solvated jellium parameters.
    sj=sj,
    # Implicit solvent parameters.
    cavity=cavity,
    dielectric=dielectric,
    interactions=interactions)
atoms.calc = calc

# Run the calculation.
opt = BFGS(atoms, trajectory='qn-Au111-H-seq.traj',
           logfile='qn-Au111-H-seq.log')
opt.run()
```

输出应类似常规弛豫，只是控制量从“电子数固定”变为“电势固定”。

---

#### 6.5 更快策略：结构与电势同时（simultaneous）优化

思路：不必每次都把电势严格收敛到很小容差才走离子步；可以在离子步的同时粗略调节电子数。做法：设置较松的 `tol`，并 `always_adjust=True`。

```python
from ase import Atom
import ase.io
from ase.units import Pascal, m
from ase.optimize import BFGS
from ase.constraints import FixAtoms

from gpaw.solvation.sjm import SJM, SJMPower12Potential
from gpaw.solvation import (
    EffectivePotentialCavity,
    LinearDielectric,
    GradientSurface,
    SurfaceInteraction)

# Add an H adsorbate.
atoms = ase.io.read('Au111.traj')
atoms.append(Atom('H', atoms[2].position + (0., 0., 1.5)))

# Fix some atoms.
atoms.set_constraint(FixAtoms(indices=[0, 1]))

# Solvated jellium parameters.
sj = {'target_potential': 4.2,
      'tol': 0.2,
      'always_adjust': True}

# Implicit solvent parameters (to SolvationGPAW).
epsinf = 78.36  # dielectric constant of water at 298 K
gamma = 18.4 * 1e-3 * Pascal * m
cavity = EffectivePotentialCavity(
    effective_potential=SJMPower12Potential(
        H2O_layer={'style': 'ghost_atoms'}),
    temperature=298.15,  # K
    surface_calculator=GradientSurface())
dielectric = LinearDielectric(epsinf=epsinf)
interactions = [SurfaceInteraction(surface_tension=gamma)]

# The calculator
calc = SJM(
    # General GPAW parameters.
    mode='fd',
    txt='Au111-H-sim.txt',
    gpts=(16, 16, 136),
    kpts=(9, 9, 1),
    xc='PBE',
    maxiter=1000,
    # Solvated jellium parameters.
    sj=sj,
    # Implicit solvent parameters.
    cavity=cavity,
    dielectric=dielectric,
    interactions=interactions)
atoms.calc = calc

# Run the calculation.
opt = BFGS(atoms, trajectory='qn-Au111-H-sim.traj',
           logfile='qn-Au111-H-sim.log')
opt.run()

# Tighten the tolerances again.
sj['tol'] = 0.01
sj['always_adjust'] = False
sj['slope'] = None
calc.set(sj=sj)
opt = BFGS(atoms, trajectory='qn-Au111-H-sim-1.traj',
           logfile='qn-Au111-H-sim.log')
opt.run()
```

**务必注意的经验：**

* 最后“收紧容差再跑一段”是很好的习惯：可确保最终结构确实在目标电势附近。
* 结构优化器建议用 **BFGS**，而不是 `BFGSLineSearch`：因为在优化结束前，能量与力在恒电势框架下可能存在轻微不一致，依赖能量的 line search 可能被误导。拿不准就用 BFGS。

---

#### 6.6 恒电势势垒：NEB / DyNEB

恒电势下找势垒的意义：避免 canonical 计算中功函数沿反应路径显著漂移的问题。SJM 让所有 NEB images 的功函数在用户容差内一致。

这里用一个简单扩散反应（吸附 H 从一个位点到另一个位点）演示。先准备初态和终态，然后用 DyNEB（更高效）：

```python
import ase.io
from ase.units import Pascal, m
from ase.optimize import BFGS
from ase.mep import interpolate, DyNEB

from gpaw.solvation.sjm import SJM, SJMPower12Potential
from gpaw.solvation import (
    EffectivePotentialCavity,
    LinearDielectric,
    GradientSurface,
    SurfaceInteraction)


def make_calculator(index):
    # Solvated jellium parameters.
    sj = {'target_potential': 4.2,
          'tol': 0.01,
          'excess_electrons': -0.01232,  # guess from previous
          'slope': -50.}  # guess from previous
    # Implicit solvent parameters (to SolvationGPAW).
    epsinf = 78.36  # dielectric constant of water at 298 K
    gamma = 18.4 * 1e-3 * Pascal * m
    cavity = EffectivePotentialCavity(
        effective_potential=SJMPower12Potential(
            H2O_layer={'style': 'ghost_atoms'}),
        temperature=298.15,  # K
        surface_calculator=GradientSurface())
    dielectric = LinearDielectric(epsinf=epsinf)
    interactions = [SurfaceInteraction(surface_tension=gamma)]
    calc = SJM(
        # General GPAW parameters.
        mode='fd',
        txt=f'gpaw-{index:d}.txt',
        gpts=(16, 16, 136),
        kpts=(9, 9, 1),
        xc='PBE',
        maxiter=1000,
        # Solvated jellium parameters.
        sj=sj,
        # Implicit solvent parameters.
        cavity=cavity,
        dielectric=dielectric,
        interactions=interactions)
    return calc


# Create the band of images, attaching a calc to each.
initial = ase.io.read('qn-Au111-H-sim-1.traj')
final = ase.io.read('qn-Au111-H-hollow-1.traj')
images = [initial]
for index in range(5):
    images += [initial.copy()]
    images[-1].calc = make_calculator(index + 1)
images += [final]
interpolate(images)

# Create and relax the DyNEB.
neb = DyNEB(images)
opt = BFGS(neb, logfile='neb.log', trajectory='neb.traj')
opt.run(fmax=0.05)
```

**重要实践建议：**

* `DyNEB` 比 `NEB` 更省：它会跳过已满足力阈值的 images 的重复计算。
* 对 grand-canonical（恒电势）方法，一般不推荐按 images 并行（parallel over images）：因为某一步可能只有一个 image 需要额外电势平衡，其他 images 会“干等”。串行反而能更好利用 DyNEB 的跳过机制。

---

#### 6.7 常电荷（constant-charge）模式

SJM 也支持固定过剩电子数 `excess_electrons`，用于快速扫一段电势范围（不强求精确电势值）。用法：

```python
sj = {'excess_electrons': 0.5}
calc = SJM(sj=sj, ...)
```

如果你确实想在**常电荷系综**下工作（例如电荷外推方案），通常希望输出 canonical 能量（与常电荷力一致）。可设置：

```python
sj = {'excess_electrons': ...,
      'grand_output': False} :
calc = SJM(sj=sj, ...)
```

---

## 小结：本章在实际材料计算中的“可直接复用要点”

* **表面功函数**：务必关注真空区电势是否平坦；非对称 slab 常用 `dipolelayer='xy'` 修正。
* **jellium**：`background_charge` 是“正背景”；表面 slab 用 `JelliumSlab` 并用 `eps` 避免边界落在网格点。
* **氢 ae（裸库仑）**：PW 收敛慢是正常现象；PAW/fd 会快很多；球对称问题可用 `gpaw atom -p`。
* **muon site**：用电势体数据（cube）找极大值，等值面是最直观的定位方式；演示参数不等于可发表参数。
* **CSM/SJM**：SJM 的 `target_potential` 用绝对标尺（SHE 约 4.4 V）；建议检查 traces；结构优化优先 BFGS；NEB 推荐 DyNEB 串行流程。

如果你希望我把这份译文进一步加工成“中文文档规范版”（例如：统一术语表、增加常见坑/FAQ、把每段代码配上中文注释与参数选择建议、给出功函数/电势读取的推荐后处理函数模板），我也可以直接在这份 Markdown 上继续迭代。




