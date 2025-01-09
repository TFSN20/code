# ASE
## 最佳体验：使用挂载
- cd 到挂载目标目录（注意D盘符小写）
  ```
  cd /mnt/d/OneDrive/Codes/school/计算
  ```
## GUI 操作
- 单独ase可在windows上运行，gpaw需要在linux系统上运行（如Ubuntu）。
- Ubuntu可视化需要ssh远程连接桌面，在Ubuntu上保存traj文件后，可在windows平台上使用```ase gui myatoms.traj```运行查看3D文件。
- 安装gpaw数据库
  ```
  gpaw install-data /home/tfsn20/paw_datasets
  ```
  等待下载，按下y
## 使用colab
- gpaw前奏
  ```
  !sudo apt update
  !sudo apt install -y build-essential python3-dev python3-pip libopenblas-dev libfftw3-dev libxc-dev libscalapack-mpi-dev
  ```
- 下载库
  ```
  !pip install ase gpaw
  ```
- 下载数据库
  ```
  !gpaw install-data /paw_datasets
  ```
- 查看cpu配置
  ```
  !lscpu
  ```
- 查看显卡配置
  ```
  !nvidia-smi
  ```
- root强制多核运行
  ```
  !mpiexec --oversubscribe -np 4 --allow-run-as-root gpaw python '/content/drive/MyDrive/cal/graphene_test/cal_graphene_colab.py'
  ```
- trajectory路径
  ```
  BFGS(system, trajectory='/home/tfsn20/python/'+f'{file_name_prefix}_{cutoff_energy}_{fmax}_{other_info}_optiming.traj')
  ```
- 当前用户或~目录
  colob是虚拟环境
  ```
  import getpass
  print(getpass.getuser())
  ```

### 问题
- 有些函数官方未实现cuda版
  ```
  /usr/local/lib/python3.10/dist-packages/gpaw/new/c.py in pwlfc_expand_gpu(f_Gs, emiGR_Ga, Y_GL, l_s, a_J, s_J, cc, f_GI, I_J)
     70                      l_s, a_J, s_J,
     71                      cc, f_GI, I_J):
  ---> 72     print(xp)
     73     f_GI = xp.empty((G2 - G1, self.nI), complex)
     74     I1 = 0

  NotImplementedError: 
  ```
## 例子
- 金线
  ```
  from ase import Atoms
  d = 2.9
  L = 10.0
  wire = Atoms('Au',
               positions=[[0, L / 2, L / 2]],
               cell=[d, L, L],
               pbc=[1, 0, 0])
  ```
  金金属为fcc堆积，晶胞长度为4.078 Å。固原子间之间最近的距离为面对角线上的4.078 Å*√2/2约等于2.9 Å.
  此时晶胞不可按照原来晶胞想象。
## 常用功能
- 晶胞，3D文件，动态优化3D文件，继续优化计算    
  - 设置虚拟晶胞
    ```
    water_cluster.set_cell([20, 20, 20])  # 足够大的立方晶胞
    water_cluster.center()  # 将体系置于晶胞中心
    water_cluster.set_pbc(False)  # 关闭周期性边界条件
    ```
  - 保存未优化的traj3D文件
    ```
    from ase.io import write
    write('water_cluster_no_optim.traj', water_cluster)
    ```
  - 设置 GPAW 计算器，保存计算器数据
    ```
    from gpaw import GPAW, PW
    calc = GPAW(mode=PW(300), xc='PBE', charge=0.0,txt='water_cluster_output.txt')
    water_cluster.calc = calc
    ```
  - 优化结构，保存优化中的动态traj3D文件
    ```
    from ase.optimize import BFGS
    dyn = BFGS(water_cluster, trajectory='water_cluster_optiming.traj')
    dyn.run(fmax=0.1)
    ```
  - 获取能量
    ```
    water_cluster.get_potential_energy()
    ```
  - 继续计算
    ```
    from ase.io import Trajectory
    from ase.optimize import BFGS
    from gpaw import GPAW, PW
    from ase.visualize import view
    from ase.constraints import Hookean, FixAtoms
    path_cal_res = os.path.dirname(os.path.abspath(__file__))
    
    # 设置 GPAW 计算器
    calc = GPAW(mode=PW(300), xc='PBE', charge=0.0, txt='zn2_plus_ion_output_keeping.txt')
    # 加载轨迹文件并获取最后一步的结构
    traj = Trajectory(path_cal_res / Path('g_C-COOH_constraint_zn2+_no_constraint_optiming_keeping - 副本.traj'), 'r')  # 以只读模式加载轨迹
    system = traj[-1]  # 获取轨迹文件中的最后一步结构
    
    # 重新设置计算器
    system.calc=calc
    print(system)
    
    # 继续优化并保存到新的轨迹文件
    dyn = BFGS(system, trajectory='zn2_plus_ion_optiming_keeping.traj')
    dyn.run(fmax=0.05)
    ```
  - 终端查看
    ```
    ase gui zn_h2o_cluster_no_optim.traj
    ```
  - 代码里查看
    ```
    from ase.visualize import view
    view(atoms)
    ```
- 并行计算    
  - 计算步骤中有保存文件时可能需要等待一会才会保存
  - 打印也有延迟，bug？？？ 
  - 使用逻辑cpu
    ```
    mpiexec --use-hwthread-cpus -np 4 gpaw python cal_zn2+_6h2o.py
    ```
  - 强制使用
    ```
    mpiexec --oversubscribe -np 4 gpaw python cal_zn2+_6h2o.py
    ```
- gpu计算
  - 仅仅支持PW
  - 示例代码
    ```
    from gpaw.new.ase_interface import GPAW
    GPAW(mode=PW(cutoff_energy), xc='PBE', charge=2.0, parallel={'domain': 1, 'gpu':True},txt=f'output.txt')
    ```
- 修改模型
  - 运行ase gui 文件名.traj，直接修改，注意增加原子时坐标单位为nm.可基于全局坐标或选中原子坐标进行偏移
- 保存图片
  - 不能保存为jpg
  - 格式为名称.png@1:23
## 常见问题
### 粒子之间的距离
### 采用哪种优化方法
- 经典EMT potential很快，但对于分子不准确。
### 平面波截断能量
- 越高的平面波截断能量需要的计算时间越长
- 过渡元素一半需要较高的平面波截断能量400，450eV，轻元素，如H2O需要较低的即可。
### k点网格
- 用于描述布里渊区（Brillouin Zone）的取样密度，主要用于周期性体系的电子结构计算（如能带、总能量、态密度计算等）。
- 周期性 vs 非周期性总结    
  ![image](https://github.com/user-attachments/assets/59f31f50-0d7a-4bdc-a7ad-70f75f9d254f)    
  手动设置 vs 默认 kpts 的区别    
  ![image](https://github.com/user-attachments/assets/5a5d4854-0319-4100-8d8f-17ca846071c4)
- 石墨烯这种假如是5x5x1，加了kZ值计算是真的稳定，但是计算时间变长。
### 测试cutoff_energy时不需要优化函数吗
- 在测试 cutoff energy 时，不需要优化函数，通常是直接计算体系的总能量，因为此时的目的是评估平面波截断能量对总能量的影响，而不是对结构进行优化。
- 目标：确定平面波截断能量𝐸cutoff的收敛性，即找到一个足够高的 Ecutoff，使总能量收敛到稳定值。
- 方法：测试不同 𝐸cutoff 下的总能量，并观察能量变化趋势。
- 原因：
  - 结构优化（如使用 BFGS 或 FIRE）会改变原子位置，从而影响能量结果。
  - 结构优化会引入其他变量，使得 cutoff energy 的测试结果不单纯。
  - 因此，截断能量的测试应在固定的几何结构下进行，确保能量变化仅与𝐸cutoff相关。
- 2. 截断能量测试的步骤
  - 固定几何结构：使用初始结构（可以是合理猜测的结构）进行总能量计算，不改变原子位置。
  - 逐步增加 cutoff energy：设置一系列逐渐增大的 𝐸cutoff。
  - 计算总能量：对每个 𝐸cutoff值计算体系的总能量。
  - 分析收敛性：观察能量变化，当能量变化小于设定阈值（如 <1meV/atom）时，认为收敛。
# 机器学习原子间势
- 大多数机器学习原子间MLIP架构，没有明确地包含原子电荷和磁性。
- 对于广泛的材料科学问题，准确描述这些长程相互作用，是必要的，包括反应系统和界面、电子驱动相变、磁性材料等。
- 将电荷和磁性结合到机器学习原子间势MLIPS中的主要挑战，不在于架构开发，而在于训练数据的生成和计算效率。
## MACE
- mace-torch
- 优点是几乎不需要设置什么参数，只需要指定原子类型和位置即可，且计算速度极快。
- 缺点是计算不准确，且无法计算含有电荷的体系。
- 凭借其优点，我们可以快速确定使用GPAW计算时的一些参数，如原子收敛位置。
# Bader
## 下载安装
- 链接：https://theory.cm.utexas.edu/henkelman/code/bader/
- 解压
  ```
  tar -zxvf bader_lnx_64.tar.gz
  ```
- 设置为可执行文件
  ```
  chmod +x bader
  ```
- 添加到command可调用，在~/.bashrc文件，目录是文件的上级目录，即不包含文件每次，只是上级目录，即最后没有/
  ```
  export PATH=~/bader:${PATH}（一半新建bader目录将bader文件放进去）或
  export PATH=~:${PATH}（一半不用，因为这相当于把整个用户目录放入环境变量中）
  ```


