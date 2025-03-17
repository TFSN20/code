# 安装
## 下载软件
- vscode最新版本
- Ubuntu 22.04TLS任意版本
- 电脑开启开发者模式，开启windows功能中的适用于Linux的子系统
- 更改apt源，将Uunbutu系统根目录发送到桌面，方便定位。
- Ubuntu系统自带了Python 3，无需下载，下载pip
  ```
  sudo apt install python3-pip
  ```
- pip换源
  ```
  pip config set global.index-url https://mirrors.aliyun.com/pypi/simple/
  ```
- 下载ase
  ```
  pip install --upgrade ase
  ```
- 下载gpaw
  ```
  pip install gpaw
  ```
  下载出错多试几次
  - fatal error: xc.h: No such file or directory  
         14 | #include <xc.h> // If this file is not found, install libxc https://gpaw.readthedocs.io/install.html#libxc-installation   
    出现这个说明需要
    ```
    sudo apt-get install libxc-dev
    ```
  - /usr/bin/ld: cannot find -lblas: No such file or directory  
      collect2: error: ld returned 1 exit status  
      error: command '/usr/bin/x86_64-linux-gnu-gcc' failed with exit code 1
    说明需要
    ```
    sudo apt-get install libblas-dev libatlas-base-dev
    sudo apt install gcc
    ```
    
- 下载mace
  ```
  pip install mace-torch
  ```
  - 出现
    ERROR: Program or command 'python3' not found or not executable  
    metadata-generation-failed  
    需要在python.exe目录将python.exe复制为python3.exe
    完成之后还会
    nstalling build dependencies ... done  
    Getting requirements to build wheel ... done   
    Preparing metadata (pyproject.toml) ... error   
    error: subprocess-exited-with-error  

    × Preparing metadata (pyproject.toml) did not run successfully.   
     ..\meson.build:6:13: ERROR: Command `E:\Softwares\Python312\python3.EXE discover_version.py` failed with status 1.
    说明python版本可能过高比如3.13对于mace0.3.12就不行，有上面问题。
    
  

# 计算常识
## 单位
- 1 ev = 96 kJ/mol = 23.06 kcal/mol
- 1 kJ/mol = 0.238846 kcal/mol
- xyz单位为Å体积水溶液中分子数：xyz乘以0.999 g/ml / 18 g/mol ✖6.02✖10^23 = 615.8  （A是埃,xyz is 18✖32✖32）
- 水分子间距3.1Å，石墨烯六元环边长是1.42Å
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
- bader和电子密度
  - 电子密度cube
    ```
    import os
    from ase.build import molecule
    from ase.io import write
    from ase.units import Bohr
    from gpaw import GPAW, PW
    from gpaw.analyse.hirshfeld import HirshfeldPartitioning
    from ase.io import Trajectory
    from pathlib import Path
    
    traj_num = 19
    
    path_big_file = '/mnt/d/cal'
    path_cal_res = os.path.dirname(os.path.abspath(__file__))
    new_dir_path = os.path.join(path_big_file, os.path.basename(path_cal_res))
    os.makedirs(new_dir_path, exist_ok=True)
    
    traj = Trajectory(path_cal_res / Path("g_C-OH_constraint_optiming_keeping.traj"), 'r') 
    system = traj[traj_num]
    print(system)
    # calc = GPAW(mode=PW(400), xc='PBE', charge=+2.0, kpts=(2,2,1), txt=path_cal_res / Path('system1_output1.txt'))
    calc = GPAW(mode='fd', xc='PBE', kpts=(2,2,1), txt=path_cal_res / Path('system1_output1_fd1.txt'))
    system.calc = calc
    print(system.get_potential_energy())
    calc.write(new_dir_path / Path('system1_fd1.paw'), mode='all')
    
    
    # write Hirshfeld charges out
    hf = HirshfeldPartitioning(system.calc)
    for atom, charge in zip(system, hf.get_charges()):
        atom.charge = charge
    # system.write('Hirshfeld.traj') # XXX Trajectory writer needs a fix
    system.copy().write(new_dir_path / Path('Hirshfeld_fd1.traj'))
    
    # create electron density cube file ready for bader
    rho = system.calc.get_all_electron_density(gridrefinement=4)
    write(new_dir_path / Path('density.cube'), system, data=rho * Bohr**3)
    ```
  - 分析密度立方体文件
    ```
    ~/bader -p all_atom -p atom_index density.cube
    ```
  - 中间切面数据by density.cube and AtIndex.cube
    ```
    import os
    import pickle
    import numpy as np
    import matplotlib.pyplot as plt
    from ase.io.cube import read_cube_data
    from pathlib import Path
    
    section = 'x'
    atom_index_of_interest = 1 #图像中心是 原子index为1
    min_den = 0.1
    max_den = 1
    split_num = 200
    xmin=-5
    xmax=5
    ymin=-5
    ymax=5
    
    path_cal_res = os.path.dirname(os.path.abspath(__file__))
    dens, atoms = read_cube_data(path_cal_res / Path('density.cube'))
    bader, atoms = read_cube_data(path_cal_res / Path('AtIndex.cube'))
    
    x0, y0, z0 = atoms.positions[atom_index_of_interest]
    if section == 'x':
        x = len(dens) // 2
        dens = dens[x]
        bader = bader[x]
        y = np.linspace(0, atoms.cell[1, 1], len(dens), endpoint=False) - y0
        z = np.linspace(0, atoms.cell[2, 2], len(dens[0]), endpoint=False) - z0
    
        print(y.shape, z.shape, dens.shape, bader.shape)
        print(atoms.positions)
        print(dens.min(), dens.mean(), dens.max())
        print(bader.min(), bader.mean(), bader.max())
        plt.figure(figsize=(5, 5))
    
        # Add color bar
        cbar = plt.colorbar(plt.contourf(z, y, dens, np.linspace(min_den, max_den, split_num)))
        cbar.ax.set_ylabel('Density', rotation=-90, va="bottom")
        plt.contour(z, y, bader, [1.75], colors='k')
    
        # 确保 x 和 y 坐标轴等比例
        plt.gca().set_aspect('equal', adjustable='box')
        plt.axis(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
        plt.show()
    elif section == 'y':
        # 取 y 中间切面
        y = len(dens[0]) // 2
        dens = dens[:, y, :]
        bader = bader[:, y, :]
        
        # 定义 x 和 z 坐标
        x = np.linspace(0, atoms.cell[0, 0], len(dens), endpoint=False) - x0
        z = np.linspace(0, atoms.cell[2, 2], len(dens[0]), endpoint=False) -z0
        
        print(x.shape, z.shape, dens.shape, bader.shape)
        plt.figure(figsize=(5, 5))
        
        # 绘制密度图
        cbar = plt.colorbar(plt.contourf(z, x, dens, np.linspace(min_den, max_den, split_num)))  # 注意：x 和 z 对应坐标
        cbar.ax.set_ylabel('Density', rotation=-90, va="bottom")
        # 绘制等值线
        plt.contour(z, x, bader, [1.5], colors='k')
    
        # 确保 x 和 y 坐标轴等比例
        plt.gca().set_aspect('equal', adjustable='box')
        # 设置坐标轴范围
        plt.axis(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
        plt.show()
        # plt.savefig('h2o-bader-ycut.png')
    elif section == 'z':
        z = len(dens[0][0]) // 2
        dens = dens[:, :, z]
        bader = bader[:, :, z]
    
        # 定义 x 和 y 坐标
        x = np.linspace(0, atoms.cell[0, 0], len(dens), endpoint=False) - x0
        y = np.linspace(0, atoms.cell[1, 1], len(dens[0]), endpoint=False) -y0
    
        print(x.shape, y.shape, dens.shape, bader.shape)
        plt.figure(figsize=(5, 5))
    
        # 绘制密度图
        cbar = plt.colorbar(plt.contourf(x, y, dens.T, np.linspace(min_den, max_den, split_num)) ) # 转置数据以适应坐标系
        cbar.ax.set_ylabel('Density', rotation=-90, va="bottom")
        # 绘制等值线
        plt.contour(x, y, bader.T, [1.5], colors='k')
    
        # 确保 x 和 y 坐标轴等比例
        plt.gca().set_aspect('equal', adjustable='box')
    
        # 设置坐标轴范围
        plt.axis(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)
        plt.show()
    ```
  - 导出三维x, y, z, den数据
    ```
    import numpy as np
    import pandas as pd
    from ase.io.cube import read_cube_data
    from pathlib import Path
    import os
    
    # 读取数据
    path_cal_res = os.path.dirname(os.path.abspath(__file__))
    dens, atoms = read_cube_data(path_cal_res / Path('density.cube'))
    
    # 获取数组的形状
    nx, ny, nz = dens.shape  # (160, 208, 176)
    
    # 生成 x, y, z 坐标
    x = np.linspace(0, atoms.cell[0, 0], nx, endpoint=False)
    y = np.linspace(0, atoms.cell[1, 1], ny, endpoint=False)
    z = np.linspace(0, atoms.cell[2, 2], nz, endpoint=False)
    
    # 创建网格并展平
    x_grid, y_grid, z_grid = np.meshgrid(x, y, z, indexing='ij')
    
    # 展平网格和强度数据
    x_flat = x_grid.ravel()
    y_flat = y_grid.ravel()
    z_flat = z_grid.ravel()
    dens_flat = dens.ravel()  # 强度值（密度）
    
    # 筛选强度在范围 [0.01, 1.2] 内的点
    mask = (dens_flat >= 0.01) & (dens_flat <= 1.2)
    
    # 只保留符合条件的坐标和强度
    x_filtered = x_flat[mask]
    y_filtered = y_flat[mask]
    z_filtered = z_flat[mask]
    dens_filtered = dens_flat[mask]
    
    # 创建 DataFrame 并保存为 CSV
    data = np.column_stack((x_filtered, y_filtered, z_filtered, dens_filtered))
    df = pd.DataFrame(data, columns=['x', 'y', 'z', 'intensity'])
    
    # 保存为 CSV 文件
    output_file = 'filtered_density_data.csv'
    df.to_csv(output_file, index=False)
    
    print(f"筛选后的数据已保存至 {output_file}")
    ```
  - 绘制三维密度图
    ```
    import numpy as np
    import pandas as pd
    from mayavi import mlab
    
    # 加载筛选后的数据
    data = pd.read_csv('filtered_density_data.csv')
    
    # 提取坐标和强度
    x = data['x'].values
    y = data['y'].values
    z = data['z'].values
    intensity = data['intensity'].values
    
    """ color """
    def generate_lut(start_rgba, end_rgba, n_colors=256):
        """
        生成颜色查找表 (LUT)。
        
        参数：
        - start_rgba: list or tuple，起始颜色和透明度 [R, G, B, A]，每个值范围 0-255。
        - end_rgba: list or tuple，末尾颜色和透明度 [R, G, B, A]，每个值范围 0-255。
        - n_colors: int，生成的渐变分段数，默认 256。
    
        返回：
        - lut: numpy.ndarray，大小为 (n_colors, 4) 的 LUT 数组。
        """
        # 创建线性渐变的 LUT
        start_rgba = np.array(start_rgba, dtype=np.float32)
        end_rgba = np.array(end_rgba, dtype=np.float32)
        lut = np.zeros((n_colors, 4), dtype=np.uint8)
    
        # 线性插值计算每个通道的值
        for i in range(4):  # 对 RGBA 四个通道分别计算
            lut[:, i] = np.linspace(start_rgba[i], end_rgba[i], n_colors)
    
        return lut
    
    
    
    # 绘制三维点云，点的大小和颜色根据强度值变化
    mlab.figure('Density Visualization', bgcolor=(1, 1, 1), size=(800, 600))
    _ = mlab.points3d(x, y, z, intensity, scale_mode='none', scale_factor=0.2, colormap='viridis')
    
    
    # 定义起始和末尾颜色以及透明度
    start_color = [68, 1, 84, 0]   # 起始 [R, G, B, A]
    end_color = [253, 231, 37, 255]   # 结束 [R, G, B, A]
    n_steps = abs(20)  # 分段数
    # 生成 LUT
    lut = generate_lut(start_color, end_color, n_steps)
    
    # 修改lut
    _.module_manager.scalar_lut_manager.lut.table = lut
    
    # 添加坐标轴
    mlab.axes(xlabel='X', ylabel='Y', zlabel='Z', color=(0, 0, 0))
    
    # 添加标题
    mlab.title('3D Density Plot', size=0.5)
    
    # 显示图形
    mlab.show()
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
- gpaw的get_potential_energy在第一次计算时输出的是Free energy，以后再次访问traj获取的是Extrapolated能量。有时这两者不相等。
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
### 正交坐标系与石墨，石墨烯
- 正交转石墨烯坐标系：参考https://wiki.fysik.dtu.dk/ase/_modules/ase/visualize/mlab.html，https://github.com/enthought/mayavi/issues/1268
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
- 错误：TypeError: unhashable type: 'PolyData'
  ```
  pip install https://github.com/enthought/mayavi/zipball/main
  ```
## 电子密度,bader电荷，静电势
### density.cube代码（h=0.2和gridrefinement=2对应0.2/2=0.1，对应静电势）
  ```
  import os
  from ase.build import molecule
  from ase.io import write
  from ase.units import Bohr
  from gpaw import GPAW, PW
  from gpaw.analyse.hirshfeld import HirshfeldPartitioning
  from ase.io import Trajectory
  from pathlib import Path
  
  # 选择轨迹名称和index
  traj_name = 'g_C-F_zn2+_constraint_optiming_keeping2.traj'
  traj_index = -1
  # 有时一个目录下可能有许多traj文件都需要电子密度 需要设置文件夹flag
  flag='zn'
  # 设置计算器的总电荷，除此之外还需要设置计算器的kpts
  charge=+2.0
  
  
  path_big_file = '/mnt/d/cal'
  path_cal_res = os.path.dirname(os.path.abspath(__file__))
  # 有时一个目录下可能有许多traj文件都需要电子密度 则os.path.basename(path_cal_res)+flag
  new_dir_path = os.path.join(path_big_file, os.path.basename(path_cal_res)+flag)
  os.makedirs(new_dir_path, exist_ok=True)
  
  traj = Trajectory(path_cal_res / Path(traj_name), 'r') 
  system = traj[traj_index]
  
  # system.pbc = False
  # calc = GPAW(mode=PW(400), xc='PBE', h=0.2, charge=+2.0, kpts=(2,2,1), txt=path_cal_res / Path('system1_output1.txt'))
  calc = GPAW(mode=PW(400), xc='PBE', h=0.2, charge=charge, kpts=(2,2,1), txt=new_dir_path / Path('system_output.txt'))
  # 'fd' 模式下 设置pcb false 对于原本是orthogonal的系统会有问题，而pw模式下无法对电荷求解，各有劣势
  # 由于一般是pw模式，所以需要注释掉电荷相关
  system.calc = calc
  print(system.calc)
  print(system)
  print(system.get_potential_energy())
  # calc.write(new_dir_path / Path('system.paw'), mode='all')
  
  
  # write Hirshfeld charges out
  # hf = HirshfeldPartitioning(system.calc)
  # for atom, charge in zip(system, hf.get_charges()):
  #     atom.charge = charge
  # # system.write('Hirshfeld.traj') #  Trajectory writer needs a fix
  system.copy().write(new_dir_path / Path('system.traj'))
  
  # create electron density cube file ready for bader
  rho = system.calc.get_all_electron_density(gridrefinement=2)
  write(new_dir_path / Path('den.cube'), system, data=rho * Bohr**3)
  ```
### esp.cube代码
  ```
  from ase.build import molecule
  from ase.io import write
  from ase.units import Bohr
  from gpaw import GPAW,PW
  from ase.io import Trajectory
  
  from pathlib import Path
  import os
  
  
  # 选择轨迹名称和index
  traj_name = 'g_C-F_zn2+_constraint_optiming_keeping2.traj'
  traj_index = -1
  # 有时一个目录下可能有许多traj文件都需要电子密度 需要设置文件夹flag
  flag='zn'
  # 设置计算器的总电荷，除此之外还需要设置计算器的kpts
  charge=+2.0
  
  
  path_big_file = '/mnt/d/cal'
  path_cal_res = os.path.dirname(os.path.abspath(__file__))
  # 有时一个目录下可能有许多traj文件都需要电子密度 则os.path.basename(path_cal_res)+flag
  new_dir_path = os.path.join(path_big_file, os.path.basename(path_cal_res)+flag)
  os.makedirs(new_dir_path, exist_ok=True)
  
  traj = Trajectory(path_cal_res / Path(traj_name), 'r') 
  system = traj[traj_index]
  
  system.calc = GPAW(xc='PBE',
                   mode=PW(400),
                   charge=charge,
                   h=0.2,
                   kpts=(2,2,1),
                   txt=new_dir_path / Path('pw.txt'))
  print(system)
  print(system.get_potential_energy())
  system.calc.write(new_dir_path / Path('pw.gpw'))
  v = system.calc.get_electrostatic_potential()
  write(new_dir_path / Path('esp.cube'), system, data=v)
  ```
### esp-show代码
  ```
  import numpy as np
  from ase.io.cube import read_cube_data
  from pathlib import Path
  import os
  from mayavi import mlab
  from ase.data import covalent_radii
  from ase.data.colors import cpk_colors
  
  atoms_magnification = 1
  
  # Function to extract points and ESP values on a specific density isosurface
  def extract_points_and_esps_on_isosurface(data, esps, iso_value):
    """
    Extract points and ESP values on a specific density isosurface.

    Parameters:
        data: 3D ndarray
            The density data.
        esps: 3D ndarray
            The electrostatic potential data.
        iso_value: float
            The density isosurface value.

    Returns:
        points: ndarray
            The coordinates of the points on the isosurface.
        esp_values: ndarray
            The ESP values at the points on the isosurface.
    """
    # Find points on the isosurface (allow for a small tolerance)
    tolerance = 1e-4
    tolerance = 10**(-3.5)
    # print(10**(-4),10**(-3.6),10**(-3.5),10**(-3.4),10**(-3.3))
    isosurface_indices = np.where((data >= iso_value - tolerance) & (data <= iso_value + tolerance))

    # Extract coordinates of the points
    points = np.column_stack(isosurface_indices)

    # # Filter points based on z-coordinate
    # valid_indices = points[:, 2] <= 200
    # points = points[valid_indices]

    # # Update isosurface_indices accordingly
    # isosurface_indices = tuple(index[valid_indices] for index in isosurface_indices)

    # Extract corresponding ESP values
    esp_values = esps[isosurface_indices]

    return points, esp_values
  def plot_with_surface(atoms, data, points, esp_values, iso_value):
    """
    Plot atoms, unit-cell, iso-surfaces, and ESP surface using Mayavi.

    Parameters:
        atoms: Atoms object
            Positions, atomic numbers, and unit-cell.
        data: 3D ndarray of float
            Data for iso-surfaces.
        points: ndarray
            Coordinates of the points on the isosurface.
        esp_values: ndarray
            ESP values at the points on the isosurface.
        iso_value: float
            The density isosurface value.
    """
    mlab.figure(1, bgcolor=(1, 1, 1))  # Make a white figure
    # Plot the atoms as spheres:
    for pos, Z in zip(atoms.positions, atoms.numbers):
        mlab.points3d(*pos,
                      scale_factor=covalent_radii[Z]*atoms_magnification,
                      resolution=20,
                      color=tuple(cpk_colors[Z]))
    # Plot the density isosurface
    cp_density = mlab.contour3d(data, contours=[iso_value], transparent=True,
                                opacity=0.0, colormap='viridis')
    A = atoms.cell
    transformed_points = np.dot(points - 1, A / np.array(data.shape)[:, np.newaxis])

    # Use transformed points in mlab.points3d
    cp_surface = mlab.points3d(transformed_points[:, 0], transformed_points[:, 1], transformed_points[:, 2],
                               esp_values,
                               scale_mode='none', 
                               scale_factor=0.3,
                               colormap='RdBu', 
                               mode="sphere",
                               resolution=12,
                               opacity=0.03)

    # cp_surface.module_manager.scalar_lut_manager.data_range = [-15, -5]
    print(esp_values.min(),esp_values.max())
    colorbar = mlab.colorbar(cp_surface, title='esp', orientation='vertical', nb_labels=11)
    # 设置刻度数字颜色为黑色
    colorbar.scalar_bar.unconstrained_font_size = True  # 确保字体大小可调
    colorbar.label_text_property.font_size = 14  # Increase font size
    colorbar.label_text_property.color = (0, 0, 0)  # RGB格式，黑色为 (0, 0, 0)
    colorbar.title_text_property.color = (0, 0, 0)  # Black color for the title
    colorbar.title_text_property.font_size = 0  # Increase font size
    # Combine cp_density and cp_surface for tvtk processing
    cp = [cp_density]

    # Do some tvtk magic in order to allow for non-orthogonal unit cells
    A = atoms.cell
    for c in cp:
        polydata = c.actor.actors[0].mapper.input
        pts = np.array(polydata.points) - 1
        polydata.points = np.dot(pts, A / np.array(data.shape)[:, np.newaxis])

    # Show the 3D plot
    mlab.show()
  
  
  def main():
    # Define the path to the current script directory
    path_cal_res = Path(os.path.dirname(os.path.abspath(__file__)))

    # Read density and ESP data
    dens, atoms = read_cube_data(path_cal_res / 'den.cube')
    esps, _ = read_cube_data(path_cal_res / 'esp.cube')
    esps = esps * -1
    # Extract points and ESP values on the density isosurface (den = 0.001)
    iso_value = 0.001
    points, esp_values = extract_points_and_esps_on_isosurface(dens, esps, iso_value)

    # Print the points and ESP values
    print(f"Points on the density isosurface (den = {iso_value}):")
    print(points)
    print(f"ESP values on the density isosurface (den = {iso_value}):")
    print(esp_values)
    print(len(points))
    print(len(esp_values))

    # Plot with surface
    plot_with_surface(atoms, dens, points, esp_values, iso_value)
  
  if __name__ == "__main__":
      main()
  ```
### 分析密度立方体文件
  ```
  ~/bader -p all_atom -p atom_index density.cube
  ```
### 三维多等值面电子密度图像
  ```
  import numpy as np
  from ase.io.cube import read_cube_data
  from mayavi import mlab
  from pathlib import Path
  import os
  
  from collections import namedtuple
  
  ElementInfo = namedtuple('ElementInfo', ['radius', 'color'])
  element_data = {
      "H":    ElementInfo(radius=0.32, color=(1.00, 1.00, 1.00)),
      "C":    ElementInfo(radius=0.77, color=(0.56, 0.56, 0.56)),
      "O":    ElementInfo(radius=0.66, color=(1.00, 0.05, 0.05)),
      "F":    ElementInfo(radius=0.64, color=(0.56, 0.88, 0.31)),
      "Zn":   ElementInfo(radius=1.25, color=(0.49, 0.50, 0.69)),
      # 其他元素在这里添加
  }
  def get_element_info(symbol):
      symbol = symbol.capitalize()  # 确保符号首字母大写
      if symbol in element_data:
          return element_data[symbol]
      else:
          return ElementInfo(radius=None, color=None)  # 如果元素不存在，返回None值
  
  # 读取数据
  path_cal_res = os.path.dirname(os.path.abspath(__file__))
  dens, atoms = read_cube_data(path_cal_res / Path('density.cube'))
  
  
  print(dens.max())
  # 限制 dens 到 [0.001, 1.2]
  dens = np.clip(dens, 0.001, 1.9)
  
  
  
  # 获取网格坐标
  nx, ny, nz = dens.shape
  x = np.linspace(0, atoms.cell[0, 0], nx, endpoint=False)
  y = np.linspace(0, atoms.cell[1, 1], ny, endpoint=False)
  z = np.linspace(0, atoms.cell[2, 2], nz, endpoint=False)
  
  # 创建三维网格
  x_grid, y_grid, z_grid = np.meshgrid(x, y, z, indexing='ij')
  
  # 绘制等值面
  mlab.figure('Density Isosurface', bgcolor=(1, 1, 1), size=(800, 600))
  
  # 选择等值面
  contours = list(np.linspace(0.01,1.9,10))
  # contours = [0.001, 0.1, 0.5, 1]
  # contours = [0.001]
  contour = mlab.contour3d(x_grid, y_grid, z_grid, dens, contours=contours, colormap='viridis')
  
  
  # 设置透明度
  contour.actor.property.opacity = 0.08  # 透明度范围 0.0（完全透明）到 1.0（完全不透明）
  
  
  # 原子数据 xyzr
  atoms_space_info = atoms.get_positions()
  atoms_sysmbol = atoms.get_chemical_symbols()
  
  
  # Create spheres
  for i, atom in enumerate(atoms_space_info):
      x, y, z = atom
      element_info = get_element_info(atoms_sysmbol[i])
      # 默认使用element_info.radius/2的半径渲染原子，防止等值面看不清晰
      mlab.points3d(x, y, z, scale_factor=element_info.radius/2, color=element_info.color, resolution=32, opacity=0.9)
  
  
  # 显示'viridis'色带，对应数值dens的0.001到1.2
  colorbar = mlab.colorbar(contour, title='Density', orientation='vertical', nb_labels=13)
  
  # 设置刻度数字颜色为黑色
  colorbar.scalar_bar.unconstrained_font_size = True  # 确保字体大小可调
  colorbar.label_text_property.font_size = 10  # Increase font size
  colorbar.label_text_property.color = (0, 0, 0)  # RGB格式，黑色为 (0, 0, 0)
  colorbar.title_text_property.color = (0, 0, 0)  # Black color for the title
  colorbar.title_text_property.font_size = 14  # Increase font size
  
  
  # 显示等值面
  mlab.show()
  ```
### 静电势3d图像代码：
  ```
  import numpy as np
  from ase.io.cube import read_cube_data
  from mayavi import mlab
  from pathlib import Path
  import os
  import matplotlib.pyplot as plt
  from scipy.spatial import Delaunay
  
  from collections import namedtuple
  
  ElementInfo = namedtuple('ElementInfo', ['radius', 'color'])
  element_data = {
      "H":    ElementInfo(radius=0.32, color=(1.00, 1.00, 1.00)),
      "C":    ElementInfo(radius=0.77, color=(0.56, 0.56, 0.56)),
      "O":    ElementInfo(radius=0.66, color=(1.00, 0.05, 0.05)),
      "F":    ElementInfo(radius=0.64, color=(0.56, 0.88, 0.31)),
      "Zn":   ElementInfo(radius=1.25, color=(0.49, 0.50, 0.69)),
      # 其他元素在这里添加
  }
  def get_element_info(symbol):
      symbol = symbol.capitalize()  # 确保符号首字母大写
      if symbol in element_data:
          return element_data[symbol]
      else:
          return ElementInfo(radius=None, color=None)  # 如果元素不存在，返回None值
  
  
  # 读取数据
  path_cal_res = os.path.dirname(os.path.abspath(__file__))
  dens, atoms0 = read_cube_data(path_cal_res / Path('density.cube'))
  esps, atoms2 = read_cube_data(path_cal_res / Path('esp.cube'))
  print(esps.min(),esps.max())
  dens = np.clip(dens, 0.0001, 1)
  
  
  # 获取网格坐标
  nx, ny, nz = dens.shape
  x = np.linspace(0, atoms0.cell[0, 0], nx, endpoint=False)
  y = np.linspace(0, atoms0.cell[1, 1], ny, endpoint=False)
  z = np.linspace(0, atoms0.cell[2, 2], nz, endpoint=False)
  
  # 创建三维网格
  x_grid, y_grid, z_grid = np.meshgrid(x, y, z, indexing='ij')
  
  # 绘制等值面
  mlab.figure('Density Isosurface with ESP', bgcolor=(1, 1, 1), size=(800, 600))
  
  # 选择等值面
  contours = [0.001]
  contour = mlab.contour3d(x_grid, y_grid, z_grid, dens, contours=contours, colormap='viridis')
  
  
  
  # 设置透明度
  contour.actor.property.opacity = 0.0
  
  # 提取等值面点
  points = contour.actor.mapper.input.points.to_array()
  
  # 使用 ESP 数据重新映射颜色
  # First, we need to interpolate ESP values at the points of the isosurface
  from scipy.interpolate import RegularGridInterpolator
  
  # Create an interpolator for the ESP data
  esp_interpolator = RegularGridInterpolator((x, y, z), esps)
  
  # Interpolate ESP at each point of the isosurface
  esp_values = esp_interpolator(points)
  
  print(min(esp_values),max(esp_values))
  
  
  # Create Delaunay triangulation for surface
  tri = Delaunay(points)
  
  # Create a surface from the triangulation
  surface = mlab.pipeline.surface(mlab.pipeline.delaunay3d(mlab.pipeline.scalar_scatter(points[:,0], points[:,1], points[:,2], esp_values)))
  
  # Adjust surface properties
  if abs(max(esp_values) - min(esp_values)) >1:
    surface.module_manager.scalar_lut_manager.data_range = (min(esp_values), max(esp_values))
  elif abs(max(esp_values) - min(esp_values)) <=1:
    surface.module_manager.scalar_lut_manager.data_range = (esps.min(),esps.max())
  surface.module_manager.scalar_lut_manager.lut_mode = 'RdBu'
  surface.actor.property.opacity = 0.6
  
  # 显示色带
  surface_ = mlab.colorbar(surface, title='ESP', orientation='vertical', nb_labels=11)
  
  # Adjust colorbar text
  surface.module_manager.scalar_lut_manager.label_text_property.color = (0, 0, 0)
  surface_.label_text_property.font_size = 1
  surface.module_manager.scalar_lut_manager.title_text_property.color = (0, 0, 0)
  surface_.title_text_property.font_size = 2
  
  
  
  # 原子数据
  atoms_sysmbol = atoms0.get_chemical_symbols()
  
  
  # Create spheres
  for i, atom in enumerate(atoms0.get_positions()):
      x, y, z = atom
      element_info = get_element_info(atoms_sysmbol[i])
      # 默认使用element_info.radius/2的半径渲染原子，防止等值面看不清晰
      mlab.points3d(x, y, z, scale_factor=element_info.radius/2, color=element_info.color, resolution=32, opacity=0.9)
  mlab.show()
  ```

# 常用代码片段
## 使用机器学习势如MACE快速优化结构
  ```
  import os
  from pathlib import Path
  import random
  from ase.build import molecule
  import numpy as np
  
  from ase.visualize import view
  from ase.build import graphene
  from ase import Atoms
  from ase.optimize import BFGS
  from ase.io import write
  from ase.constraints import Hookean, FixAtoms
  
  X_Y_bond_length=1.51 # Å 1.3 1.5
  fmax=0.05
  file_name_prefix='g'
  other_info=fr'oh'
  
  mace_model_path_windows = r"d:\Downloads\mace-mpa-0-medium.model"
  mace_model_path_linux = r"mnt/d/Downloads/2023-12-03-mace-128-L1_epoch-199.model"
  
  path_cal_res = os.path.dirname(os.path.abspath(__file__))
  # path_cal_res=r'd:\UbuntuFiles\cal'
  environment = None
  
  
  if os.name == "nt":  # 检查是否是 Windows
      environment = 'Windows'
  elif os.name == "posix":  # 检查是否是类 Unix 系统
      system_name = platform.system()  # 获取系统名称
      if system_name == "Linux":  
          # 检查是否是 WSL 环境
          if "microsoft" in platform.uname().release.lower():
              environment = 'Ubuntu (WSL)'
              # path_cal_res = r'/mnt/d/UbuntuFiles/cal'
              path_cal_res = ''
          else:
              environment = 'Linux'
      elif system_name == "Darwin":
          environment = 'macOS'
      else:
          environment = "Unknown POSIX System"
  else:
      environment = "Unknown OS"
  
  
  def get_system():
      system = Atoms()
      return system
  
  
  def main(mode):
      """ 
      mode: 'mace' or 'chgnet'
      """
      system = get_system()
  
      print(system)
      view(system)
  
      # 保存未优化的建模文件
      write(path_cal_res /
          Path(rf'{file_name_prefix}_{other_info}_no_optim.traj'), system)
      calc = None
      if mode == 'mace':
          from mace.calculators import mace_mp
          if environment == 'Windows':
              # calc = mace_mp(model=r"d:\Downloads\2023-12-03-mace-128-L1_epoch-199.model")
              calc = mace_mp(model=mace_model_path_windows)
          elif environment == 'Ubuntu (WSL)':
              calc = mace_mp(model=mace_model_path_linux)
      elif mode == 'chgnet':
          from chgnet.model.model import CHGNet
          from chgnet.model.dynamics import CHGNetCalculator
          chgnet = CHGNet.load(model_name="0.3.0")  # 您可以更改 model_name 来选择不同的模型
          calc = CHGNetCalculator(chgnet)
  
      system.calc = calc
      # 优化几何结构
      dyn = BFGS(system, trajectory=str(path_cal_res /
                 Path(fr'{file_name_prefix}_{other_info}_optiming_{mode}.traj')))
      dyn.run(fmax=fmax)
      # 获取总能量
      e_system = system.get_potential_energy()
      print(f"Binding energy of {file_name_prefix}: {e_system} eV")
  
  main('mace')
  # main('chgnet')
  ```
## 使用GPAW进行二次优化
```

import os
from pathlib import Path
from ase.build import graphene
from ase import Atoms
from ase.optimize import BFGS
from gpaw import GPAW, PW
from ase.visualize import view
from ase.constraints import Hookean
from ase.constraints import Hookean, FixAtoms
# from ase.build.surface import add_adsorbate
# from ase.build import molecule

X_Y_bond_length = 1.5  # Å 1.3 1.5
cutoff_energy = 400   # 截止能量 eV
fmax = 0.01  # 最大受力的阈值 eV/Å
file_name_prefix = 'g'
other_info = fr'oh'

path_cal_res = os.path.dirname(os.path.abspath(__file__))


def main(mode='mace'):
    keeping('mace')



def keeping(mode):
    from ase.io import Trajectory
    # 设置 GPAW 计算器
    _ = path_cal_res / \
        Path(rf'{file_name_prefix}_{other_info}_output_keeping.txt')
    calc = GPAW(mode=PW(cutoff_energy), xc='PBE', charge =0.0 , kpts=(2,2,1), txt=str(_))
    # 加载轨迹文件并获取最后一步的结构
    traj = Trajectory(str(
        path_cal_res / Path(rf'{file_name_prefix}_{other_info}_optiming_{mode}.traj')), 'r')
    system = traj[-1]  # 获取轨迹文件中的最后一步结构

    # 对system的处理

    # 重新设置计算器
    system.calc = calc

    # 继续优化并保存到新的轨迹文件
    dyn = BFGS(system, trajectory=str(path_cal_res /
               Path(rf'{file_name_prefix}_{other_info}_optiming_keeping.traj')))
    dyn.run(fmax=fmax)

    # 获取总能量
    e_system = system.get_potential_energy()

    print(f"Binding energy of {file_name_prefix}: {e_system} eV")


main('mace')
```
## 电子密度
### 获取电子密度
```
import os
from ase.build import molecule
from ase.io import write
from ase.units import Bohr
from gpaw import GPAW, PW
from gpaw.analyse.hirshfeld import HirshfeldPartitioning
from ase.io import Trajectory
from pathlib import Path

# 选择轨迹名称和index
traj_name = 'g_C-COOH_constraint_optiming_keeping.traj'
traj_index = 13


path_big_file = '/mnt/d/cal' # 由于cube文件较大，可以设置额外的目录
path_cal_res = os.path.dirname(os.path.abspath(__file__))
# 有时一个目录下可能有许多traj文件都需要电子密度 则os.path.basename(path_cal_res)+flag
new_dir_path = os.path.join(path_big_file, os.path.basename(path_cal_res)+'onlyg-OH')
os.makedirs(new_dir_path, exist_ok=True)

traj = Trajectory(path_cal_res / Path(traj_name), 'r') 
system = traj[traj_index]

calc = GPAW(mode=PW(400), xc='PBE', h=0.2, charge=0.0, kpts=(2,2,1), txt=new_dir_path / Path('system_output.txt'))
# 'fd' 模式下 设置pcb false 对于原本是orthogonal的系统会有问题，而pw模式下无法对电荷求解，各有劣势
# 由于一般是pw模式，所以需要注释掉电荷相关
system.calc = calc
print(system)
print(system.get_potential_energy())
# calc.write(new_dir_path / Path('system1_fd.paw'), mode='all')


# write Hirshfeld charges out
# hf = HirshfeldPartitioning(system.calc)
# for atom, charge in zip(system, hf.get_charges()):
#     atom.charge = charge
# # system.write('Hirshfeld.traj') #  Trajectory writer needs a fix
system.copy().write(new_dir_path / Path('system.traj'))

# create electron density cube file ready for bader
rho = system.calc.get_all_electron_density(gridrefinement=2)
write(new_dir_path / Path('den.cube'), system, data=rho * Bohr**3)
```
### 显示3D电子密度
```
# fmt: off

import optparse

import numpy as np

from ase.calculators.calculator import get_calculator_class
from ase.data import covalent_radii
from ase.data.colors import jmol_colors
from ase.io.cube import read_cube_data


cpk_colors = jmol_colors

atoms_magnification = 1.8

def plot(atoms, data, contours):
    """Plot atoms, unit-cell and iso-surfaces using Mayavi.

    Parameters:

    atoms: Atoms object
        Positions, atomiz numbers and unit-cell.
    data: 3-d ndarray of float
        Data for iso-surfaces.
    countours: list of float
        Contour values.
    """

    # Delay slow imports:
    import os

    from mayavi import mlab

    # mayavi GUI bug fix for remote access via ssh (X11 forwarding)
    if "SSH_CONNECTION" in os.environ:
        f = mlab.gcf()
        f.scene._lift()

    mlab.figure(1, bgcolor=(1, 1, 1))  # make a white figure

    # Plot the atoms as spheres:
    for pos, Z in zip(atoms.positions, atoms.numbers):
        mlab.points3d(*pos,
                      scale_factor=covalent_radii[Z]*atoms_magnification,
                      resolution=20,
                      color=tuple(cpk_colors[Z]))

    # Draw the unit cell:
    A = atoms.cell
    # for i1, a in enumerate(A):
    #     i2 = (i1 + 1) % 3
    #     i3 = (i1 + 2) % 3
    #     for b in [np.zeros(3), A[i2]]:
    #         for c in [np.zeros(3), A[i3]]:
    #             p1 = b + c
    #             p2 = p1 + a
    #             mlab.plot3d([p1[0], p2[0]],
    #                         [p1[1], p2[1]],
    #                         [p1[2], p2[2]],
    #                         tube_radius=0.1)

    cp = mlab.contour3d(data, contours=contours, transparent=True,
                        opacity=0.1, colormap='viridis')
    cp.module_manager.scalar_lut_manager.data_range = [0, 1]

    # Do some tvtk magic in order to allow for non-orthogonal unit cells:
    polydata = cp.actor.actors[0].mapper.input
    pts = np.array(polydata.points) - 1
    # Transform the points to the unit cell:
    polydata.points = np.dot(pts, A / np.array(data.shape)[:, np.newaxis])
    colorbar = mlab.colorbar(cp, title='Density', orientation='vertical', nb_labels=11)
    # 设置刻度数字颜色为黑色
    colorbar.scalar_bar.unconstrained_font_size = True  # 确保字体大小可调
    colorbar.label_text_property.font_size = 10  # Increase font size
    colorbar.label_text_property.color = (0, 0, 0)  # RGB格式，黑色为 (0, 0, 0)
    colorbar.title_text_property.color = (0, 0, 0)  # Black color for the title
    colorbar.title_text_property.font_size = 14  # Increase font size

    # 添加坐标轴
    axes = mlab.axes(cp, color=(0, 0, 0), xlabel='X', ylabel='Y', zlabel='Z')

    # 自定义刻度
    axes.label_text_property.font_size = 12
    axes.label_text_property.color = (0, 0, 0)  # Black labels
    axes.title_text_property.color = (0, 0, 0)  # Black axis titles
    axes.title_text_property.font_size = 14
    # Apparently we need this to redraw the figure, maybe it can be done in
    # another way?
    mlab.view(azimuth=155, elevation=70, distance='auto')
    # Show the 3d plot:
    mlab.show()


from pathlib import Path
import os

# 读取数据
path_cal_res = os.path.dirname(os.path.abspath(__file__))
dens, atoms = read_cube_data(path_cal_res / Path('den.cube'))

# 限制 dens 到 [0.001, 1.2]
dens = np.clip(dens, 0.0001, 1.0)

contours = [0.001]
# contours = list(np.linspace(0.1,1,10))
plot(atoms,dens,contours)
```
## 静电势
### 获取静电势
```
from ase.build import molecule
from ase.io import write
from ase.units import Bohr
from gpaw import GPAW,PW
from ase.io import Trajectory

from pathlib import Path
import os


# 选择轨迹名称和index
traj_name = 'g_C-COOH_constraint_optiming_keeping.traj'
traj_index = 13


path_big_file = '/mnt/d/cal'
path_cal_res = os.path.dirname(os.path.abspath(__file__))
# 有时一个目录下可能有许多traj文件都需要电子密度 则os.path.basename(path_cal_res)+flag
new_dir_path = os.path.join(path_big_file, os.path.basename(path_cal_res)+'onlyg-OH')
os.makedirs(new_dir_path, exist_ok=True)

traj = Trajectory(path_cal_res / Path(traj_name), 'r') 
system = traj[traj_index]

system.calc = GPAW(xc='PBE',
                 mode=PW(400),
                 charge=0.0,
                 h=0.2,
                 kpts=(2,2,1),
                 txt=new_dir_path / Path('pw.txt'))
print(system)
print(system.get_potential_energy())
system.calc.write(new_dir_path / Path('pw.gpw'))
v = system.calc.get_electrostatic_potential()
write(new_dir_path / Path('esp.cube'), system, data=v)
```
### 显示某个电子密度等值面上的静电势
```
import numpy as np
from ase.io.cube import read_cube_data
from pathlib import Path
import os
from mayavi import mlab
from ase.data import covalent_radii
from ase.data.colors import jmol_colors
cpk_colors = jmol_colors
atoms_magnification = 1.8

# Function to extract points and ESP values on a specific density isosurface
def extract_points_and_esps_on_isosurface(data, esps, iso_value):
    """
    Extract points and ESP values on a specific density isosurface.

    Parameters:
        data: 3D ndarray
            The density data.
        esps: 3D ndarray
            The electrostatic potential data.
        iso_value: float
            The density isosurface value.

    Returns:
        points: ndarray
            The coordinates of the points on the isosurface.
        esp_values: ndarray
            The ESP values at the points on the isosurface.
    """
    # Find points on the isosurface (allow for a small tolerance)
    tolerance = 1e-4
    tolerance = 10**(-3.5)
    # print(10**(-4),10**(-3.6),10**(-3.5),10**(-3.4),10**(-3.3))
    isosurface_indices = np.where((data >= iso_value - tolerance) & (data <= iso_value + tolerance))

    # Extract coordinates of the points
    points = np.column_stack(isosurface_indices)

    # # Filter points based on z-coordinate
    # valid_indices = points[:, 2] <= 200
    # points = points[valid_indices]

    # # Update isosurface_indices accordingly
    # isosurface_indices = tuple(index[valid_indices] for index in isosurface_indices)

    # Extract corresponding ESP values
    esp_values = esps[isosurface_indices]

    return points, esp_values
def plot_with_surface(atoms, data, points, esp_values, iso_value):
    """
    Plot atoms, unit-cell, iso-surfaces, and ESP surface using Mayavi.

    Parameters:
        atoms: Atoms object
            Positions, atomic numbers, and unit-cell.
        data: 3D ndarray of float
            Data for iso-surfaces.
        points: ndarray
            Coordinates of the points on the isosurface.
        esp_values: ndarray
            ESP values at the points on the isosurface.
        iso_value: float
            The density isosurface value.
    """
    mlab.figure(1, bgcolor=(1, 1, 1))  # Make a white figure
    # Plot the atoms as spheres:
    for pos, Z in zip(atoms.positions, atoms.numbers):
        mlab.points3d(*pos,
                      scale_factor=covalent_radii[Z]*atoms_magnification,
                      resolution=20,
                      color=tuple(cpk_colors[Z]))
    # Plot the density isosurface
    cp_density = mlab.contour3d(data, contours=[iso_value], transparent=True,
                                opacity=0.0, colormap='viridis')
    A = atoms.cell
    transformed_points = np.dot(points - 1, A / np.array(data.shape)[:, np.newaxis])

    # Use transformed points in mlab.points3d
    cp_surface = mlab.points3d(transformed_points[:, 0], transformed_points[:, 1], transformed_points[:, 2],
                               esp_values,
                               scale_mode='none', 
                               scale_factor=0.3,
                               colormap='RdBu', 
                               mode="sphere",
                               resolution=12,
                               opacity=0.03)

    cp_surface.module_manager.scalar_lut_manager.data_range = [esp_values.min(), esp_values.max()]
    print(esp_values.min(),esp_values.max())
    colorbar = mlab.colorbar(cp_surface, title='esp', orientation='vertical', nb_labels=2)
    # 设置刻度数字颜色为黑色
    colorbar.scalar_bar.unconstrained_font_size = True  # 确保字体大小可调
    colorbar.label_text_property.font_size = 24  # Increase font size
    colorbar.label_text_property.color = (0, 0, 0)  # RGB格式，黑色为 (0, 0, 0)
    colorbar.title_text_property.color = (0, 0, 0)  # Black color for the title
    colorbar.title_text_property.font_size = 0  # Increase font size

    # Combine cp_density and cp_surface for tvtk processing
    cp = [cp_density]

    # Do some tvtk magic in order to allow for non-orthogonal unit cells
    A = atoms.cell
    for c in cp:
        polydata = c.actor.actors[0].mapper.input
        pts = np.array(polydata.points) - 1
        polydata.points = np.dot(pts, A / np.array(data.shape)[:, np.newaxis])

    # Show the 3D plot
    mlab.show()


def main():
    # Define the path to the current script directory
    path_cal_res = Path(os.path.dirname(os.path.abspath(__file__)))

    # Read density and ESP data
    dens, atoms = read_cube_data(path_cal_res / 'den.cube')
    esps, _ = read_cube_data(path_cal_res / 'esp.cube')
    esps = esps * -1
    # Extract points and ESP values on the density isosurface (den = 0.001)
    iso_value = 0.001
    points, esp_values = extract_points_and_esps_on_isosurface(dens, esps, iso_value)

    # Print the points and ESP values
    print(f"Points on the density isosurface (den = {iso_value}):")
    print(points)
    print(f"ESP values on the density isosurface (den = {iso_value}):")
    print(esp_values)
    print(len(points))
    print(len(esp_values))

    # Plot with surface
    plot_with_surface(atoms, dens, points, esp_values, iso_value)

if __name__ == "__main__":
    main()
```
## 电子密度差（电荷差）
### 吸附/反应后的电子密度
```
import os
from ase.build import molecule
from ase.io import write
from ase.units import Bohr
from gpaw import GPAW, PW

from ase.io import Trajectory
from pathlib import Path
from ase.visualize import view

# 选择轨迹名称和index
traj_name = 'g_C-COOH_constraint_zn2+_constraint_38-39-4_36-39-5_optiming_keeping.traj'
traj_index = 84
# 差分电荷 + flag 用于区分den.cube文件
flag = 'whole'

path_big_file = '/mnt/d/cal'
path_cal_res = os.path.dirname(os.path.abspath(__file__))
new_dir_path = os.path.join(path_big_file, os.path.basename(path_cal_res)+'-ccd')
os.makedirs(new_dir_path, exist_ok=True)

traj = Trajectory(path_cal_res / Path(traj_name), 'r') 
system = traj[traj_index]
# del system[[atom.index for atom in system if atom.symbol=='H' or atom.symbol=='O']] #删除所有指定原子
# view(system)

# calc = GPAW(mode=PW(400), xc='PBE', h=0.2, charge=+2.0, kpts=(2,2,1), txt=path_cal_res / Path('system1_output1.txt'))
calc = GPAW(mode=PW(400), xc='PBE', h=0.2, charge=+2.0, kpts=(2,2,1), txt=new_dir_path / Path(f'system-{flag}.txt'))

system.calc = calc
print(system.calc)
print(system)
print(system.get_potential_energy())
# calc.write(new_dir_path / Path('system1_fd.paw'), mode='all')

system.copy().write(new_dir_path / Path('system.traj'))

# create electron density cube file ready for bader
rho = system.calc.get_all_electron_density(gridrefinement=2)
write(new_dir_path / Path(f'den-{flag}.cube'), system, data=rho * Bohr**3)
```
### 载体/平面物
```
import os
from ase.build import molecule
from ase.io import write
from ase.units import Bohr
from ase.constraints import Hookean, FixAtoms
from gpaw import GPAW, PW

from ase.io import Trajectory
from pathlib import Path
from ase.visualize import view

# 选择轨迹名称和index
traj_name = 'g_C-COOH_constraint_zn2+_constraint_38-39-4_36-39-5_optiming_keeping.traj'
traj_index = 84
# 差分电荷 + flag 用于区分den.cube文件
flag = 'A'

path_big_file = '/mnt/d/cal'
path_cal_res = os.path.dirname(os.path.abspath(__file__))
new_dir_path = os.path.join(path_big_file, os.path.basename(path_cal_res)+'-ccd')
os.makedirs(new_dir_path, exist_ok=True)

traj = Trajectory(path_cal_res / Path(traj_name), 'r') 
system = traj[traj_index]
del system.constraints
k=15
layer_distance_constraint = [
                            Hookean(a1=17, a2=36, rt=1.3825931989763465, k=k),
                            Hookean(a1=17, a2=38, rt=1.3871286855967784, k=k),
                             FixAtoms([atom.index for atom in system if atom.symbol == 'C' and atom.index not in [17,16,11]])
                             ]
del system[[atom.index for atom in system if atom.symbol=='Zn']] #删除所有指定原子
system.set_constraint(layer_distance_constraint)
# view(system)
print(system)

# calc = GPAW(mode=PW(400), xc='PBE', h=0.2, charge=+2.0, kpts=(2,2,1), txt=path_cal_res / Path('system1_output1.txt'))
calc = GPAW(mode=PW(400), xc='PBE', h=0.2, charge=+0.0, kpts=(2,2,1), txt=new_dir_path / Path(f'system-{flag}.txt'))

system.calc = calc
print(system.calc)
print(system)
print(system.get_potential_energy())
# calc.write(new_dir_path / Path('system1_fd.paw'), mode='all')

system.copy().write(new_dir_path / Path('system.traj'))

# create electron density cube file ready for bader
rho = system.calc.get_all_electron_density(gridrefinement=2)
write(new_dir_path / Path(f'den-{flag}.cube'), system, data=rho * Bohr**3)
```
### 被吸附物的电子密度
```
import os
from ase.build import molecule
from ase.io import write
from ase.units import Bohr
from gpaw import GPAW, PW

from ase.io import Trajectory
from pathlib import Path
from ase.visualize import view

# 选择轨迹名称和index
traj_name = 'g_C-COOH_constraint_zn2+_constraint_38-39-4_36-39-5_optiming_keeping.traj'
traj_index = 84
# 差分电荷 + flag 用于区分den.cube文件
flag = 'B'

path_big_file = '/mnt/d/cal'
path_cal_res = os.path.dirname(os.path.abspath(__file__))
new_dir_path = os.path.join(path_big_file, os.path.basename(path_cal_res)+'-ccd')
os.makedirs(new_dir_path, exist_ok=True)

traj = Trajectory(path_cal_res / Path(traj_name), 'r') 
system = traj[traj_index]
del system.constraints
del system[[atom.index for atom in system if atom.symbol!='Zn']] #删除所有指定原子
# view(system)
print(system)

# calc = GPAW(mode=PW(400), xc='PBE', h=0.2, charge=+2.0, kpts=(2,2,1), txt=path_cal_res / Path('system1_output1.txt'))
calc = GPAW(mode=PW(400), xc='PBE', h=0.2, charge=+2.0, kpts=(2,2,1), txt=new_dir_path / Path(f'system-{flag}.txt'))

system.calc = calc
print(system.calc)
print(system)
print(system.get_potential_energy())
# calc.write(new_dir_path / Path('system1_fd.paw'), mode='all')

system.copy().write(new_dir_path / Path('system.traj'))

# create electron density cube file ready for bader
rho = system.calc.get_all_electron_density(gridrefinement=2)
write(new_dir_path / Path(f'den-{flag}.cube'), system, data=rho * Bohr**3)
```
### 显示3D电子密度差和电荷位移曲线
```
# fmt: off

import optparse

import numpy as np

from ase.calculators.calculator import get_calculator_class
from ase.data import covalent_radii
from ase.data.colors import jmol_colors
from ase.io.cube import read_cube_data

cpk_colors = jmol_colors
atoms_magnification = 1

save_path = r'C:\Users\Administrator\Desktop' 
middle = 200
axis = 'z'
是否为石墨烯体系=True
# need_plot_3d=True
need_plot_3d=False

def plot(atoms, data, contours,axis=None):
    """Plot atoms, unit-cell and iso-surfaces using Mayavi.

    Parameters:

    atoms: Atoms object
        Positions, atomiz numbers and unit-cell.
    data: 3-d ndarray of float
        Data for iso-surfaces.
    countours: list of float
        Contour values.
    """

    # Delay slow imports:
    import os

    from mayavi import mlab

    # mayavi GUI bug fix for remote access via ssh (X11 forwarding)
    if "SSH_CONNECTION" in os.environ:
        f = mlab.gcf()
        f.scene._lift()

    mlab.figure(1, bgcolor=(1, 1, 1))  # make a white figure

    # Plot the atoms as spheres:
    for pos, Z in zip(atoms.positions, atoms.numbers):
        mlab.points3d(*pos,
                      scale_factor=covalent_radii[Z]*atoms_magnification,
                      resolution=20,
                      color=tuple(cpk_colors[Z]))

    # Draw the unit cell:
    A = atoms.cell
    for i1, a in enumerate(A):
        i2 = (i1 + 1) % 3
        i3 = (i1 + 2) % 3
        for b in [np.zeros(3), A[i2]]:
            for c in [np.zeros(3), A[i3]]:
                p1 = b + c
                p2 = p1 + a
                mlab.plot3d([p1[0], p2[0]],
                            [p1[1], p2[1]],
                            [p1[2], p2[2]],
                            tube_radius=0.1,
                            opacity=0.4)
    cp = mlab.contour3d(data, contours=contours, transparent=True,
                        opacity=0.4, colormap='coolwarm')

    # Do some tvtk magic in order to allow for non-orthogonal unit cells:
    polydata = cp.actor.actors[0].mapper.input
    pts = np.array(polydata.points) - 1
    # Transform the points to the unit cell:
    polydata.points = np.dot(pts, A / np.array(data.shape)[:, np.newaxis])
    cp.module_manager.scalar_lut_manager.data_range = [data.min(), data.max()]
    # colorbar = mlab.colorbar(cp, title='Density', orientation='vertical', nb_labels=21)
    colorbar = mlab.colorbar(cp, title='Density', orientation='vertical', nb_labels=2)
    # 设置刻度数字颜色为黑色
    colorbar.scalar_bar.unconstrained_font_size = True  # 确保字体大小可调
    colorbar.label_text_property.font_size = 24  # Increase font size
    colorbar.label_text_property.color = (0, 0, 0)  # RGB格式，黑色为 (0, 0, 0)
    colorbar.title_text_property.color = (0, 0, 0)  # Black color for the title
    colorbar.title_text_property.font_size = 0  # Increase font size

    # Add the lines at axis=10 and axis=30:
    if axis != None:
        axis_mapping = {'x': 0, 'y': 1, 'z': 2}
        axis_index=[i for i in range(3) if i != axis_mapping[axis]]
        x, y = np.meshgrid(np.linspace(0, A[axis_index[0], axis_index[0]], data.shape[axis_index[0]]), np.linspace(0, A[axis_index[1], axis_index[1]], data.shape[axis_index[1]]))

        # Create the lines for the planes at Z=100 and Z=300
        z1 = 10 * np.ones_like(x)  # Z=100 plane
        z2 = 30 * np.ones_like(x)  # Z=300 plane

        # Plot lines at Z=100 and Z=300
        mlab.plot3d(x.flatten(), y.flatten(), z1.flatten(), color=(0, 0, 0), tube_radius=0.05)  # Line at Z=100
        mlab.plot3d(x.flatten(), y.flatten(), z2.flatten(), color=(0, 0, 0), tube_radius=0.05)  # Line at Z=300

    

    mlab.view(azimuth=155, elevation=70, distance='auto')
    # Show the 3d plot:
    mlab.show()


from pathlib import Path
import os

# 读取数据
path_cal_res = os.path.dirname(os.path.abspath(__file__))
dens_whole, atoms_whole = read_cube_data(path_cal_res / Path('den-whole.cube'))
dens_A, atoms_A = read_cube_data(path_cal_res / Path('den-A.cube'))
dens_B, atoms_B = read_cube_data(path_cal_res / Path('den-B.cube'))

# 计算差分电荷密度
dens_diff = dens_whole - (dens_A + dens_B) 
print(dens_diff.min(),dens_diff.max(),dens_diff.mean())

# cbb_range_list=[0.0002, 0.002, 0.02]
cbb_range_list=[0.002]

# 绘制差分电荷密度图
for e in cbb_range_list:
    # 限制 dens_diff 以避免过小或过大的数值
    dens_diff_ = np.clip(dens_diff, -e, e)
    # 选择轮廓线值
    # contours = [-e/3*2, -e/3, e/3, e/3*2]  # 或者使用 np.linspace 生成一系列轮廓值
    contours = [-e/2, e/2]  # 或者使用 np.linspace 生成一系列轮廓值
    if 是否为石墨烯体系 and axis == 'z' and need_plot_3d:
        plot(atoms_whole, dens_diff_, contours, 'z')
    elif need_plot_3d:
        plot(atoms_whole, dens_diff_, contours)




""" 一系列曲线 """

def local_integral_curve(dens_diff, axis='z'):
    """
    计算密度差的局部积分曲线（沿指定轴的积分）

    参数：
    dens_diff : ndarray
        电荷密度差数据
    axis : str, 可选, 默认 'z'
        沿哪个方向计算局部积分，支持 'x', 'y', 'z' 方向
    """
    axis_mapping = {'x': 0, 'y': 1, 'z': 2}  # 映射轴名称到数组轴索引
    if axis not in axis_mapping:
        raise ValueError("Axis must be 'x', 'y', or 'z'.")
    
    # 对除指定轴外的其他轴求和，然后对指定轴进行积分
    axis_index = axis_mapping[axis]
    integral_curve = np.sum(dens_diff, axis=tuple(i for i in range(3) if i != axis_index))  # 除去指定轴的其他轴求和
    return integral_curve

def charge_displacement_curve(dens_diff, axis='z', z_ini=0):
    """
    计算电荷位移曲线，即对局部积分曲线进行进一步积分

    参数：
    dens_diff : ndarray
        电荷密度差数据
    axis : str, 可选, 默认 'z'
        沿哪个方向计算电荷位移曲线，支持 'x', 'y', 'z' 方向
    z_ini : int, 可选, 默认 0
        积分起点，表示从该 Z 层次开始积分
    """
    axis_mapping = {'x': 0, 'y': 1, 'z': 2}  # 映射轴名称到数组轴索引
    if axis not in axis_mapping:
        raise ValueError("Axis must be 'x', 'y', or 'z'.")
    
    # 对除指定轴外的其他轴求和，然后对指定轴进行积分
    axis_index = axis_mapping[axis]

    # 计算局部积分曲线
    local_integral = local_integral_curve(dens_diff, axis)

    # 进行积分
    # 假设我们对局部积分曲线进行数值积分
    z_values = np.arange(dens_diff.shape[axis_index])  # Z 轴上的网格点
    z_values = z_values[z_values >= z_ini]  # 只考虑从 z_ini 开始的积分部分
    
    # 对局部积分曲线进行累积积分，得到电荷位移曲线
    displacement_curve = np.cumsum(local_integral[z_values])  # 累积分
    

    return z_values, displacement_curve

def plane_averaged_curve(dens_diff, axis='z', 是否是石墨烯体系=False):
    """
    计算平面平均密度差曲线（每个截面的平均密度差）
    
    参数：
    dens_diff : ndarray
        电荷密度差数据
    axis : str, 可选, 默认 'z'
        沿哪个方向计算平面平均，支持 'x', 'y', 'z' 方向
    """
    axis_mapping = {'x': 0, 'y': 1, 'z': 2}  # 映射轴名称到数组轴索引
    if axis not in axis_mapping:
        raise ValueError("Axis must be 'x', 'y', or 'z'.")
    
    # 对除指定轴外的其他轴求和
    axis_index = axis_mapping[axis]
    summed_density = np.sum(dens_diff, axis=tuple(i for i in range(3) if i != axis_index))  # 求和
    
    axis_index = [i for i in range(3) if i != axis_index]
    # 计算每个截面的平均密度差
    if 是否是石墨烯体系:
        area_xy = dens_diff.shape[axis_index[0]] * dens_diff.shape[axis_index[1]] / 2 * 3**(1/2)
    else:
        area_xy = dens_diff.shape[axis_index[0]] * dens_diff.shape[axis_index[1]]  # XY 平面的面积
    plane_avg_density = summed_density / area_xy  # 除以XY平面的面积，得到平均值
    return plane_avg_density


def plot_curve1(dens_diff,axis='z'):
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter, MultipleLocator
    plt.figure(figsize=(2, 5))  # 设置图像长宽比为1:7
    axis_mapping = {'x': 0, 'y': 1, 'z': 2}  # 映射轴名称到数组轴索引
    # 计算电荷位移曲线
    # axis_values, displacement_curve = charge_displacement_curve(dens_diff, axis, 0)

    # 绘制局部积分曲线
    axis_values = np.arange(dens_diff.shape[axis_mapping[axis]])  # 获取指定轴的网格点
    local_integral = local_integral_curve(dens_diff, axis)
    print(local_integral.min(),local_integral.max(),local_integral.mean())
    # 绘制电荷位移曲线
    plt.plot(local_integral, axis_values)
    # plt.ylabel(f'{axis.upper()}-axis (Grid Points)')
    # plt.xlabel('Charge Displacement')
    # plt.title(f'Charge Displacement along {axis.upper()}-axis')
    # plt.grid(True)
    plt.ylim(axis_values.min()-1,axis_values.max()+1)
    plt.gca().invert_yaxis()
    if  是否为石墨烯体系:
        plt.ylim(301,99)
    # 格式化函数
    def divide_by_10(x, pos):
        if x==middle:
            return '0'
        return f'{-(x-middle)/10:.1f}'  # 显示1位小数，可改成.0f（无小数）或.2f（2位小数）
    plt.gca().yaxis.set_major_formatter(FuncFormatter(divide_by_10))
    plt.tight_layout()
    plt.savefig(os.path.join(save_path, 'ccd1'), dpi=500)
    plt.show()

def plot_curve2(dens_diff,axis='z'):
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter, MultipleLocator
    plt.figure(figsize=(2, 5))  # 设置图像长宽比为1:7
    axis_mapping = {'x': 0, 'y': 1, 'z': 2}  # 映射轴名称到数组轴索引
    # 计算电荷位移曲线
    axis_values, displacement_curve = charge_displacement_curve(dens_diff, axis, 0)

    print(displacement_curve.min(),displacement_curve.max(),displacement_curve.mean())
    # 绘制电荷位移曲线
    plt.plot(displacement_curve, axis_values)
    # plt.ylabel(f'{axis.upper()}-axis (Grid Points)')
    # plt.xlabel('Charge Displacement')
    # plt.title(f'Charge Displacement along {axis.upper()}-axis')
    # plt.grid(True)
    plt.ylim(axis_values.min()-1,axis_values.max()+1)
    plt.gca().invert_yaxis()
    if  是否为石墨烯体系:
        plt.ylim(301,99)
    # 格式化函数
    def divide_by_10(x, pos):
        if x==middle:
            return '0'
        return f'{-(x-middle)/10:.1f}'  # 显示1位小数，可改成.0f（无小数）或.2f（2位小数）
    plt.gca().yaxis.set_major_formatter(FuncFormatter(divide_by_10))
    plt.tight_layout()
    plt.savefig(os.path.join(save_path, 'ccd2'), dpi=500)
    plt.show()

def plot_curve3(dens_diff,axis='z',是否是石墨烯体系=False):
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter, MultipleLocator
    plt.figure(figsize=(2, 5))  # 设置图像长宽比为1:7
    axis_mapping = {'x': 0, 'y': 1, 'z': 2}  # 映射轴名称到数组轴索引
    # 计算电荷位移曲线
    # axis_values, displacement_curve = charge_displacement_curve(dens_diff, axis, 0)

    axis_values = np.arange(dens_diff.shape[axis_mapping[axis]])  # 获取指定轴的网格点
    # 计算平面平均密度差曲线
    plane_avg_density = plane_averaged_curve(dens_diff, axis)
    print(plane_avg_density.min(),plane_avg_density.max(),plane_avg_density.mean())
    # 绘制电荷位移曲线
    plt.plot(plane_avg_density, axis_values)
    # plt.ylabel(f'{axis.upper()}-axis (Grid Points)')
    # plt.xlabel('Charge Displacement')
    # plt.title(f'Charge Displacement along {axis.upper()}-axis')
    # plt.grid(True)
    plt.ylim(axis_values.min()-1,axis_values.max()+1)
    plt.gca().invert_yaxis()
    if  是否为石墨烯体系:
        plt.ylim(301,99)
    # 格式化函数
    def divide_by_10(x, pos):
        if x==middle:
            return '0'
        return f'{-(x-middle)/10:.1f}'  # 显示1位小数，可改成.0f（无小数）或.2f（2位小数）
    plt.gca().yaxis.set_major_formatter(FuncFormatter(divide_by_10))
    plt.tight_layout()
    plt.savefig(os.path.join(save_path, 'ccd3'), dpi=500)
    plt.show()

plot_curve1(dens_diff,axis)
plot_curve2(dens_diff,axis)
plot_curve3(dens_diff,axis,是否为石墨烯体系)
```
## HOMO LUMO
### 获取gpw（方便多次获取HOMO信息）
```
import os
from ase.units import Bohr
import numpy as np
from gpaw import GPAW, PW, restart
from ase.io import Trajectory
from pathlib import Path
from ase.io import write

# 选择轨迹名称和index
traj_name = 'g_C-F_optiming_keeping.traj'
traj_index = 8

path_big_file = '/mnt/d/cal'
path_cal_res = os.path.dirname(os.path.abspath(__file__))
# 有时一个目录下可能有许多traj文件都需要 则os.path.basename(path_cal_res)+flag
new_dir_path = os.path.join(path_big_file, os.path.basename(path_cal_res)+'onlyg-C-F')
os.makedirs(new_dir_path, exist_ok=True)



traj = Trajectory(path_cal_res / Path(traj_name), 'r') 
system = traj[traj_index]

calc = GPAW(mode=PW(400), xc='PBE', h=0.2, charge=0.0, kpts=(2,2,1), txt=None)
system.calc = calc
print(system)
print(system.get_potential_energy())
system.calc.write(new_dir_path / Path('pw_all.gpw'), mode='all')
```
### 获取homo lomo波函数
```
import os
from ase.units import Bohr
import numpy as np
from gpaw import GPAW, PW, restart
from ase.io import Trajectory
from pathlib import Path
from ase.io import write


path_big_file = '/mnt/d/cal'
path_cal_res = os.path.dirname(os.path.abspath(__file__))
# 有时一个目录下可能有许多traj文件都需要 则os.path.basename(path_cal_res)+flag
new_dir_path = os.path.join(path_big_file, os.path.basename(path_cal_res)+'onlyg-OH')
os.makedirs(new_dir_path, exist_ok=True)


atoms, calc = restart(new_dir_path / Path('pw_all.gpw'))
# print(calc.wfs.kpt_qs)
# print(calc.wfs.kpt_qs[0][0].psit.array.shape)
# print(len(calc.wfs.kpt_qs))

# 获取轨道和能量
eigenvalues = calc.get_eigenvalues()  # 获取所有轨道的能量
homo,lumo = atoms.calc.get_homo_lumo()
print(eigenvalues)
print(len(eigenvalues))
print('nband',calc.get_number_of_bands())
print('occ occupation_numbers')
print(calc.get_occupation_numbers())
print('homo 一般情况: ',(calc.get_occupation_numbers() >= 1.0).sum() - 1)
print(homo,lumo)
homo_band = np.abs(eigenvalues - homo).argmin()
lumo_band = np.abs(eigenvalues - lumo).argmin()
print(homo_band,lumo_band)
if lumo_band == homo_band and lumo_band != len(eigenvalues)-1 :
    lumo_band = homo_band + 1
print(homo_band,lumo_band)
# 提取HOMO和LUMO波函数
homo_wave_function = calc.get_pseudo_wave_function(band=homo_band, kpt = 1, periodic=True)
lumo_wave_function = calc.get_pseudo_wave_function(band=lumo_band, kpt = 1, periodic=True)

# print('homo_wave_function: ', homo_wave_function)
# print('lumo_wave_function: ', lumo_wave_function)
# write(path_cal_res / Path('homo_wave_function.cube'), atoms, data=homo_wave_function * Bohr**1.5)
# write(path_cal_res / Path('lumo_wave_function.cube'), atoms, data=lumo_wave_function * Bohr**1.5)
```
### 显示能级
```
import optparse
import numpy as np
from ase.calculators.calculator import get_calculator_class
from ase.data import covalent_radii
from ase.data.colors import cpk_colors
from ase.io.cube import read_cube_data

atoms_magnification = 1

def plot(atoms, data, contours):
    """Plot atoms, unit-cell and iso-surfaces using Mayavi.

    Parameters:

    atoms: Atoms object
        Positions, atomiz numbers and unit-cell.
    data: 3-d ndarray of float
        Data for iso-surfaces.
    countours: list of float
        Contour values.
    """

    # Delay slow imports:
    import os

    from mayavi import mlab


    mlab.figure(1, bgcolor=(1, 1, 1))  # make a white figure

    # Plot the atoms as spheres:
    for pos, Z in zip(atoms.positions, atoms.numbers):
        mlab.points3d(*pos,
                      scale_factor=covalent_radii[Z]*atoms_magnification,
                      resolution=20,
                      color=tuple(cpk_colors[Z]))

    # Draw the unit cell:
    A = atoms.cell
    # for i1, a in enumerate(A):
    #     i2 = (i1 + 1) % 3
    #     i3 = (i1 + 2) % 3
    #     for b in [np.zeros(3), A[i2]]:
    #         for c in [np.zeros(3), A[i3]]:
    #             p1 = b + c
    #             p2 = p1 + a
    #             mlab.plot3d([p1[0], p2[0]],
    #                         [p1[1], p2[1]],
    #                         [p1[2], p2[2]],
    #                         tube_radius=0.1)

    cp = mlab.contour3d(data, contours=contours, transparent=True,
                        opacity=0.7, colormap='viridis')
    # cp.module_manager.scalar_lut_manager.data_range = [0,0.002]

    # Do some tvtk magic in order to allow for non-orthogonal unit cells:
    polydata = cp.actor.actors[0].mapper.input
    pts = np.array(polydata.points) - 1
    # Transform the points to the unit cell:
    polydata.points = np.dot(pts, A / np.array(data.shape)[:, np.newaxis])
    colorbar = mlab.colorbar(cp, title='Density', orientation='vertical', nb_labels=11)
    # 设置刻度数字颜色为黑色
    colorbar.scalar_bar.unconstrained_font_size = True  # 确保字体大小可调
    colorbar.label_text_property.font_size = 10  # Increase font size
    colorbar.label_text_property.color = (0, 0, 0)  # RGB格式，黑色为 (0, 0, 0)
    colorbar.title_text_property.color = (0, 0, 0)  # Black color for the title
    colorbar.title_text_property.font_size = 14  # Increase font size


    # Apparently we need this to redraw the figure, maybe it can be done in
    # another way?
    mlab.view(azimuth=155, elevation=70, distance='auto')
    # Show the 3d plot:
    mlab.show()


from pathlib import Path
import os

# 读取数据
path_cal_res = os.path.dirname(os.path.abspath(__file__))
dens_homo, atoms = read_cube_data(path_cal_res / Path('homo_wave_function.cube'))
dens_lumo, _ = read_cube_data(path_cal_res / Path('lumo_wave_function.cube'))
print(dens_homo.min(),dens_homo.max())
print(dens_lumo.min(),dens_lumo.max())

# 检查每个元素的符号
sign = np.sign(dens_homo)

# 对波函数平方
dens_homo = np.abs(dens_homo)**2

# 如果原始波函数为负，重新调整符号
dens_homo *= sign  # 保留符号

# 检查每个元素的符号
sign = np.sign(dens_lumo)

# 对波函数平方
dens_lumo = np.abs(dens_lumo)**2

# 如果原始波函数为负，重新调整符号
dens_lumo *= sign  # 保留符号


# dens_homo = np.abs(dens_homo) ** 2
# dens_lumo = np.abs(dens_lumo) ** 2

contours = [0.001]
# contours = list(np.linspace(dens.min(),dens.max(),10))
# plot(atoms,dens_homo,list(np.linspace(dens_homo.min(),dens_homo.max(),5)))
# plot(atoms,dens_lumo,list(np.linspace(dens_lumo.min(),dens_lumo.max(),5)))
# plot(atoms,dens_homo,[0.001,dens_homo.max()/5])
# plot(atoms,dens_lumo,[0.001,dens_lumo.max()/5])

plot(atoms,dens_homo,contours)
plot(atoms,dens_lumo,contours)
```



