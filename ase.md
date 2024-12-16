- 单独ase可在windows上运行，gpaw需要在linux系统上运行（如Ubuntu）。
- Ubuntu可视化需要ssh远程连接桌面，在Ubuntu上保存traj文件后，可在windows平台上使用```ase gui myatoms.traj```运行查看3D文件。
- 安装gpaw数据库
  ```
  gpaw install-data /home/tfsn20/paw_datasets
  ```
  等待下载，按下y
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
  设置虚拟晶胞
  ```
  water_cluster.set_cell([20, 20, 20])  # 足够大的立方晶胞
  water_cluster.center()  # 将体系置于晶胞中心
  water_cluster.set_pbc(False)  # 关闭周期性边界条件
  ```
  保存未优化的traj3D文件
  ```
  from ase.io import write
  write('water_cluster_no_optim.traj', water_cluster)
  ```
  设置 GPAW 计算器，保存计算器数据
  ```
  from gpaw import GPAW, PW
  calc = GPAW(mode=PW(300), xc='PBE', txt='water_cluster_output.txt')
  water_cluster.calc = calc
  ```
  优化结构，保存优化中的动态traj3D文件
  ```
  from ase.optimize import BFGS
  dyn = BFGS(water_cluster, trajectory='water_cluster_optiming.traj')
  dyn.run(fmax=0.1)
  ```
  获取能量
  ```
  water_cluster.get_potential_energy()
  ```
  继续计算
  ```
  from ase.io import Trajectory
  from ase.optimize import BFGS
  from gpaw import GPAW, PW
  
  # 设置 GPAW 计算器
  calc = GPAW(mode=PW(300), xc='PBE', txt='zn2_plus_ion_output_keeping.txt')
  # 加载轨迹文件并获取最后一步的结构
  traj = Trajectory('zn2_plus_ion_optiming.traj', 'r')  # 以只读模式加载轨迹
  system = traj[-1]  # 获取轨迹文件中的最后一步结构
  
  # 重新设置计算器
  system.calc=calc
  
  # 继续优化并保存到新的轨迹文件
  dyn = BFGS(system, trajectory='zn2_plus_ion_optiming_keeping.traj')
  dyn.run(fmax=0.05)
  ```
  终端查看
  ```
  ase gui zn_h2o_cluster_no_optim.traj
  ```
  代码里查看
  ```
  from ase.visualize import view
  view(atoms)
  ```
- 创建一个离子
  ```
  from ase import Atoms
  zn2_plus = Atoms('Zn', charges=[2.0])
  ```
## 常见问题
