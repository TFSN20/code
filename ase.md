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
- 保存traj3D文件
  ```
  from ase.io import write
  write('zn_h2o_cluster.traj', zn_h2o_cluster)
  ```
  终端查看
  ```
  ase gui zn_h2o_cluster.traj
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
- 创建一个分子
