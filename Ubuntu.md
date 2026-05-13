# Ubuntu安装
- 不同版本有不同的安装文件夹名称，自行搜索，以下为Ubuntu22.04.3 LTS版本，D:/mklink/Ubuntu要先创建
  ```
  mklink /D "C:\Users\tfsn20\AppData\Local\Packages\CanonicalGroupLimited.Ubuntu22.04LTS_79rhkp1fndgsc" "D:/mklink/Ubuntu"
  ```
- 微软商店搜索下载（不要挂梯子），下载完成后打开（报错自行搜索错误代码，比如我的是0x8007019e，就管理员打开powershell，输入以下命令回车），新建用户名，回车，随便输入密码（界面不会有任何变化），回车。
  ```
  Enable-WindowsOptionalFeature -Online -FeatureName Microsoft-Windows-Subsystem-Linu
  ```
- 更改apt源（用的阿里源，清华源也可以，jammy是Ubuntu22.04的代号，每个大版本都不一样）
  ```
  sudo cp /etc/apt/sources.list /etc/apt/sources.list.backup #root备份源
  ```
  ```
  deb http://mirrors.aliyun.com/ubuntu/ jammy main restricted universe multiverse
  deb-src http://mirrors.aliyun.com/ubuntu/ jammy main restricted universe multiverse

  deb http://mirrors.aliyun.com/ubuntu/ jammy-updates main restricted universe multiverse
  deb-src http://mirrors.aliyun.com/ubuntu/ jammy-updates main restricted universe multiverse

  deb http://mirrors.aliyun.com/ubuntu/ jammy-backports main restricted universe multiverse
  deb-src http://mirrors.aliyun.com/ubuntu/ jammy-backports main restricted universe multiverse

  deb http://mirrors.aliyun.com/ubuntu/ jammy-security main restricted universe multiverse
  deb-src http://mirrors.aliyun.com/ubuntu/ jammy-security main restricted universe multiverse
  ```
  ```
  sudo apt update
  ```
  ```
  sudo apt upgrade
  ```
# 一些库的安装
- ```
  sudo apt-get install software-properties-common
  ```
## python库
### fenics
- 添加私人软件源
  ```
  sudo add-apt-repository ppa:fenics-packages/fenics
  ```
  下载
  ```
  sudo apt install fenics
  ```
# 常用命令
- 查看文件和文件夹
  ```
  ls
  ```
- 查看文件和文件夹的权限
  ```
  ls -la
  ```
- 创建文件
  ```
  touch filename.txt
  ```
- 删除文件夹下所有内容 先cd到目标文件夹的父目录
  ```
  rm -rf /目标文件夹
  ```
- 重命名
  ```
  mv name.py name2.py
  ```
- cd到一个文件夹a下，将所有cal开头的文件移动到文件夹a下的b文件夹
  ```
  mv cal*.* b文件夹/
  ```
- 创建一个文件
  ```
  touch xx.txt
  ```
- 编辑器编辑文本,光标移动，ctrl+x退出，y保存，回车确认文件名保存
  ```
  nano xxx.txt
  ```
- 代理全局
  ```
  export host_ip=$(ip route | grep default | awk '{print $3}')
  export http_proxy="http://${host_ip}:10808"
  export https_proxy="http://${host_ip}:10808"
  export all_proxy="socks5:/${host_ip}:10808"
  export no_proxy="localhost,127.0.0.1,*.local"
  ```
- 禁止dns自动生成 
  ```
  nano /etc/wsl.conf
  ```
  ```
  [network]
  generateResolvConf = false
  [boot]
  systemd=true
  ```
  dns修改:
  ```
  nano /etc/resolv.conf
  ```
  ```
  nameserver 223.5.5.5
  nameserver 180.76.76.76
  ```
- 取消代理
  ```
  unset http_proxy
  unset https_proxy
  unset all_proxy
  ```
  
  
# 常见问题
- -rw-r--r-- 1 tfsn20 tfsn20    0 Dec 16 10:32 my_script.py 能执行python3 my_script.py
  Python脚本不需要执行权限才能被运行。关键是Python解释器需要执行权限,而不是脚本本身。
- sudo密码忘记
方法一：从 Windows PowerShell 重置（推荐，最简单）
步骤 1：确认发行版名称
在 Windows PowerShell 或 CMD 中运行：
wsl -l -v
输出示例：
  NAME      STATE           VERSION
* Ubuntu    Running         2
记住你的发行版名称（如 Ubuntu、Ubuntu-20.04 等）。
步骤 2：以 root 用户进入 WSL2
在 Windows PowerShell 中运行：
wsl -d Ubuntu -u root
将 Ubuntu 替换为你的实际发行版名称
此时你会直接进入 WSL2 的 root shell，无需密码。
步骤 3：重置用户密码
在 WSL2 内运行：
```
# 查看当前用户列表（确认你的用户名）
ls /home
# 重置密码（将 username 替换为你的实际用户名）
passwd username
```
系统会提示你输入新密码（输入时不显示字符），确认两次即可。
- 查看wsl1或wsl2`wsl -l -v`
## 安装
### 前置工作
- 控制面板里启用linux和子系统
- powershell管理员运行
  ```
  dism.exe /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux /all /norestart
  dism.exe /online /enable-feature /featurename:VirtualMachinePlatform /all /norestart
  ```
  重启
- 查看并复制一个发行版名称(NAME)
  ```
  wsl --list --online
  # 或
  wsl -l -o
  ```
- 确保C盘至少剩余10GB+，因为默认安装到C:\Users\<用户名>\AppData\Local\Packages\下，安装完成后重启电脑
  ```
  wsl --install -d <发行版名称>
  ```
- 命令安装默认wsl2，如果想转换到wsl1
  ```
  wsl --set-version <发行版名称> 1
  ```
  或在安装前指定全局wsl1
  ```
  wsl --set-default-version 1
  ```
- 查看已安装的发行版名称、状态、wsl版本
  ```
  wsl --list --verbose
  # 或
  wsl -l -v
  ```
- 停止运行中的 WSL
  ```
  wsl --shutdown
  ```
- 导出/备份指定发行版
  ```
  wsl --export <发行版名称> C:\backup\ubuntu2204-backup.tar
  ```
- 从文件导入安装/改名
  ```
  # 语法：wsl --import <新名字> <安装路径> <备份文件> --version 1
  wsl --import Ubuntu-22.04-WSL1 D:\WSL\Ubuntu-WSL1 C:\backup\ubuntu2204-backup.tar --version 1
  ```
- 通过 --import 导入的系统，默认登录用户会变成 root，若想自动以普通用户登录
  ```
  # 查看用户名
  ls /home
  # 编辑配置文件
  nano /etc/wsl.conf
  ```
  输入
  ```
  [user]
  default=<你的用户名>
  ```
- 卸载某个发行版（会删除所有数据）
  ```
  wsl --unregister <发行版名称>
  ```
- 设置默认启动的发行版
  ```
  wsl --set-default <发行版名称>
  ```
- 以特定用户身份启动并进入指定目录（~是用户目录）
  ```
  wsl -d <发行版名称> -u <用户名> --cd ~
  ```
## 源
- 官方源
  ```
  # 备份
  sudo cp /etc/apt/sources.list /etc/apt/sources.list.bak.$(date +%F)
  # 打开修改
  sudo nano /etc/apt/sources.list
  # 清空原有，写入下面
  deb https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ jammy main restricted universe multiverse
  deb https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ jammy-updates main restricted universe multiverse
  deb https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ jammy-backports main restricted universe multiverse
  deb https://mirrors.tuna.tsinghua.edu.cn/ubuntu/ jammy-security main restricted universe multiverse
  # 更新索引
  sudo apt update
  ```
## 安装gaussian
- 下载bash
  ```
  sudo apt update
  sudo apt install csh tcsh
  ```
- 选择下载，我选择了g16-C01-AVX优化版本
- 解压并赋予权限，添加scr临时文件目录，初始化
  ```
  cd
  mkdir gaussian
  cd gaussian
  tar xvf G16-C01-AVX.tbJ
  cd g16
  chmod -R 755 *
  mkdir -p ~/gaussian/scr
  ./bsd/install
  ```
- 设置变量
  ```
  nano $HOME/.bashrc
  ```
  写入以下
  ```
  export g16root=$HOME/gaussian
  export GAUSS_EXEDIR=$g16root/g16
  export GAUSS_SCRDIR=$g16root/scr
  source $g16root/g16/bsd/g16.profile
  ```
  加载bashrc
  ```
  source $HOME/.bahsrc
  ```
  
## uv
### win
- powershell临时变量
  ```
  $env:MY_VAR = "Value"
  ```
- 安装uv到指定目录（无需python/pip环境，无需管理员权限，win上还会自动将"D:\Tools\uv"添加到用户Path变量）
  ```
  $env:UV_INSTALL_DIR = "D:\Tools\uv"
  powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
  ```
- 更新uv
  ```
  uv self update
  ```
- uv版本
  ```
  uv -v
  ```
- uv位置（powershell会把where识别为Where-Object，因此需要全名）
  ```
  where.exe uv
  ```
- python下载位置，包下载位置，包下载源查看和更改（通过环境变量）
  ```
  # 默认C:\Users\用户名\AppData\Roaming\uv\python\
  UV_PYTHON_INSTALL_DIR

  # 默认C:\Users\用户名\AppData\Local\uv\cache
  UV_CACHE_DIR

  UV_INDEX_URL
  https://pypi.tuna.tsinghua.edu.cn/simple
  ```
### linux
- 安装uv
  ```
  # 1. 确保系统基础工具已更新（可选，但推荐）
  sudo apt update && sudo apt install -y curl

  # 2. 执行官方安装脚本
  # 原理：下载 install.sh 脚本 -> 脚本检测架构 -> 下载对应 uv 二进制文件到 ~/.local/bin -> 自动尝试修改 shell 配置文件
  curl -LsSf https://astral.sh/uv/install.sh | sh

  # 3. 验证安装
  uv --version
- 更改python下载源和pip源
  ```
  # 1. 打开配置文件 (根据你的 shell 选择 .bashrc 或 .zshrc)
  nano ~/.bashrc

  # 2. 在文件末尾添加以下内容 (以npmmirror为例)
  # 原理：每次启动终端时自动设置该环境变量
  export UV_PYTHON_INSTALL_MIRROR="https://registry.npmmirror.com/-/binary/python-build-standalone/"
  export UV_INDEX_URL="https://pypi.tuna.tsinghua.edu.cn/simple"
  # 3. 保存退出 (Ctrl+O, Enter, Ctrl+X)，然后刷新配置
  source ~/.bashrc
  ```
### 使用
- 查看可安装和已安装
  ```
  uv python list
  ```
- 下载python
  ```
  uv python install 3.12
  uv python install 3.10
  ```
- 指定python版本创建项目（不指定时默认使用最新python），下载项目需要的库
  ```
  # 初始化 uv 项目
  uv init gaussian_project_python3.10 --python 3.10
  cd gaussian_project
  # 添加计算化学常用包
  uv add numpy pandas scipy ase rdkit cclib
  ```
- 创建项目文件夹后进入
  ```
  uv init
  uv python pin 3.12
  ```
- 修改项目使用的python版本（会重新选择和下载项目使用的库）
  ```
  uv python pin 3.8
  uv sync --reinstall
  ```
- 项目库信息
  ```
  uv pip list
  # 或
  uv pip freeze
  ```
  库间依赖
  ```
  uv pip tree
  ```
## 其他工具
### pdf可搜索（中文为主，混合英文专有名词）
- 使用ocrmypdf
- 文件夹内批量处理
  ```
  mkdir output

  for f in *.pdf; do
    ocrmypdf \
    -l chi_sim+eng \
    --image-dpi 300 \
    --rotate-pages \
    --deskew \
    --clean \
    --skip-text \
    --clean-final \
    "$f" \
    "output/$f"
  done
  ```
