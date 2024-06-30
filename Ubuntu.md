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
